"""Count-matrix tools: fetch a processed count matrix from GEO instead of FASTQs.

THE CONTRACT (same as tools.py):
  Every tool returns a structured, JSON-serializable summary dict. URLs, sample
  metadata and the parsed matrix are saved to disk and read back by later tools —
  the LLM refers to files by name and to samples by GSM ID, never relays URLs.

The LLM makes the judgement calls (which file holds raw counts, which column is the
gene ID, which column is which GSM, which characteristic is the condition). The tools
parse, validate and refuse anything inconsistent.

Supported (step 1): single-matrix series supplementary files (csv/tsv/txt, optionally
gzipped) and NCBI-generated raw counts. Per-sample files, tar archives and xlsx are
listed but not yet parsed.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import re
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

import pandas as pd

from core import config
from core.session import SESSION

Summary = dict[str, Any]

_GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
_GEO_DOWNLOAD = "https://www.ncbi.nlm.nih.gov/geo/download/"
_GEO_ACC = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
_GENE_INFO = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia"
_GENE_INFO_FILES = {"Homo sapiens": "Homo_sapiens.gene_info.gz", "Mus musculus": "Mus_musculus.gene_info.gz"}

_MAX_BYTES = 500 * 1024**2
_MAX_RETRIES = 3
_TABULAR = re.compile(r"\.(txt|tsv|csv|tab)(\.gz)?$", re.IGNORECASE)
_ANNOTATION_NAMES = {
    "chr", "start", "end", "strand", "length", "gene_name", "genename", "symbol", "gene_symbol",
    "description", "biotype", "gene_biotype", "gene_type", "entrezid", "entrez_id",
}
_MAX_SAMPLES_SHOWN = 100


class BlockedByNCBI(RuntimeError):
    """NCBI returned an HTML page (captcha / rate limit) where data was expected."""


# --- HTTP ---------------------------------------------------------------------


def _fetch(url: str, *, timeout: int = 60) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "agentic-rnaseq-workflow/0.1"})
    for attempt in range(_MAX_RETRIES):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read(_MAX_BYTES + 1)
            if len(data) > _MAX_BYTES:
                raise ValueError(f"File exceeds {_MAX_BYTES // 1024**2} MB limit: {url}")
            return data
        except urllib.error.HTTPError as e:
            if e.code == 429 and attempt < _MAX_RETRIES - 1:
                time.sleep(2 ** (attempt + 1))
                continue
            raise
    return b""


def _fetch_data(url: str) -> bytes:
    """Fetch a data file, refusing HTML (NCBI serves a captcha page when rate limiting)."""
    data = _fetch(url)
    head = data[:200].lstrip().lower()
    if head.startswith((b"<!doctype html", b"<html")):
        raise BlockedByNCBI(
            "NCBI returned a web page instead of the file (likely rate limiting / captcha). "
            "Wait a few minutes and retry."
        )
    return data


# --- GEO parsing --------------------------------------------------------------


def _series_dir(gse: str) -> str:
    digits = gse[3:]
    stub = f"GSE{digits[:-3]}nnn" if len(digits) > 3 else "GSEnnn"
    return f"{_GEO_FTP}/{stub}/{gse}"


def _parse_ftp_index(html: str) -> list[dict]:
    """Parse an NCBI FTP-over-HTTPS directory listing into [{name, size}]."""
    files = []
    for m in re.finditer(r'<a href="([^"?/][^"]*)">[^<]*</a>\s+\S+\s+\S+\s+(\S+)', html):
        files.append({"name": m.group(1), "size": m.group(2)})
    return files


def _parse_filelist(text: str) -> list[dict]:
    """Parse a GEO suppl filelist.txt (contents of _RAW.tar)."""
    out = []
    for line in text.splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) >= 5 and parts[0] == "File":
            out.append({"name": parts[1], "bytes": int(parts[3]) if parts[3].isdigit() else None})
    return out


def _parse_samples_soft(text: str) -> dict[str, dict]:
    """Parse GEO SOFT (targ=gsm, view=brief) into {GSM: {title, organism, strategy, characteristics}}."""
    samples: dict[str, dict] = {}
    current: dict | None = None
    for line in text.splitlines():
        if line.startswith("^SAMPLE"):
            gsm = line.split("=", 1)[1].strip()
            current = {"title": "", "organism": "", "library_strategy": "", "characteristics": {}}
            samples[gsm] = current
        elif current is None or " = " not in line:
            continue
        else:
            key, value = line.split(" = ", 1)
            value = value.strip()
            if key == "!Sample_title":
                current["title"] = value
            elif key == "!Sample_organism_ch1":
                current["organism"] = value
            elif key == "!Sample_library_strategy":
                current["library_strategy"] = value
            elif key == "!Sample_characteristics_ch1" and ":" in value:
                k, v = value.split(":", 1)
                current["characteristics"][_snake(k)] = v.strip()
    return samples


def _parse_ncbi_counts_links(html: str, gse: str) -> list[dict]:
    """Find NCBI-generated raw count files on the GEO download page."""
    out = []
    for m in re.finditer(r'href="(/geo/download/\?[^"]*type=rnaseq_counts[^"]*)"', html):
        href = m.group(1).replace("&amp;", "&")
        fm = re.search(r"file=([^&]+)", href)
        if fm and "_raw_counts_" in fm.group(1) and gse in fm.group(1):
            out.append({"name": fm.group(1), "url": "https://www.ncbi.nlm.nih.gov" + href})
    return out


def _snake(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", text.strip().lower()).strip("_") or "field"


def _sample_name(title: str, gsm: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", title).strip("_") or gsm


def _normalise(text: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(text).lower())


# --- Table parsing ------------------------------------------------------------


def _read_table(raw: bytes) -> tuple[pd.DataFrame, dict]:
    """Parse a delimited count table. Returns (all-string DataFrame, format info)."""
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    if raw[:4] == b"PK\x03\x04":
        raise ValueError("Excel/zip files are not supported yet.")
    text = raw.decode("utf-8-sig", errors="replace")
    lines = text.splitlines()
    skip = 0
    while skip < len(lines) and lines[skip].startswith("#"):  # e.g. featureCounts header
        skip += 1
    body = [ln for ln in lines[skip:skip + 3] if ln.strip()]
    if len(body) < 2:
        raise ValueError("File has fewer than two non-comment lines.")
    if "\t" in body[0]:
        sep, sep_name = "\t", "tab"
    elif "," in body[0]:
        sep, sep_name = ",", "comma"
    else:
        sep, sep_name = r"\s+", "whitespace"

    df = pd.read_csv(io.StringIO(text), sep=sep, skiprows=skip, dtype=str, engine="python")
    row_names = not isinstance(df.index, pd.RangeIndex)  # R write.table: header lacks the index name
    if row_names:
        df = df.reset_index().rename(columns={"index": "row_names"})
    df.columns = [str(c).strip().strip('"') for c in df.columns]
    return df, {"separator": sep_name, "comment_lines_skipped": skip, "unnamed_first_column": row_names}


def _numeric_profile(col: pd.Series) -> dict:
    values = pd.to_numeric(col, errors="coerce")
    numeric = bool(values.notna().mean() > 0.99)
    out: dict[str, Any] = {"numeric": numeric}
    if numeric:
        v = values.dropna()
        out["integer"] = bool((v == v.round()).all())
    return out


def _value_type(counts: pd.DataFrame) -> tuple[str, float]:
    values = counts.to_numpy(dtype=float)
    pct_non_int = float((values != values.round()).mean())
    col_sums = counts.sum(axis=0)
    if pct_non_int < 0.005:
        return "raw_integer_counts", pct_non_int
    if values.max() < 30:
        return "log_scale", pct_non_int
    if ((col_sums - 1e6).abs() / 1e6 < 0.01).all():
        return "tpm_like", pct_non_int
    return "fractional (expected counts, or FPKM/normalised)", pct_non_int


def _gene_id_type(ids: pd.Series) -> str:
    sample = ids.dropna().astype(str).head(500)
    if sample.str.fullmatch(r"\d+").mean() > 0.9:
        return "entrez"
    if sample.str.match(r"ENS[A-Z]*G\d+").mean() > 0.9:
        return "ensembl"
    return "symbol"


def _symbol_maps(organism: str) -> tuple[dict[str, str], dict[str, str]] | None:
    """(entrez -> symbol, ensembl -> symbol) from NCBI gene_info, cached in REFERENCE_DIR."""
    fname = _GENE_INFO_FILES.get(organism)
    if not fname:
        return None
    cache = config.REFERENCE_DIR / fname
    if not cache.is_file():
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_bytes(_fetch_data(f"{_GENE_INFO}/{fname}"))
    info = pd.read_csv(cache, sep="\t", usecols=["GeneID", "Symbol", "dbXrefs"], dtype=str)
    entrez = dict(zip(info["GeneID"], info["Symbol"]))
    ensembl = {}
    for symbol, xrefs in zip(info["Symbol"], info["dbXrefs"]):
        for m in re.finditer(r"Ensembl:(ENS[A-Z]*G\d+)", xrefs or ""):
            ensembl[m.group(1)] = symbol
    return entrez, ensembl


def _load_sources() -> dict | None:
    path = SESSION.require_paths().geo_sources
    return json.loads(path.read_text()) if path.is_file() else None


def _get_file(sources: dict, filename: str) -> tuple[bytes, dict]:
    """Return (bytes, file entry) for a listed file, downloading into geo/ on first use."""
    entry = next((f for f in sources["files"] if f["name"] == filename), None)
    if entry is None:
        raise KeyError(filename)
    local = SESSION.require_paths().geo_dir / filename
    if not local.is_file():
        local.parent.mkdir(parents=True, exist_ok=True)
        local.write_bytes(_fetch_data(entry["url"]))
    return local.read_bytes(), entry


# --- Agent tools --------------------------------------------------------------


def list_geo_count_sources(accession: str) -> Summary:
    """List candidate count files for a GEO series, plus per-sample metadata.

    Saves URLs and sample metadata to geo_sources.json for the other tools.
    """
    gse = accession.strip().upper()
    if not re.fullmatch(r"GSE\d+", gse):
        return {"error": "bad_accession", "message": "Only GSE accessions are supported for count matrices."}

    series = _series_dir(gse)
    try:
        soft = _fetch(f"{_GEO_ACC}?acc={gse}&targ=gsm&form=text&view=brief").decode("utf-8", "replace")
        samples = _parse_samples_soft(soft)
        if not samples:
            return {"error": "not_found", "message": f"No GEO samples found for {gse}. Check the accession."}
        try:
            listing = _parse_ftp_index(_fetch(f"{series}/suppl/").decode("utf-8", "replace"))
        except urllib.error.HTTPError as e:
            if e.code != 404:
                raise
            listing = []  # series has no supplementary files
        download_page = _fetch(f"{_GEO_DOWNLOAD}?acc={gse}").decode("utf-8", "replace")
    except (urllib.error.URLError, BlockedByNCBI, ValueError) as exc:
        return {"error": "fetch_failed", "message": f"Could not query GEO: {exc}"}

    files = []
    for f in listing:
        if f["name"] == "filelist.txt":
            continue
        entry = {"name": f["name"], "size": f["size"], "source": "author", "url": f"{series}/suppl/{f['name']}"}
        if _TABULAR.search(f["name"]):
            entry["kind"], entry["supported"] = "table", True
        elif f["name"].lower().endswith(".tar"):
            entry["kind"], entry["supported"] = "archive (per-sample files)", False
        elif re.search(r"\.xlsx?(\.gz)?$", f["name"], re.IGNORECASE):
            entry["kind"], entry["supported"] = "excel", False
        else:
            entry["kind"], entry["supported"] = "other", False
        files.append(entry)

    if any(f["kind"].startswith("archive") for f in files):
        try:
            contents = _parse_filelist(_fetch(f"{series}/suppl/filelist.txt").decode("utf-8", "replace"))
            for f in files:
                if f["kind"].startswith("archive"):
                    f["contents"] = contents
        except urllib.error.URLError:
            pass

    for link in _parse_ncbi_counts_links(download_page, gse):
        files.append({**link, "size": None, "source": "ncbi_generated", "kind": "table", "supported": True})

    organisms = sorted({s["organism"] for s in samples.values() if s["organism"]})
    paths = SESSION.require_paths()
    paths.geo_sources.write_text(json.dumps({"accession": gse, "files": files, "samples": samples}, indent=2) + "\n")

    shown = dict(list(samples.items())[:_MAX_SAMPLES_SHOWN])
    return {
        "accession": gse,
        "organisms": organisms,
        "n_samples": len(samples),
        "samples": {
            gsm: {k: v for k, v in s.items() if k != "organism" or len(organisms) > 1}
            for gsm, s in shown.items()
        },
        "samples_truncated": len(samples) > _MAX_SAMPLES_SHOWN,
        "files": [
            {k: v for k, v in f.items() if k != "url" and (k != "contents" or len(v) <= 50)}
            for f in files
        ],
        "ncbi_generated_available": any(f["source"] == "ncbi_generated" for f in files),
    }


def preview_geo_file(filename: str) -> Summary:
    """Download (cached) and preview a listed count file: format, columns, first rows."""
    sources = _load_sources()
    if sources is None:
        return {"error": "no_sources", "message": "Call list_geo_count_sources first."}
    try:
        raw, entry = _get_file(sources, filename)
    except KeyError:
        return {"error": "unknown_file", "message": f"'{filename}' is not in the listed files."}
    except (urllib.error.URLError, BlockedByNCBI, ValueError) as exc:
        return {"error": "fetch_failed", "message": str(exc)}
    if not entry["supported"]:
        return {"error": "unsupported", "message": f"'{filename}' ({entry['kind']}) is not supported yet."}
    try:
        df, fmt = _read_table(raw)
    except (ValueError, pd.errors.ParserError) as exc:
        return {"error": "parse_failed", "message": str(exc)}

    titles = {gsm: s["title"] for gsm, s in sources["samples"].items()}
    by_norm = {_normalise(gsm): gsm for gsm in titles} | {_normalise(t): gsm for gsm, t in titles.items() if t}
    columns = []
    for col in df.columns:
        info: dict[str, Any] = {"name": col, **_numeric_profile(df[col])}
        if _normalise(col) in by_norm:
            info["matches_gsm"] = by_norm[_normalise(col)]
        if col.lower() in _ANNOTATION_NAMES:
            info["looks_like_annotation"] = True
        columns.append(info)

    head = df.head(5)
    return {
        "filename": filename,
        "source": entry["source"],
        "format": fmt,
        "n_rows": int(df.shape[0]),
        "n_columns": int(df.shape[1]),
        "columns": columns[:80],
        "columns_truncated": len(columns) > 80,
        "first_rows": [
            {c: str(v)[:40] for c, v in row.items()} for row in head.iloc[:, :10].to_dict("records")
        ],
    }


def fetch_geo_counts(filename: str, gene_id_column: str, sample_map: dict[str, str]) -> Summary:
    """Parse a listed count file into counts.tsv, mapping file columns to GSM IDs.

    sample_map is {file column: GSM ID}. Every value must be a GSM of this series, used
    once. Samples are renamed to their (sanitised) GEO titles. Duplicate gene IDs are
    summed. Gene symbols are added from NCBI gene_info for human/mouse Entrez/Ensembl IDs.
    """
    sources = _load_sources()
    if sources is None:
        return {"error": "no_sources", "message": "Call list_geo_count_sources first."}
    try:
        raw, entry = _get_file(sources, filename)
    except KeyError:
        return {"error": "unknown_file", "message": f"'{filename}' is not in the listed files."}
    except (urllib.error.URLError, BlockedByNCBI, ValueError) as exc:
        return {"error": "fetch_failed", "message": str(exc)}
    if not entry["supported"]:
        return {"error": "unsupported", "message": f"'{filename}' ({entry['kind']}) is not supported yet."}
    try:
        df, _ = _read_table(raw)
    except (ValueError, pd.errors.ParserError) as exc:
        return {"error": "parse_failed", "message": str(exc)}

    samples = sources["samples"]
    if gene_id_column not in df.columns:
        return {"error": "bad_gene_id_column", "message": f"'{gene_id_column}' is not a column. Columns: {list(df.columns)[:20]}"}
    if not sample_map:
        return {"error": "empty_sample_map", "message": "sample_map must map at least one column to a GSM."}
    missing_cols = [c for c in sample_map if c not in df.columns]
    if missing_cols:
        return {"error": "unknown_columns", "message": f"Not columns in the file: {missing_cols[:10]}"}
    bad_gsms = [g for g in sample_map.values() if g not in samples]
    if bad_gsms:
        return {"error": "unknown_gsm", "message": f"Not samples of {sources['accession']}: {bad_gsms[:10]}"}
    dupes = sorted({g for g in sample_map.values() if list(sample_map.values()).count(g) > 1})
    if dupes:
        return {"error": "duplicate_gsm", "message": f"Each GSM must map to exactly one column; repeated: {dupes}"}

    counts = df[list(sample_map)].apply(pd.to_numeric, errors="coerce")
    bad = [c for c in counts.columns if counts[c].isna().any()]
    if bad:
        return {"error": "non_numeric", "message": f"Columns with non-numeric or missing values: {bad[:10]}"}
    if (counts < 0).any().any():
        return {"error": "negative_values", "message": "Negative values found — this is not a count matrix (log ratios?)."}

    gene_ids = df[gene_id_column].astype(str).str.strip()
    keep = gene_ids.ne("") & gene_ids.ne("nan")
    id_type = _gene_id_type(gene_ids[keep])
    if id_type == "ensembl":
        gene_ids = gene_ids.str.replace(r"\.\d+$", "", regex=True)
    counts = counts[keep].set_axis(gene_ids[keep])
    n_before = len(counts)
    counts = counts.groupby(level=0, sort=False).sum()
    n_duplicates_summed = n_before - len(counts)

    names = {col: _sample_name(samples[gsm]["title"], gsm) for col, gsm in sample_map.items()}
    if len(set(names.values())) < len(names):
        names = {col: f"{names[col]}_{gsm}" for col, gsm in sample_map.items()}
    counts = counts.rename(columns=names)
    counts.index.name = "gene_id"

    organisms = sorted({samples[g]["organism"] for g in sample_map.values()})
    gene_names = None
    symbol_note = None
    if id_type == "symbol":
        gene_names = counts.index.to_series()
    elif len(organisms) == 1:
        try:
            maps = _symbol_maps(organisms[0])
        except (urllib.error.URLError, BlockedByNCBI, ValueError) as exc:
            maps, symbol_note = None, f"gene_info download failed: {exc}"
        if maps:
            lookup = maps[0] if id_type == "entrez" else maps[1]
            gene_names = counts.index.to_series().map(lookup)
        elif symbol_note is None:
            symbol_note = f"No gene_info mapping for organism {organisms[0]!r}."

    value_type, pct_non_int = _value_type(counts)
    out = counts.reset_index()
    if gene_names is not None:
        out.insert(1, "gene_name", gene_names.fillna(counts.index.to_series()).to_numpy())

    paths = SESSION.require_paths()
    out.to_csv(paths.counts_matrix, sep="\t", index=False)
    dropped = [c for c in df.columns if c not in sample_map and c != gene_id_column]
    unmapped = [g for g in samples if g not in sample_map.values()]
    sample_table = [
        {"sample": names[col], "gsm": gsm, "title": samples[gsm]["title"], "file_column": col}
        for col, gsm in sample_map.items()
    ]
    metadata = {
        "accession": sources["accession"],
        "source": entry["source"],
        "filename": filename,
        "url": entry["url"],
        "md5": hashlib.md5(raw).hexdigest(),
        "gene_id_column": gene_id_column,
        "gene_id_type": id_type,
        "value_type": value_type,
        "pct_non_integer": round(pct_non_int, 4),
        "n_genes": int(counts.shape[0]),
        "n_duplicates_summed": n_duplicates_summed,
        "samples": sample_table,
        "dropped_columns": dropped,
        "unmapped_gsms": unmapped,
    }
    paths.counts_metadata.write_text(json.dumps(metadata, indent=2) + "\n")

    return {
        "counts_path": str(paths.counts_matrix),
        "n_genes": metadata["n_genes"],
        "n_samples": len(sample_table),
        "samples": {s["sample"]: s["gsm"] for s in sample_table},
        "gene_id_type": id_type,
        "gene_names_added": gene_names is not None,
        "symbol_note": symbol_note,
        "value_type": value_type,
        "pct_non_integer": metadata["pct_non_integer"],
        "n_duplicates_summed": n_duplicates_summed,
        "dropped_columns": dropped[:20],
        "unmapped_gsms": [{"gsm": g, "title": samples[g]["title"]} for g in unmapped[:20]],
    }


def save_geo_design(condition_field: str | None = None, conditions: dict[str, str] | None = None) -> Summary:
    """Write design.csv for the fetched samples from GEO characteristics.

    Give exactly one of: condition_field (a characteristics key, e.g. "treatment") or
    conditions ({GSM: label}, when the condition is only in the titles). Every other
    characteristic that varies across samples is kept as a covariate column; constant
    ones (e.g. time point) are kept as columns too, as facts for the report.
    """
    paths = SESSION.require_paths()
    sources = _load_sources()
    if sources is None or not paths.counts_metadata.is_file():
        return {"error": "no_counts", "message": "Call fetch_geo_counts first."}
    if (condition_field is None) == (conditions is None):
        return {"error": "bad_arguments", "message": "Give exactly one of condition_field or conditions."}

    sample_table = json.loads(paths.counts_metadata.read_text())["samples"]
    chars = {s["gsm"]: sources["samples"][s["gsm"]]["characteristics"] for s in sample_table}

    if condition_field is not None:
        missing = [g for g, c in chars.items() if condition_field not in c]
        if missing:
            keys = sorted({k for c in chars.values() for k in c})
            return {"error": "bad_field", "message": f"'{condition_field}' missing for {missing[:5]}. Keys: {keys}"}
        labels = {g: c[condition_field] for g, c in chars.items()}
    else:
        missing = [g for g in chars if g not in conditions]
        if missing:
            return {"error": "missing_conditions", "message": f"No condition given for: {missing[:10]}"}
        labels = {g: str(conditions[g]) for g in chars}
    if len(set(labels.values())) < 2:
        return {"error": "single_condition", "message": "All samples have the same condition — nothing to compare."}

    keys = sorted({k for c in chars.values() for k in c} - {condition_field})
    varying = [k for k in keys if len({c.get(k, "") for c in chars.values()}) > 1]
    # Constant characteristics (e.g. time_point, cell type) are kept too: they're facts the
    # report needs, and constant columns are never used as covariates or plotted.
    constant = [k for k in keys if k not in varying]
    rows = [
        {"sample": s["sample"], "condition": labels[s["gsm"]], "gsm": s["gsm"], "title": s["title"],
         **{k: chars[s["gsm"]].get(k, "") for k in varying + constant}}
        for s in sample_table
    ]
    pd.DataFrame(rows).to_csv(paths.design, index=False)

    counts = pd.Series(list(labels.values())).value_counts().to_dict()
    return {
        "design_path": str(paths.design),
        "conditions": counts,
        "covariates": varying,
        "constant_characteristics": {k: next(iter(chars.values())).get(k, "") for k in constant},
        "low_replication": [c for c, n in counts.items() if n < 3],
    }

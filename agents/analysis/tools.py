"""Downstream analysis tools for bulk RNA-seq.

THE CONTRACT:
  Every tool is a deterministic Python function that does real work AND returns a
  structured, JSON-serializable summary dict. That dict is the *only* thing the LLM
  sees — it decides the next step entirely from it. A tool that returns None is a bug.

Tools operate on DataFrames held by the running session (see core/session.py). Each
tool validates the state it expects and returns a clear error summary the agent can
react to.
"""

from __future__ import annotations

import urllib.request
import urllib.error
import json as _json
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from core.session import SESSION

Summary = dict[str, Any]

_HEATMAP_PER_DIRECTION = 25
_DE_PADJ = 0.05  # significance threshold used for DE summaries, tables and figures
_TOP_TABLE_ROWS = 10
_ENRICHR_LIBRARIES = {
    "human": ["GO_Biological_Process_2023", "KEGG_2021_Human"],
    "mouse": ["GO_Biological_Process_2023", "KEGG_2019_Mouse"],
    "yeast": ["GO_Biological_Process_2023", "KEGG_2019"],
}

_QUANTILES = {"min": 0.0, "p25": 0.25, "median": 0.5, "p75": 0.75, "p95": 0.95, "max": 1.0}


def _quantile_summary(values) -> dict[str, float]:
    arr = np.asarray(values, dtype=float)
    return {name: round(float(np.quantile(arr, q)), 4) for name, q in _QUANTILES.items()}


_NCBI_TIMEOUT = 15


def _ncbi_get(url: str) -> str:
    """Fetch a URL from NCBI with a short timeout."""
    req = urllib.request.Request(url, headers={"User-Agent": "agentic-rnaseq-workflow/0.1"})
    with urllib.request.urlopen(req, timeout=_NCBI_TIMEOUT) as resp:
        return resp.read().decode("utf-8")


def fetch_geo_metadata(accession: str) -> Summary:
    """Fetch GEO series metadata via NCBI E-utilities.

    Returns title, summary, organism, sample descriptions, and PubMed IDs.
    The agent should call this early to understand the experimental context.
    """
    accession = accession.strip().upper()
    if not accession.startswith("GSE"):
        return {"error": "bad_accession", "message": "Only GSE accessions are supported."}

    try:
        search_url = (
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            f"?db=gds&term={accession}[ACCN]&retmode=json"
        )
        search_data = _json.loads(_ncbi_get(search_url))
        id_list = search_data.get("esearchresult", {}).get("idlist", [])
        if not id_list:
            return {"error": "not_found", "message": f"No GEO entry found for {accession}."}

        geo_id = id_list[0]
        summary_url = (
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            f"?db=gds&id={geo_id}&retmode=json"
        )
        summary_data = _json.loads(_ncbi_get(summary_url))
        record = summary_data.get("result", {}).get(str(geo_id), {})

        samples = []
        for s in record.get("samples", []):
            samples.append({"accession": s.get("accession", ""), "title": s.get("title", "")})

        return {
            "accession": accession,
            "title": record.get("title", ""),
            "summary": record.get("summary", ""),
            "organism": record.get("taxon", ""),
            "n_samples": record.get("n_samples", len(samples)),
            "samples": samples,
            "pubmed_ids": record.get("pubmedids", []),
            "platform": record.get("gpl", ""),
        }
    except (urllib.error.URLError, OSError, KeyError, _json.JSONDecodeError) as exc:
        return {"error": "fetch_failed", "message": f"Could not fetch GEO metadata: {exc}"}


def fetch_abstract(pmid: str) -> Summary:
    """Fetch a PubMed abstract via NCBI E-utilities.

    Returns title, authors, journal, year, and abstract text. The agent can use
    this to cite the original paper and contextualise findings.
    """
    pmid = str(pmid).strip()
    if not pmid.isdigit():
        return {"error": "bad_pmid", "message": "PMID must be a numeric string."}

    try:
        url = (
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            f"?db=pubmed&id={pmid}&retmode=xml"
        )
        xml_text = _ncbi_get(url)
        root = ET.fromstring(xml_text)
        article = root.find(".//PubmedArticle/MedlineCitation/Article")
        if article is None:
            return {"error": "not_found", "message": f"No PubMed article found for PMID {pmid}."}

        title_el = article.find("ArticleTitle")
        title = title_el.text if title_el is not None and title_el.text else ""

        abstract_parts = []
        abstract_el = article.find("Abstract")
        if abstract_el is not None:
            for part in abstract_el.findall("AbstractText"):
                label = part.get("Label", "")
                text = "".join(part.itertext()).strip()
                if label:
                    abstract_parts.append(f"{label}: {text}")
                else:
                    abstract_parts.append(text)

        authors = []
        author_list = article.find("AuthorList")
        if author_list is not None:
            for author in author_list.findall("Author"):
                last = author.findtext("LastName", "")
                initials = author.findtext("Initials", "")
                if last:
                    authors.append(f"{last} {initials}".strip())

        journal_el = article.find("Journal")
        journal = ""
        year = ""
        if journal_el is not None:
            journal = journal_el.findtext("Title", "") or journal_el.findtext("ISOAbbreviation", "")
            pub_date = journal_el.find("JournalIssue/PubDate")
            if pub_date is not None:
                year = pub_date.findtext("Year", "")

        SESSION.references[pmid] = {"title": title, "authors": authors, "journal": journal, "year": year,
                                    "abstract": "\n".join(abstract_parts)}
        return {
            "pmid": pmid,
            "title": title,
            "authors": authors[:10],
            "journal": journal,
            "year": year,
            "abstract": "\n".join(abstract_parts) if abstract_parts else "",
            "cite_as": f"{{{{cite:{pmid}}}}}",
        }
    except (urllib.error.URLError, OSError, ET.ParseError) as exc:
        return {"error": "fetch_failed", "message": f"Could not fetch abstract: {exc}"}


_COUNT_KEYWORDS = {"count", "counts", "gene_counts", "raw_counts", "featurecounts"}
_DESIGN_KEYWORDS = {"design", "metadata", "coldata", "sample_info", "conditions", "phenotype"}
_TABULAR_EXTENSIONS = {".csv", ".tsv", ".txt"}


def scan_results(results_dir: str) -> Summary:
    """Scan a directory for analysis inputs: count matrices, design files, MultiQC.

    Detects both nf-core/rnaseq output structure and loose user-provided files.
    Non-mutating; the agent calls this first to orient itself.
    """
    rdir = Path(results_dir)
    if not rdir.is_dir():
        return {"error": "not_a_directory", "message": f"'{results_dir}' is not a directory."}

    counts_path = None
    tpm_path = None
    design_path = None
    source = "unknown"

    # --- Try nf-core structure first ---
    star_salmon = rdir / "star_salmon"
    salmon_dir = star_salmon if star_salmon.is_dir() else rdir / "salmon"

    if salmon_dir.is_dir():
        nf_counts = salmon_dir / "salmon.merged.gene_counts.tsv"
        nf_tpm = salmon_dir / "salmon.merged.gene_tpm.tsv"
        if nf_counts.is_file():
            counts_path = nf_counts
            source = "nf-core"
        if nf_tpm.is_file():
            tpm_path = nf_tpm


    # --- Fall back to scanning for loose tabular files ---
    if counts_path is None:
        candidates = []
        for p in sorted(rdir.iterdir()):
            if not p.is_file():
                continue
            stem_lower = p.stem.lower()
            ext = p.suffix.lower()
            if ext not in _TABULAR_EXTENSIONS:
                continue
            if any(kw in stem_lower for kw in _COUNT_KEYWORDS):
                candidates.append(p)
        if len(candidates) == 1:
            counts_path = candidates[0]
            source = "user-provided"
        elif len(candidates) > 1:
            counts_path = candidates[0]
            source = "user-provided"

    # --- Look for design files ---
    for search_dir in (rdir, rdir.parent):
        for p in sorted(search_dir.iterdir()):
            if not p.is_file():
                continue
            stem_lower = p.stem.lower()
            if p.suffix.lower() in _TABULAR_EXTENSIONS and any(kw in stem_lower for kw in _DESIGN_KEYWORDS):
                design_path = p
                break
        if design_path:
            break
    if design_path is None:
        candidate = rdir.parent / "design.csv"
        if candidate.is_file():
            design_path = candidate

    # --- List all tabular files for the agent to inspect if needed ---
    tabular_files = []
    for p in sorted(rdir.iterdir()):
        if p.is_file() and p.suffix.lower() in _TABULAR_EXTENSIONS:
            tabular_files.append(p.name)

    SESSION.results_dir = rdir
    mqc_files = _multiqc_stats_files()

    return {
        "results_dir": str(rdir),
        "source": source,
        "counts_found": counts_path is not None,
        "counts_path": str(counts_path) if counts_path else None,
        "tpm_found": tpm_path is not None,
        "tpm_path": str(tpm_path) if tpm_path else None,
        "multiqc_files": [str(p) for p in mqc_files],
        "design_found": design_path is not None,
        "design_path": str(design_path) if design_path else None,
        "tabular_files": tabular_files,
    }


_FEATURECOUNTS_META = {"Geneid", "Chr", "Start", "End", "Strand", "Length"}


def _detect_separator(path: Path) -> str:
    """Guess CSV vs TSV from the first line."""
    with path.open() as f:
        first_line = f.readline()
    if "\t" in first_line:
        return "\t"
    return ","


def _load_design(design_path: str) -> tuple[pd.DataFrame | None, str | None]:
    """Load a design CSV/TSV. Returns (design_df, error_message)."""
    dp = Path(design_path)
    if not dp.is_file():
        return None, f"Design file '{design_path}' does not exist."
    sep = _detect_separator(dp)
    design = pd.read_csv(dp, sep=sep)
    if "sample" not in design.columns or "condition" not in design.columns:
        return None, "Design file must have 'sample' and 'condition' columns."
    design = design.set_index("sample")
    return design, None


def load_counts(counts_path: str, design_path: str | None = None) -> Summary:
    """Load a gene count matrix and optional design file.

    Handles multiple formats:
    - nf-core/rnaseq: TSV with gene_id, gene_name, then sample columns
    - featureCounts: TSV with Geneid, Chr, Start, End, Strand, Length, then samples
    - Generic CSV/TSV: first column as gene index, remaining as samples
    Auto-detects separator (tab vs comma).
    """
    path = Path(counts_path)
    if not path.is_file():
        return {"error": "not_a_file", "message": f"'{counts_path}' does not exist."}

    sep = _detect_separator(path)
    df = pd.read_csv(path, sep=sep)

    if df.shape[1] < 2:
        return {"error": "bad_format", "message": "Count matrix has fewer than 2 columns."}

    # Identify format and extract the count matrix
    gene_names = None
    gene_id_col = None

    if "gene_id" in df.columns:
        # nf-core format
        gene_id_col = "gene_id"
        if "gene_name" in df.columns:
            gene_names = df["gene_name"]
        meta_cols = {"gene_id", "gene_name"} & set(df.columns)

    elif "Geneid" in df.columns:
        # featureCounts format
        gene_id_col = "Geneid"
        meta_cols = _FEATURECOUNTS_META & set(df.columns)

    else:
        # Generic: first column is gene IDs, rest are samples
        gene_id_col = df.columns[0]
        meta_cols = {gene_id_col}

    gene_ids = df[gene_id_col]
    counts = df.drop(columns=list(meta_cols))

    # Check that remaining columns look numeric
    non_numeric = [c for c in counts.columns if not pd.api.types.is_numeric_dtype(counts[c])]
    if non_numeric:
        return {
            "error": "non_numeric_columns",
            "message": f"Columns {non_numeric[:5]} are not numeric. Check the file format — "
            "the count matrix should have gene IDs in one column and numeric values in the rest.",
        }

    counts = counts.set_index(gene_ids)

    SESSION.counts_df = counts
    SESSION.counts_path = path.resolve()
    SESSION.gene_names = gene_names.set_axis(gene_ids) if gene_names is not None else None

    # Load TPM if available (nf-core convention: same directory)
    nf_tpm = path.parent / "salmon.merged.gene_tpm.tsv"
    if nf_tpm.is_file():
        tpm_df = pd.read_csv(nf_tpm, sep="\t")
        tpm_meta = {"gene_id", "gene_name"} & set(tpm_df.columns)
        tpm_id_col = "gene_id" if "gene_id" in tpm_df.columns else tpm_df.columns[0]
        tpm = tpm_df.drop(columns=list(tpm_meta)).set_index(tpm_df[tpm_id_col])
        SESSION.tpm_df = tpm

    samples = list(counts.columns)
    lib_sizes = {s: round(float(counts[s].sum()), 2) for s in samples}

    # Load design if provided
    conditions = None
    if design_path:
        design, err = _load_design(design_path)
        if err:
            return {"error": "bad_design", "message": err}
        SESSION.design_df = design
        conditions = sorted(design["condition"].unique().tolist())

    detected_format = "nf-core" if "gene_id" in df.columns else (
        "featureCounts" if "Geneid" in df.columns else "generic"
    )

    return {
        "n_genes": int(counts.shape[0]),
        "n_samples": len(samples),
        "samples": samples,
        "library_sizes": lib_sizes,
        "has_gene_names": gene_names is not None,
        "tpm_loaded": SESSION.tpm_df is not None,
        "design_loaded": SESSION.design_df is not None,
        "conditions": conditions,
        "design_columns": list(SESSION.design_df.columns) if SESSION.design_df is not None else None,
        "detected_format": detected_format,
    }


def inspect_counts() -> Summary:
    """Summarise the loaded count matrix so the agent can assess normalisation status.

    Returns per-sample and global statistics: value range, fraction of non-integer
    values, fraction of zeros, and distribution quantiles. These let the agent
    determine whether the data is raw counts, TPM/FPKM, or log-transformed.
    """
    if SESSION.counts_df is None:
        return {"error": "counts_not_loaded", "message": "Call load_counts before inspect_counts."}

    counts = SESSION.counts_df
    values = counts.values.ravel()
    n_values = values.size
    n_zero = int((values == 0).sum())
    n_non_integer = int(np.sum(values != np.floor(values)))

    global_stats = {
        "min": round(float(np.nanmin(values)), 4),
        "max": round(float(np.nanmax(values)), 4),
        "mean": round(float(np.nanmean(values)), 4),
        "median": round(float(np.nanmedian(values)), 4),
        "fraction_zero": round(n_zero / n_values, 4) if n_values else 0,
        "fraction_non_integer": round(n_non_integer / n_values, 4) if n_values else 0,
        "n_genes": int(counts.shape[0]),
        "n_samples": int(counts.shape[1]),
    }

    per_sample = {}
    for col in counts.columns:
        col_vals = counts[col].values
        per_sample[col] = {
            "mean": round(float(np.nanmean(col_vals)), 2),
            "median": round(float(np.nanmedian(col_vals)), 2),
            "max": round(float(np.nanmax(col_vals)), 2),
            "fraction_zero": round(float((col_vals == 0).sum() / len(col_vals)), 4),
        }

    hints = []
    salmon_fractional = False
    if global_stats["fraction_non_integer"] < 0.01 and global_stats["max"] > 100:
        hints.append("Values are almost entirely integers with a wide range — likely raw counts.")
    elif 0.01 <= global_stats["fraction_non_integer"] <= 0.5 and global_stats["max"] > 100:
        salmon_fractional = True
        hints.append(
            "Mostly integers with some fractional values (typical of Salmon/kallisto "
            "probabilistic quantification). These are raw counts suitable for DESeq2 — "
            "fractional values will be rounded to integers before DE analysis."
        )
    elif global_stats["fraction_non_integer"] > 0.5 and global_stats["max"] < 25:
        hints.append("Mostly non-integer values in a narrow range — likely log-transformed. "
                      "DESeq2 requires raw counts; do not run it on this data.")
    elif global_stats["fraction_non_integer"] > 0.5 and global_stats["max"] > 100:
        hints.append("Non-integer values with a wide range — likely TPM/FPKM normalised. "
                      "DESeq2 requires raw counts; do not run it on this data.")

    return {
        "global": global_stats,
        "per_sample": per_sample,
        "hints": hints,
        "salmon_fractional": salmon_fractional,
    }


def set_design(rows: list[dict]) -> Summary:
    """Define the sample-to-condition mapping when no design.csv exists.

    The agent infers conditions from sample names and the user's prompt, then calls
    this to set the mapping. Validates that sample names match the loaded count matrix.
    """
    if SESSION.counts_df is None:
        return {"error": "counts_not_loaded", "message": "Call load_counts before set_design."}
    if not rows:
        return {"error": "empty_design", "message": "Design must have at least one row."}

    for i, row in enumerate(rows):
        if "sample" not in row:
            return {"error": "missing_sample", "message": f"Row {i}: missing 'sample' key."}
        if "condition" not in row:
            return {"error": "missing_condition", "message": f"Row {i}: missing 'condition' key."}

    design = pd.DataFrame(rows).set_index("sample")
    count_samples = set(SESSION.counts_df.columns)
    design_samples = set(design.index)

    missing = count_samples - design_samples
    if missing:
        return {
            "error": "samples_missing_from_design",
            "message": f"Samples in counts but not in design: {sorted(missing)}",
        }

    extra = design_samples - count_samples
    if extra:
        design = design.loc[design.index.isin(count_samples)]

    SESSION.design_df = design

    # Write to disk for traceability
    paths = SESSION.require_paths()
    design.reset_index().to_csv(paths.design, index=False)

    conditions = sorted(design["condition"].unique().tolist())
    return {
        "n_samples": len(design),
        "conditions": conditions,
        "design_path": str(paths.design),
    }


def compute_qc() -> Summary:
    """Library sizes, gene detection rates, and PCA for outlier detection.

    Performs PCA on log2(counts + 1) so the agent can see whether samples cluster
    by condition or whether there are outliers. Parses MultiQC general stats if the
    results directory was scanned.
    """
    if SESSION.counts_df is None:
        return {"error": "counts_not_loaded", "message": "Call load_counts before compute_qc."}

    counts = SESSION.counts_df
    samples = list(counts.columns)

    lib_sizes = {s: int(counts[s].sum()) for s in samples}
    genes_detected = {s: int((counts[s] > 0).sum()) for s in samples}
    SESSION.qc_snapshot = {"library_sizes": pd.Series(lib_sizes), "genes_detected": pd.Series(genes_detected),
                           "filtered": SESSION.filter_settings is not None}

    # PCA on log2(counts + 1)
    log_counts = np.log2(counts.values.astype(float).T + 1)  # samples x genes
    centered = log_counts - log_counts.mean(axis=0)
    U, S, _Vt = np.linalg.svd(centered, full_matrices=False)
    n_comps = min(2, len(S))
    pca_coords = (U[:, :n_comps] * S[:n_comps]).tolist()
    var_explained = ((S ** 2) / (S ** 2).sum())[:n_comps].tolist()
    pca_summary = {
        samples[i]: {"PC1": round(pca_coords[i][0], 4), "PC2": round(pca_coords[i][1], 4) if n_comps > 1 else 0.0}
        for i in range(len(samples))
    }

    return {
        "library_sizes": lib_sizes,
        "library_size_distribution": _quantile_summary(list(lib_sizes.values())),
        "genes_detected": genes_detected,
        "pca": pca_summary,
        "variance_explained": [round(v, 4) for v in var_explained],
        "pca_note": "This PCA is on unfiltered counts. The report's PCA figure and "
                    "{{qc.pca_pc1_pct}}/{{qc.pca_pc2_pct}} use the filtered matrix — cite those.",
    }


def _multiqc_stats_files() -> list[Path]:
    """Every MultiQC general-stats table under results/ — the folder layout differs
    between nf-core/rnaseq versions and aligners, so search rather than assume."""
    if not SESSION.results_dir:
        return []
    return sorted(p.resolve() for p in SESSION.results_dir.rglob("multiqc_general_stats.txt"))


def read_multiqc(path: str | None = None) -> Summary:
    """MultiQC general statistics (mapping rate, duplication, GC, trimming, ...) with
    whatever columns this nf-core version reports. Requires scan_results first."""
    if SESSION.results_dir is None:
        return {"error": "no_results_dir", "message": "Call scan_results first."}
    files = _multiqc_stats_files()
    if not files:
        return {"error": "not_found", "message": "No multiqc_general_stats.txt under the results directory."}
    if path is None and len(files) > 1:
        return {"error": "choose_file", "message": "Several MultiQC tables found; pass one as path.",
                "candidates": [str(p) for p in files]}
    chosen = files[0] if path is None else Path(path).resolve()
    if chosen not in files:
        return {"error": "unknown_file", "message": f"'{path}' is not one of the MultiQC tables found.",
                "candidates": [str(p) for p in files]}

    table = pd.read_csv(chosen, sep="\t").set_index("Sample")
    # MultiQC adds extra rows (e.g. per-read "S1 Read 1"); keep only this analysis's samples.
    samples = list(SESSION.counts_df.columns) if SESSION.counts_df is not None else list(table.index)
    matched = table.loc[[s for s in samples if s in table.index]].dropna(axis=1, how="all")
    SESSION.multiqc_stats = matched

    columns: dict[str, Any] = {}
    for col in matched.columns:
        values = matched[col].dropna()
        if pd.api.types.is_numeric_dtype(values):
            columns[col] = ({s: round(float(v), 2) for s, v in values.items()} if len(values) <= 24
                            else _quantile_summary(values.tolist()))
        else:
            columns[col] = sorted(set(map(str, values)))
    return {
        "path": str(chosen),
        "n_samples_matched": len(matched),
        "samples_missing": [s for s in samples if s not in table.index],
        "columns": columns,
        "facts": [f"multiqc.{_key(c)}.min|median|max" for c in matched.select_dtypes("number").columns],
    }


def filter_low_counts(min_count: int = 10, min_samples: int = 2) -> Summary:
    """Drop genes where fewer than min_samples have at least min_count reads.

    Must run before DE so the multiple-testing correction isn't diluted by
    uninformative genes.
    """
    if SESSION.counts_df is None:
        return {"error": "counts_not_loaded", "message": "Call load_counts before filter_low_counts."}

    counts = SESSION.counts_df
    before = int(counts.shape[0])
    keep = (counts >= min_count).sum(axis=1) >= min_samples
    SESSION.counts_df = counts.loc[keep]
    after = int(SESSION.counts_df.shape[0])

    SESSION.filter_settings = {
        "min_count": min_count,
        "min_samples": min_samples,
        "genes_before": before,
        "genes_after": after,
    }
    return {**SESSION.filter_settings, "genes_removed": before - after}


def run_deseq2(contrast: list[str], covariates: list[str] | None = None) -> Summary:
    """Run PyDESeq2 for a given contrast, optionally adjusting for covariates.

    contrast is [factor, test, reference], e.g. ["condition", "treated", "control"].
    covariates are other design columns (batch, donor, tissue) added before the factor:
    design = ~ cov1 + cov2 + factor. Covariates are treated as categorical.
    PyDESeq2 expects raw integer counts as samples (rows) x genes (cols).
    """
    if SESSION.counts_df is None:
        return {"error": "counts_not_loaded", "message": "Call load_counts before run_deseq2."}
    if SESSION.design_df is None:
        return {"error": "design_not_set", "message": "Call set_design or load_counts with a design before run_deseq2."}
    if len(contrast) != 3:
        return {"error": "bad_contrast", "message": "Contrast must be [factor, test, reference]."}

    factor, test, ref = contrast
    covariates = list(covariates or [])
    for col in [factor, *covariates]:
        if col not in SESSION.design_df.columns:
            return {"error": "bad_factor", "message": f"'{col}' is not a column in the design. "
                    f"Columns: {list(SESSION.design_df.columns)}"}
    if factor in covariates:
        return {"error": "bad_covariates", "message": f"'{factor}' is the contrast factor; don't list it as a covariate."}
    bad_names = [c for c in [factor, *covariates] if not str(c).isidentifier()]
    if bad_names:
        return {"error": "bad_column_name", "message": f"Design column names must be identifiers for the formula: {bad_names}"}

    missing_samples = [s for s in SESSION.counts_df.columns if s not in SESSION.design_df.index]
    if missing_samples:
        return {"error": "samples_missing_from_design", "message": f"Samples in counts but not design: {missing_samples[:10]}"}
    design = SESSION.design_df.loc[SESSION.counts_df.columns, [factor, *covariates]].astype(str)

    levels = sorted(design[factor].unique())
    for level in (test, ref):
        if level not in levels:
            return {"error": "bad_level", "message": f"'{level}' is not a level of '{factor}'. Levels: {levels}"}
    single = [c for c in covariates if design[c].nunique() < 2]
    if single:
        return {"error": "constant_covariate", "message": f"Covariates with a single value add nothing: {single}"}

    # Full-rank check: confounded covariates (e.g. batch == condition) make the model unfittable
    X = pd.get_dummies(design[[*covariates, factor]], drop_first=True).astype(float)
    X.insert(0, "intercept", 1.0)
    rank = int(np.linalg.matrix_rank(X.to_numpy()))
    if rank < X.shape[1]:
        return {"error": "confounded_design", "message": (
            f"Covariates {covariates} are confounded with '{factor}' (design matrix rank {rank} < "
            f"{X.shape[1]} columns). Drop the confounded covariate — its effect can't be separated.")}
    if len(design) - rank < 1:
        return {"error": "no_residual_df", "message": "Too many parameters for the number of samples."}
    formula = "~" + " + ".join([*covariates, factor])

    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    # PyDESeq2 expects samples x genes with integer counts
    counts_t = SESSION.counts_df.T

    # Salmon/kallisto produce fractional counts from probabilistic assignment.
    # Round to integers only when the data looks like raw counts with minor
    # fractional noise (not normalised data). Detection mirrors inspect_counts.
    counts_rounded = False
    frac_non_int = float(np.mean(counts_t.values != np.floor(counts_t.values)))
    if frac_non_int > 0.005 and float(counts_t.values.max()) > 100:
        counts_t = counts_t.round().astype(int)
        counts_rounded = True

    dds = DeseqDataSet(counts=counts_t, metadata=design, design=formula)
    dds.deseq2()

    stat_res = DeseqStats(dds, contrast=[factor, test, ref])
    stat_res.summary()
    results = stat_res.results_df.copy()

    SESSION.deseq_results = results
    SESSION.deseq_design = formula
    SESSION.deseq_contrast = [factor, test, ref]
    SESSION.deseq_rounded_pct = round(frac_non_int * 100, 1) if counts_rounded else None

    # Write to disk
    paths = SESSION.require_paths()
    paths.analysis_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(paths.de_results)

    sig = results[results["padj"] < 0.05].dropna(subset=["padj"])
    n_up = int((sig["log2FoldChange"] > 0).sum())
    n_down = int((sig["log2FoldChange"] < 0).sum())

    # Top genes by padj
    top_up = _ranked("up").head(5)
    top_down = _ranked("down").head(5)

    def _gene_row(gene_id, row):
        out = {"gene_id": str(gene_id), "log2FC": round(float(row["log2FoldChange"]), 4), "padj": float(row["padj"])}
        if SESSION.gene_names is not None and gene_id in SESSION.gene_names.index:
            out["gene_name"] = str(SESSION.gene_names.loc[gene_id])
        return out

    # Check replication per group
    group_sizes = design[factor].value_counts()
    min_reps = int(group_sizes.min())

    result: Summary = {
        "contrast": contrast,
        "design": formula,
        "n_tested": int(results["padj"].notna().sum()),
        "n_significant": int(len(sig)),
        "n_up": n_up,
        "n_down": n_down,
        "top_up": [_gene_row(gid, row) for gid, row in top_up.iterrows()],
        "top_down": [_gene_row(gid, row) for gid, row in top_down.iterrows()],
        "de_results_path": str(paths.de_results),
        "replicates_per_group": group_sizes.to_dict(),
    }
    if min_reps < 3:
        result["low_replication_warning"] = (
            f"Only {min_reps} replicates in the smallest group. "
            "Dispersion estimates are unreliable with < 3 replicates — "
            "treat DE results as exploratory, not confirmatory."
        )
    if counts_rounded:
        result["counts_rounded"] = True
        result["pct_non_integer_before_rounding"] = round(frac_non_int * 100, 1)
    return result


def _ranked(direction: str, padj_max: float = _DE_PADJ, lfc_min: float = 0.0) -> pd.DataFrame:
    """DE genes passing padj < padj_max and |log2FC| >= lfc_min in one direction, ranked
    by padj, then |log2FC| (largest first), then gene ID — fully deterministic, so ties
    (e.g. many genes with padj = 0) always resolve the same way. Used by every tool that
    picks "top" genes, so tables, heatmap and enrichment agree."""
    res = SESSION.deseq_results.dropna(subset=["padj"])
    lfc = res["log2FoldChange"]
    keep = (res["padj"] < padj_max) & (lfc.abs() >= lfc_min) & ((lfc > 0) if direction == "up" else (lfc < 0))
    ranked = res[keep].assign(_abs_lfc=lfc[keep].abs(), _gene_id=res.index[keep].astype(str))
    return ranked.sort_values(["padj", "_abs_lfc", "_gene_id"], ascending=[True, False, True]).drop(
        columns=["_abs_lfc", "_gene_id"])


def _gene_label(gene_id) -> str:
    if SESSION.gene_names is not None and gene_id in SESSION.gene_names.index:
        name = SESSION.gene_names.loc[gene_id]
        if isinstance(name, str) and name:
            return name
    return str(gene_id)


def get_top_genes(n: int = 20, direction: str = "both") -> Summary:
    """Return top DE genes (padj < 0.05) for the agent to inspect, ranked by padj then |log2FC|.

    direction: "up" (log2FC > 0), "down" (log2FC < 0), or "both". Informational only —
    enrichment and report tables select genes themselves.
    """
    if SESSION.deseq_results is None:
        return {"error": "no_deseq_results", "message": "Call run_deseq2 before get_top_genes."}
    if direction not in ("up", "down", "both"):
        return {"error": "bad_direction", "message": "direction must be 'up', 'down' or 'both'."}

    directions = ["up", "down"] if direction == "both" else [direction]
    top = pd.concat([_ranked(d) for d in directions])
    if direction == "both":
        top = top.sort_values("padj", kind="stable")
    genes = []
    for gene_id, row in top.head(n).iterrows():
        entry = {
            "gene_id": str(gene_id),
            "log2FC": round(float(row["log2FoldChange"]), 4),
            "padj": float(row["padj"]),
            "baseMean": round(float(row["baseMean"]), 2),
        }
        if SESSION.gene_names is not None and gene_id in SESSION.gene_names.index:
            entry["gene_name"] = str(SESSION.gene_names.loc[gene_id])
        genes.append(entry)

    return {
        "direction": direction,
        "n_returned": len(genes),
        "n_significant_total": int(len(top)),
        "genes": genes,
    }


def _gene_labels() -> pd.Series:
    """Display symbol for every gene in the DESeq2 results, indexed by gene ID."""
    res = SESSION.deseq_results
    return pd.Series([_gene_label(g) for g in res.index], index=res.index)


def _gene_status(row) -> str:
    if pd.isna(row["padj"]):
        return "not tested (removed by independent filtering)"
    if row["padj"] >= _DE_PADJ:
        return "not significant"
    return "up" if row["log2FoldChange"] > 0 else "down"


def query_genes(symbols: list[str] | None = None, prefix: str | None = None) -> Summary:
    """DE results for any genes: a list of symbols, or every gene whose symbol starts with
    a prefix (a gene family, e.g. "ITG" for integrins). Check a gene or family here before
    making a claim about it."""
    if SESSION.deseq_results is None:
        return {"error": "no_deseq_results", "message": "Call run_deseq2 before query_genes."}
    if not symbols and not prefix:
        return {"error": "nothing_to_query", "message": "Pass symbols, a prefix, or both."}
    res, labels = SESSION.deseq_results, _gene_labels()
    loaded = set(SESSION.gene_names.dropna().astype(str)) if SESSION.gene_names is not None else set()

    def row(gene_id) -> dict[str, Any]:
        r = res.loc[gene_id]
        return {"gene": labels[gene_id], "log2FC": round(float(r["log2FoldChange"]), 3),
                "padj": None if pd.isna(r["padj"]) else float(r["padj"]),
                "baseMean": round(float(r["baseMean"]), 1), "status": _gene_status(r)}

    out: Summary = {}
    if symbols:
        found, missing = [], {}
        for sym in symbols:
            hits = labels.index[labels.str.upper() == str(sym).upper()]
            if len(hits):
                found += [row(g) for g in hits]
            else:
                missing[sym] = ("filtered out before DE (low counts)" if sym in loaded or sym in SESSION.counts_df.index
                                else "not in the data (check the symbol)")
        out["genes"] = found
        out["cite_as"] = "{{gene:SYMBOL}} for each gene"
        if missing:
            out["not_found"] = missing
    if prefix:
        family = labels.index[labels.str.upper().str.startswith(prefix.upper())]
        rows = sorted((row(g) for g in family), key=lambda r: (r["padj"] is None, r["padj"] or 1.0))
        statuses = pd.Series([r["status"] for r in rows], dtype=str)
        out["family"] = {
            "prefix": prefix,
            "n_genes": len(rows),
            "n_up": int((statuses == "up").sum()),
            "n_down": int((statuses == "down").sum()),
            "n_not_significant": int((statuses == "not significant").sum()),
            "genes": rows[:100],
            "truncated": len(rows) > 100,
            "cite_as": f"{{{{genes:{prefix}*}}}}",
        }
    return out


def _enrichment_table(label: str) -> pd.DataFrame | None:
    """Enrichr's full results for one direction, as saved by run_enrichment."""
    path = SESSION.require_paths().analysis_dir / f"enrichment_{label}.csv"
    return pd.read_csv(path) if path.is_file() else None


def search_enrichment(query: str, direction: str | None = None) -> Summary:
    """Search every enriched term (not only the top 10) for a word, e.g. "platelet" or
    "hematopoietic", with its rank, statistics and the genes behind it."""
    labels = {"up": ["upregulated"], "down": ["downregulated"], None: ["upregulated", "downregulated"]}
    if direction not in labels:
        return {"error": "bad_direction", "message": "direction must be 'up', 'down' or omitted."}
    matches = []
    for label in labels[direction]:
        table = _enrichment_table(label)
        if table is None:
            continue
        for gs, sub in table.groupby("Gene_set"):
            sub = sub.sort_values("Adjusted P-value").reset_index(drop=True)
            for rank, r in sub.iterrows():
                if query.lower() in str(r["Term"]).lower():
                    genes = str(r["Genes"]).split(";")
                    matches.append({
                        "direction": "up" if label == "upregulated" else "down",
                        "library": gs, "term": r["Term"], "rank_in_library": rank + 1,
                        "overlap": r["Overlap"], "padj": float(r["Adjusted P-value"]),
                        "genes": genes[:50], "n_genes": len(genes),
                        "cite_as": f"{{{{term:{'up' if label == 'upregulated' else 'down'}:{r['Term']}}}}}",
                    })
    if not matches and not any(_enrichment_table(l) is not None for l in labels[direction]):
        return {"error": "no_enrichment", "message": "Call run_enrichment first."}
    matches.sort(key=lambda m: m["padj"])
    return {"query": query, "n_matches": len(matches), "matches": matches[:20]}


def _abstract_genes() -> list[tuple[str, list[str]]]:
    """Gene symbols named in the fetched abstracts, with the PMIDs naming them, in the
    order they first appear. Only exact matches to genes in the loaded data count, and
    only all-caps or digit-containing words (so ordinary words aren't taken for genes)."""
    if not SESSION.references or SESSION.counts_df is None:
        return []
    universe = set(map(str, SESSION.counts_df.index))
    if SESSION.gene_names is not None:
        universe |= set(SESSION.gene_names.dropna().astype(str))
    found: dict[str, list[str]] = {}
    for pmid, ref in SESSION.references.items():
        for token in re.findall(r"[A-Za-z0-9][A-Za-z0-9-]*[A-Za-z0-9]", ref.get("abstract") or ""):
            if token in universe and (token.isupper() or any(c.isdigit() for c in token)):
                found.setdefault(token, [])
                if pmid not in found[token]:
                    found[token].append(pmid)
    return list(found.items())


def _paper_genes() -> list[dict[str, str]]:
    """This analysis's result for every gene the fetched abstracts name."""
    if SESSION.deseq_results is None:
        return []
    res, labels = SESSION.deseq_results, _gene_labels()
    rows = []
    for symbol, pmids in _abstract_genes():
        hits = labels.index[labels == symbol]
        if len(hits) == 0:
            hits = [g for g in res.index if str(g) == symbol]
        base = {"gene": symbol, "pmids": ", ".join(pmids)}
        if len(hits) == 0:
            rows.append({**base, "log2FC": "", "padj": "", "result": "filtered out before DE (low counts)"})
            continue
        r = res.loc[hits[0]]
        rows.append({**base, "log2FC": f"{r['log2FoldChange']:.2f}",
                     "padj": "" if pd.isna(r["padj"]) else _fmt_p(r["padj"]), "result": _gene_status(r)})
    return rows


def run_enrichment(
    direction: str,
    padj_max: float = _DE_PADJ,
    lfc_min: float = 1.0,
    max_genes: int = 500,
    organism: str = "human",
) -> Summary:
    """Over-representation analysis via Enrichr (GO and KEGG) on DE genes selected here.

    Genes are chosen by the tool from the DESeq2 results — never passed in — using
    padj < padj_max and |log2FC| >= lfc_min in the given direction, ranked (padj, then
    |log2FC|, then ID) and capped at max_genes. IDs are converted to gene symbols. The
    exact input list and Enrichr's raw results are saved to analysis/.
    """
    if SESSION.deseq_results is None:
        return {"error": "no_deseq_results", "message": "Call run_deseq2 before run_enrichment."}
    if direction not in ("up", "down"):
        return {"error": "bad_direction", "message": "direction must be 'up' or 'down'."}
    organism = organism.lower()
    if organism not in _ENRICHR_LIBRARIES:
        return {"error": "bad_organism", "message": f"organism must be one of {sorted(_ENRICHR_LIBRARIES)}."}
    if max_genes < 1:
        return {"error": "bad_max_genes", "message": "max_genes must be at least 1."}

    passing = _ranked(direction, padj_max, lfc_min)
    selected = passing.head(max_genes)
    symbols, n_without_symbol = [], 0
    for gene_id in selected.index:
        name = _gene_label(gene_id)
        if re.match(r"^(ENS[A-Z]*G\d+(\.\d+)?|\d+)$", name):
            n_without_symbol += 1  # Ensembl/Entrez ID with no symbol: Enrichr can't use it
            continue
        symbols.append(name)
    symbols = list(dict.fromkeys(symbols))
    if not symbols:
        return {"error": "no_genes", "message": (
            f"No {direction}-regulated genes pass padj < {padj_max} and |log2FC| >= {lfc_min} "
            f"({len(passing)} before symbol mapping). Consider a lower lfc_min (e.g. 0.585 = 1.5-fold).")}

    import gseapy

    gene_sets = _ENRICHR_LIBRARIES[organism]
    try:
        enr = gseapy.enrichr(gene_list=symbols, gene_sets=gene_sets, organism=organism,
                             outdir=None, no_plot=True)
    except Exception as exc:
        return {"error": "enrichr_failed", "message": str(exc)}
    results_df = enr.results

    label = "upregulated" if direction == "up" else "downregulated"
    paths = SESSION.require_paths()
    paths.analysis_dir.mkdir(parents=True, exist_ok=True)
    (paths.analysis_dir / f"enrichment_input_{label}.txt").write_text("\n".join(symbols) + "\n")
    results_df.to_csv(paths.analysis_dir / f"enrichment_{label}.csv", index=False)

    enrichment = {}
    for gs in gene_sets:
        subset = results_df[results_df["Gene_set"] == gs].nsmallest(10, "Adjusted P-value")
        enrichment[gs] = [
            {
                "term": str(row["Term"]),
                "pval": float(row["P-value"]),
                "padj": float(row["Adjusted P-value"]),
                "overlap": str(row["Overlap"]),
                "genes": str(row["Genes"]),
            }
            for _, row in subset.iterrows()
        ]

    selection = {
        "direction": direction,
        "padj_max": padj_max,
        "lfc_min": lfc_min,
        "max_genes": max_genes,
        "n_passing": int(len(passing)),
        "capped": bool(len(passing) > max_genes),
        "n_without_symbol": n_without_symbol,
    }
    if SESSION.enrichment_results is None:
        SESSION.enrichment_results = {}
    SESSION.enrichment_results[label] = {
        "organism": organism,
        "gene_sets": gene_sets,
        "n_input_genes": len(symbols),
        "selection": selection,
        "results": enrichment,
    }

    return {
        "label": label,
        "organism": organism,
        "gene_sets_queried": gene_sets,
        "n_input_genes": len(symbols),
        "selection": selection,
        "input_genes_file": str(paths.analysis_dir / f"enrichment_input_{label}.txt"),
        "enrichment": enrichment,
    }


def _data_provenance() -> tuple[str, dict[str, Any]]:
    """Where the loaded counts came from, from files on disk — never from the LLM."""
    counts_path = SESSION.counts_path
    paths = SESSION.paths

    for subdir, method in (("star_salmon", "STAR alignment + Salmon quantification"),
                           ("salmon", "Salmon pseudo-alignment (alignment skipped)")):
        quant = (SESSION.results_dir / subdir) if SESSION.results_dir else None
        if quant and counts_path and quant.resolve() in counts_path.parents:
            return "nf-core/rnaseq", {"quantification": method, "counts_file": str(counts_path)}

    if counts_path and paths and paths.counts_metadata.is_file() and counts_path == paths.counts_matrix.resolve():
        meta = _json.loads(paths.counts_metadata.read_text())
        return "GEO count matrix", {
            "accession": meta["accession"],
            "source": "authors' supplementary file" if meta["source"] == "author"
                      else "NCBI-generated counts",
            "file": meta["filename"],
            "md5": meta["md5"],
            "value_type": meta["value_type"],
            "gene_id_type": meta["gene_id_type"],
            "n_duplicate_gene_ids_summed": meta["n_duplicates_summed"],
        }

    return "user-provided", {"counts_file": str(counts_path) if counts_path else None}


def summarize_findings() -> Summary:
    """Consolidate the final analysis state into one factual summary for the report.

    Returns source type, methods metadata, and all accumulated results so the
    agent has the facts it needs to write an accurate report.
    """
    out: Summary = {}
    out["data_source"], out["data_provenance"] = _data_provenance()
    if SESSION.filter_settings:
        out["filtering"] = SESSION.filter_settings

    if SESSION.counts_df is not None:
        out["n_genes"] = int(SESSION.counts_df.shape[0])
        out["n_samples"] = int(SESSION.counts_df.shape[1])

    if SESSION.design_df is not None:
        out["conditions"] = sorted(SESSION.design_df["condition"].unique().tolist())
        out["samples_per_condition"] = SESSION.design_df["condition"].value_counts().to_dict()

    if SESSION.deseq_results is not None:
        sig = SESSION.deseq_results.dropna(subset=["padj"])
        sig = sig[sig["padj"] < 0.05]
        out["de_summary"] = {
            "n_tested": int(SESSION.deseq_results["padj"].notna().sum()),
            "n_significant": int(len(sig)),
            "n_up": int((sig["log2FoldChange"] > 0).sum()),
            "n_down": int((sig["log2FoldChange"] < 0).sum()),
            "design_formula": SESSION.deseq_design,
        }

    if SESSION.enrichment_results:
        out["enrichment_inputs"] = {
            label: {"n_input_genes": e.get("n_input_genes"), "gene_sets": e.get("gene_sets"),
                    "organism": e.get("organism"), "selection": e.get("selection")}
            for label, e in SESSION.enrichment_results.items()
        }
        top_terms = {}
        for label, entry in SESSION.enrichment_results.items():
            results = entry.get("results", {})
            for gs, terms in results.items():
                key = f"{label}:{gs}" if label != "default" else gs
                top_terms[key] = [t["term"] for t in terms[:5]]
        out["top_enrichment_terms"] = top_terms

    if (paper := _paper_genes()):
        out["paper_genes"] = paper

    versions = _software_versions()
    pipeline = _pipeline_versions()
    if pipeline:
        out["pipeline_versions"] = pipeline
    out["software_versions"] = versions

    # What the report may cite: placeholder names with their current rendered values.
    out["facts"] = _facts()
    out["tables_available"] = sorted(_tables())
    out["tables_required"] = _required_tables()
    out["references"] = {pmid: f"{{{{cite:{pmid}}}}}" for pmid in SESSION.references}
    return out


def _software_versions() -> dict[str, str]:
    versions = {"python": __import__("sys").version.split()[0]}
    for pkg in ("numpy", "pandas", "pydeseq2", "matplotlib", "scipy", "gseapy"):
        try:
            versions[pkg] = str(__import__(pkg).__version__)
        except (ImportError, AttributeError):
            pass
    return versions


_DISPLAY_NAMES = {"python": "Python", "numpy": "NumPy", "pandas": "pandas", "pydeseq2": "PyDESeq2",
                  "matplotlib": "Matplotlib", "scipy": "SciPy", "gseapy": "GSEApy", "salmon": "Salmon"}


def _pipeline_versions() -> dict[str, str]:
    if not SESSION.results_dir:
        return {}
    versions_yml = SESSION.results_dir / "pipeline_info" / "nf_core_rnaseq_software_mqc_versions.yml"
    if not versions_yml.is_file():
        return {}
    try:
        import yaml
        data = yaml.safe_load(versions_yml.read_text()) or {}
    except Exception:
        return {}
    out = {k: str(v) for k, v in (data.get("Workflow") or {}).items()}
    salmon = (data.get("SALMON_QUANT") or {}).get("salmon")
    if salmon:
        out["salmon"] = str(salmon)
    tximeta = (data.get("TXIMETA_TXIMPORT") or {}).get("bioconductor-tximeta")
    if tximeta:
        out["tximeta"] = str(tximeta)
    return out


def _fmt_int(n) -> str:
    return f"{int(n):,}"


def _fmt_pct(part, whole) -> str:
    return f"{100 * part / whole:.1f}%" if whole else "n/a"


def _fmt_p(p) -> str:
    # padj of exactly 0 is floating-point underflow, not a true zero
    return "< 1e-300" if p == 0 else (f"{p:.2e}" if p < 0.001 else f"{p:.3f}")


def _key(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(text).lower()).strip("_")


def _library_stats() -> tuple[pd.Series, pd.Series]:
    """Per-sample total counts and genes detected, from the compute_qc snapshot (before
    low-count filtering) when available, else from the current matrix. Used by the facts,
    the QC table and the library-size figure so they always agree."""
    if SESSION.qc_snapshot is not None:
        return SESSION.qc_snapshot["library_sizes"], SESSION.qc_snapshot["genes_detected"]
    counts = SESSION.counts_df
    return counts.sum(axis=0), (counts > 0).sum(axis=0)


def _pca() -> dict[str, Any]:
    """PCA of log2(counts + 1) on the current (filtered) matrix — the single computation
    behind the PCA figures, their captions and the qc.pca_* facts."""
    counts = SESSION.counts_df
    log_counts = np.log2(counts.values.astype(float).T + 1)
    centered = log_counts - log_counts.mean(axis=0)
    U, S, _Vt = np.linalg.svd(centered, full_matrices=False)
    n_pcs = min(10, len(S))
    return {"coords": U[:, :n_pcs] * S[:n_pcs], "var_exp": (S ** 2) / (S ** 2).sum(),
            "samples": list(counts.columns), "n_pcs": n_pcs}


def _sample_correlation() -> pd.DataFrame:
    return np.log2(SESSION.counts_df.astype(float) + 1).corr(method="pearson")


def _qc_before_filtering() -> bool:
    """Library sizes and genes detected come from a compute_qc run before filtering."""
    return (SESSION.qc_snapshot is not None and not SESSION.qc_snapshot["filtered"]
            and SESSION.filter_settings is not None)


def _pca_groups() -> dict[str, Any] | None:
    """Replicate spread and group separation on PC1-PC2 of the PCA figure, so the report
    describes clustering from numbers rather than by eye."""
    counts, design = SESSION.counts_df, SESSION.design_df
    if counts is None or design is None or "condition" not in design.columns or counts.shape[1] < 3:
        return None
    pca = _pca()
    if pca["n_pcs"] < 2:
        return None
    xy = pd.DataFrame(pca["coords"][:, :2], index=pca["samples"], columns=["PC1", "PC2"])
    cond = design.loc[xy.index, "condition"].astype(str)
    if cond.nunique() < 2:
        return None
    centroids = xy.groupby(cond).mean()

    def dist(sample, group) -> float:
        return float(np.linalg.norm(xy.loc[sample] - centroids.loc[group]))

    spread = {g: float(np.mean([dist(s, g) for s in xy.index[cond == g]])) for g in centroids.index}
    pairs = [(a, b, float(np.linalg.norm(centroids.loc[a] - centroids.loc[b])))
             for i, a in enumerate(centroids.index) for b in centroids.index[i + 1:]]
    closest = min(pairs, key=lambda p: p[2])
    misclustered = {s: min(centroids.index, key=lambda g: dist(s, g)) for s in xy.index}
    misclustered = {s: g for s, g in misclustered.items() if g != cond[s]}
    r2 = {}
    for pc in ("PC1", "PC2"):
        v = xy[pc]
        ss_total = float(((v - v.mean()) ** 2).sum())
        ss_between = sum(len(v[cond == g]) * (v[cond == g].mean() - v.mean()) ** 2 for g in centroids.index)
        r2[pc] = ss_between / ss_total if ss_total else 0.0
    return {"spread": spread, "closest_pair": closest, "misclustered": misclustered, "condition_r2": r2}


def _pca_summary(groups: dict[str, Any]) -> str:
    spread = groups["spread"]
    a, b, d = groups["closest_pair"]
    what = "group centroids are" if len(spread) == 2 else f"the closest group centroids ({a}, {b}) are"
    mis = ", ".join(f"{s} (closer to {g})" for s, g in groups["misclustered"].items()) or "none"
    text = ("On PC1–PC2, mean distance of replicates to their group centroid: "
            + ", ".join(f"{g} {v:.1f}" for g, v in spread.items()) + f"; {what} {d:.1f} apart.")
    if min(spread.values()) > 0:
        text += f" The most spread group is {max(spread.values()) / min(spread.values()):.1f}× the least spread."
    return (text + f" Condition explains {groups['condition_r2']['PC1']:.0%} of PC1 and "
            f"{groups['condition_r2']['PC2']:.0%} of PC2 variance. Samples closer to another group's "
            f"centroid: {mis}.")


def _reference_provenance() -> dict[str, str]:
    """The genome and annotation nf-core actually used, from its own params file."""
    if not SESSION.results_dir:
        return {}
    params_files = sorted((SESSION.results_dir / "pipeline_info").glob("params_*.json"))
    if not params_files:
        return {}
    params = _json.loads(params_files[-1].read_text())
    out: dict[str, str] = {}
    if params.get("genome"):
        out["genome"] = str(params["genome"])
    gtf = params.get("gtf") or params.get("gff")
    if gtf:
        source = str(gtf)
        cached = Path(source).parent / "reference.json"
        if "://" not in source and cached.is_file():  # shared reference cache (agents/submission/reference.py)
            out["cached_reference"] = str(cached.parent)
            source = _json.loads(cached.read_text()).get("gtf_source") or source
        out["annotation_file"] = source
        m = re.search(r"/igenomes/+([^/]+)/([^/]+)/([^/]+)/", source)
        if m:
            out["annotation"] = f"iGenomes {m.group(1).replace('_', ' ')} {m.group(2)} {m.group(3)}"
            out["annotation_provider"] = m.group(2)
        elif params.get("gencode"):
            out["annotation"] = out["annotation_provider"] = "GENCODE"
    return out


def _methods() -> str:
    """The Methods section, written from what the tools actually did. The LLM places it
    with {{table:methods}} and never describes the methods itself."""
    f = _facts()
    source, prov = _data_provenance()
    ref = _reference_provenance()
    sections: list[tuple[str, str]] = []

    if source == "nf-core/rnaseq":
        text = (f"Reads were processed with {f.get('pipeline.nf_core_rnaseq', 'nf-core/rnaseq')}"
                + (f" ({f['pipeline.nextflow']})" if "pipeline.nextflow" in f else "") + ". "
                f"Quantification: {prov['quantification']}"
                + (f" ({f['pipeline.salmon']})" if "pipeline.salmon" in f else "") + "; the pipeline "
                "summarised transcript estimates to gene-level counts"
                + (f" with {f['pipeline.tximeta']}" if "pipeline.tximeta" in f else "")
                + f", in `{Path(prov['counts_file']).name}`.")
        if ref.get("genome"):
            text += f" Reference genome: `{ref['genome']}`."
        if ref.get("annotation_file"):
            text += (f" Gene annotation: {ref.get('annotation', 'GTF')} (`{ref['annotation_file']}`)"
                     + (", reused from the shared reference cache" if "cached_reference" in ref else "") + ".")
        sections.append(("Data and quantification", text))
    elif source == "GEO count matrix":
        sections.append(("Data", (
            f"Counts: {prov['source']} `{prov['file']}` from GEO {prov['accession']} (MD5 `{prov['md5']}`); "
            f"values: {str(prov['value_type']).replace('_', ' ')}; gene IDs: {prov['gene_id_type']}; "
            f"{prov['n_duplicate_gene_ids_summed']} duplicate gene IDs summed.")))
    else:
        sections.append(("Data", f"Counts were provided by the user (`{prov.get('counts_file')}`)."))

    if SESSION.counts_df is not None:
        when = " before low-count filtering" if _qc_before_filtering() else ""
        matrix = (f"the filtered matrix ({f['filter.genes_after']} genes)" if SESSION.filter_settings
                  else "all genes")
        text = (f"Library size is the total count per sample and genes detected the number of genes with "
                f"at least one count, both computed{when or ' on the loaded matrix'}. PCA: singular value "
                f"decomposition of gene-centred log2(counts + 1) of {matrix}, without scaling. Sample "
                f"correlation: Pearson correlation of log2(counts + 1) of {matrix}.")
        if SESSION.multiqc_stats is not None:
            text += " Read-level metrics come from the pipeline's MultiQC general statistics."
        sections.append(("Quality control", text))

    if SESSION.filter_settings:
        sections.append(("Low-count filtering", (
            f"Genes were kept when at least {f['filter.min_samples']} samples had ≥ {f['filter.min_count']} "
            f"counts: {f['filter.genes_after']} of {f['filter.genes_before']} genes kept "
            f"({f['filter.pct_removed']} removed).")))

    if SESSION.deseq_results is not None:
        factor, test, ref_level = SESSION.deseq_contrast
        covariates = [c.strip() for c in str(SESSION.deseq_design).lstrip("~").split("+")][:-1]
        text = ""
        if SESSION.deseq_rounded_pct is not None:
            text += (f"Non-integer estimated counts ({SESSION.deseq_rounded_pct}% of values) were rounded "
                     "to integers. ")
        text += (
            f"Differential expression was tested with {f.get('versions.pydeseq2', 'PyDESeq2')} using the "
            f"design `{SESSION.deseq_design}`"
            + (f" (covariates {', '.join(covariates)} treated as categorical)" if covariates else "")
            + f", comparing {test} with {ref_level} (reference level). PyDESeq2 fits a negative binomial "
            "generalised linear model per gene, with median-of-ratios size factors and dispersions "
            f"shrunk towards a fitted trend, and tests the {factor} coefficient with a Wald test. "
            "P-values were adjusted with the Benjamini–Hochberg method after PyDESeq2's default "
            f"independent filtering and Cook's-distance outlier handling; {f['de.n_tested']} genes "
            f"received an adjusted p-value. Genes with padj < {_DE_PADJ} are called significant. "
            "Log2 fold changes are unshrunken maximum-likelihood estimates.")
        sections.append(("Differential expression", text))

    for label, e in (SESSION.enrichment_results or {}).items():
        sel = e.get("selection") or {}
        direction = {"upregulated": "up-regulated", "downregulated": "down-regulated"}.get(label, label)
        text = (f"Over-representation analysis of {direction} genes used Enrichr through "
                f"{f.get('versions.gseapy', 'GSEApy')} against {', '.join(e['gene_sets'])} ({e['organism']}). ")
        if sel:
            text += (f"Genes with padj < {sel['padj_max']} and |log2FC| ≥ {sel['lfc_min']} "
                     f"({_fmt_int(sel['n_passing'])} genes) were ranked by padj, then |log2FC|, and the top "
                     f"{_fmt_int(min(sel['n_passing'], sel['max_genes']))} submitted as "
                     f"{_fmt_int(e['n_input_genes'])} gene symbols")
            text += (f" ({sel['n_without_symbol']} without a symbol excluded). " if sel.get("n_without_symbol")
                     else ". ")
        text += ("Enrichr's default background (all genes in each library) was used, not the genes tested "
                 "here. Terms are reported by Enrichr's adjusted p-value (Benjamini–Hochberg). Exact input: "
                 f"`analysis/enrichment_input_{label}.txt`.")
        sections.append((f"Enrichment ({direction})", text))

    versions = {**{_DISPLAY_NAMES.get(k, k): v for k, v in _software_versions().items()},
                **{f"nf-core: {_DISPLAY_NAMES.get(k, k)}": v for k, v in _pipeline_versions().items()}}
    sections.append(("Software", _md_table(["Software", "Version"], [[k, v] for k, v in versions.items()])))
    sections.append(("Reproducibility", (
        "Every analysis step is logged in `analysis/tool_calls.jsonl`, and `analysis/replay.py` re-runs "
        "them without the LLM. This Methods section is generated from those steps, not written by the LLM.")))
    return "\n\n".join(f"### {title}\n\n{body}" for title, body in sections)


def _library_key(gene_set: str) -> str:
    if gene_set.startswith("GO_Biological_Process"):
        return "go_bp"
    if gene_set.startswith("KEGG"):
        return "kegg"
    return _key(gene_set)


def _facts() -> dict[str, str]:
    """Every value the report may cite, rendered as text, keyed by placeholder name.

    The single source for both summarize_findings (what the LLM reads) and write_report
    (what gets printed), so the two can never disagree.
    """
    f: dict[str, str] = {}
    source, provenance = _data_provenance()
    f["provenance.data_source"] = source
    for k, v in provenance.items():
        if v is not None:
            f[f"provenance.{k}"] = str(v).replace("_", " ") if k == "value_type" else str(v)

    if SESSION.counts_df is not None:
        counts = SESSION.counts_df
        lib, detected = _library_stats()
        f["samples.n"] = _fmt_int(counts.shape[1])
        f["qc.library_size_min_millions"] = f"{lib.min() / 1e6:.1f}"
        f["qc.library_size_max_millions"] = f"{lib.max() / 1e6:.1f}"
        f["qc.library_size_median_millions"] = f"{lib.median() / 1e6:.1f}"
        f["qc.genes_detected_min"] = _fmt_int(detected.min())
        f["qc.genes_detected_max"] = _fmt_int(detected.max())
        if counts.shape[1] >= 2:
            pca = _pca()
            f["qc.pca_pc1_pct"] = f"{pca['var_exp'][0]:.1%}"
            if pca["n_pcs"] > 1:
                f["qc.pca_pc2_pct"] = f"{pca['var_exp'][1]:.1%}"
            f["qc.sample_correlation_min"] = f"{_sample_correlation().values.min():.2f}"
        groups = _pca_groups()
        if groups:
            for g, v in groups["spread"].items():
                f[f"qc.pca_spread.{_key(g)}"] = f"{v:.1f}"
            f["qc.pca_centroid_distance"] = f"{groups['closest_pair'][2]:.1f}"
            f["qc.pca_condition_r2_pc1"] = f"{groups['condition_r2']['PC1']:.0%}"
            f["qc.pca_condition_r2_pc2"] = f"{groups['condition_r2']['PC2']:.0%}"
            f["qc.pca_misclustered"] = ", ".join(groups["misclustered"]) or "none"
            f["qc.pca_summary"] = _pca_summary(groups)
    for k, v in _reference_provenance().items():
        if k in ("genome", "annotation", "annotation_provider"):
            f[f"reference.{k}"] = v
    if SESSION.multiqc_stats is not None:
        for col in SESSION.multiqc_stats.select_dtypes("number").columns:
            values = SESSION.multiqc_stats[col].dropna()
            for stat in ("min", "median", "max"):
                f[f"multiqc.{_key(col)}.{stat}"] = f"{getattr(values, stat)():,.2f}"
    if SESSION.design_df is not None and SESSION.counts_df is not None:
        design = SESSION.design_df.loc[[c for c in SESSION.counts_df.columns if c in SESSION.design_df.index]]
        per = design["condition"].value_counts()
        f["samples.conditions"] = ", ".join(map(str, per.index))
        for cond, n in per.items():
            f[f"samples.n_{_key(cond)}"] = _fmt_int(n)
        # Other design columns: constant ones are facts about every sample (e.g. time point,
        # cell type); varying ones list their levels. Identifiers are skipped.
        for col in design.columns:
            if col in ("condition", "gsm", "title") or design[col].nunique() == len(design) > 1:
                continue
            values = list(dict.fromkeys(design[col].astype(str)))
            if len(values) == 1:
                f[f"design.{_key(col)}"] = values[0]
            else:
                f[f"design.{_key(col)}_values"] = ", ".join(values)
    if SESSION.filter_settings:
        fs = SESSION.filter_settings
        f["filter.min_count"] = _fmt_int(fs["min_count"])
        f["filter.min_samples"] = _fmt_int(fs["min_samples"])
        f["filter.genes_before"] = _fmt_int(fs["genes_before"])
        f["filter.genes_after"] = _fmt_int(fs["genes_after"])
        f["filter.genes_removed"] = _fmt_int(fs["genes_before"] - fs["genes_after"])
        f["filter.pct_removed"] = _fmt_pct(fs["genes_before"] - fs["genes_after"], fs["genes_before"])

    if SESSION.deseq_results is not None:
        res = SESSION.deseq_results.dropna(subset=["padj"])
        sig = res[res["padj"] < _DE_PADJ]
        n_tested, n_sig = len(res), len(sig)
        n_up = int((sig["log2FoldChange"] > 0).sum())
        factor, test, ref = SESSION.deseq_contrast or ("condition", "?", "?")
        f.update({
            "de.contrast": f"{test} vs {ref}",
            "de.test_level": str(test),
            "de.reference_level": str(ref),
            "de.factor": str(factor),
            "de.design": str(SESSION.deseq_design),
            "de.padj_threshold": str(_DE_PADJ),
            "de.n_tested": _fmt_int(n_tested),
            "de.n_significant": _fmt_int(n_sig),
            "de.pct_significant": _fmt_pct(n_sig, n_tested),
            "de.n_up": _fmt_int(n_up),
            "de.pct_up": _fmt_pct(n_up, n_tested),
            "de.n_down": _fmt_int(n_sig - n_up),
            "de.pct_down": _fmt_pct(n_sig - n_up, n_tested),
        })

    for label, e in (SESSION.enrichment_results or {}).items():
        sel = e.get("selection") or {}
        pre = f"enrichment.{'up' if label == 'upregulated' else 'down' if label == 'downregulated' else _key(label)}"
        f[f"{pre}.n_input_genes"] = _fmt_int(e.get("n_input_genes", 0))
        f[f"{pre}.gene_sets"] = ", ".join(e.get("gene_sets", []))
        f[f"{pre}.organism"] = str(e.get("organism"))
        for gs, terms in e.get("results", {}).items():
            for rank, term in enumerate(terms, start=1):
                f[f"{pre}.{_library_key(gs)}.{rank}"] = (
                    f"{term['term']} ({term['overlap']} genes, padj {_fmt_p(term['padj'])})")
        if sel:
            f[f"{pre}.padj_max"] = str(sel["padj_max"])
            f[f"{pre}.lfc_min"] = str(sel["lfc_min"])
            f[f"{pre}.max_genes"] = _fmt_int(sel["max_genes"])
            f[f"{pre}.n_passing"] = _fmt_int(sel["n_passing"])

    # Versions render with the tool name ("PyDESeq2 0.5.4"); write_report drops the name
    # again when the prose already says it.
    for k, v in _software_versions().items():
        f[f"versions.{k}"] = f"{_DISPLAY_NAMES.get(k, k)} {v}"
    for k, v in _pipeline_versions().items():
        f[f"pipeline.{_key(k)}"] = f"{_DISPLAY_NAMES.get(k, k)} {v}"
    return f


def _md_table(header: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(str(c).replace("|", "\\|") for c in row) + " |" for row in rows]
    return "\n".join(lines)


def _tables() -> dict[str, str]:
    """Report tables generated from results (placed with {{table:<name>}})."""
    t: dict[str, str] = {}
    facts = _facts()
    if SESSION.counts_df is not None:
        counts = SESSION.counts_df
        lib, detected = _library_stats()
        rows = []
        for s_ in counts.columns:
            cond = (str(SESSION.design_df.loc[s_, "condition"])
                    if SESSION.design_df is not None and s_ in SESSION.design_df.index else "")
            rows.append([s_, cond, _fmt_int(lib[s_]), _fmt_int(detected[s_])])
        t["qc"] = _md_table(["Sample", "Condition", "Total counts", "Genes detected"], rows) + (
            "\n\n*Genes detected: genes with at least one count"
            + (", before low-count filtering" if _qc_before_filtering() else "") + ".*")
    if SESSION.deseq_results is not None:
        t["de_summary"] = _md_table(["Measure", "Value"], [
            ["Contrast", facts["de.contrast"]],
            ["Design", f"`{facts['de.design']}`"],
            ["Genes tested", facts["de.n_tested"]],
            [f"Significant (padj < {_DE_PADJ})", f"{facts['de.n_significant']} ({facts['de.pct_significant']})"],
            ["Up-regulated", f"{facts['de.n_up']} ({facts['de.pct_up']})"],
            ["Down-regulated", f"{facts['de.n_down']} ({facts['de.pct_down']})"],
        ])
        for direction in ("up", "down"):
            top = _ranked(direction).head(_TOP_TABLE_ROWS)
            show_id = any(_gene_label(g) != str(g) for g in top.index)  # only when IDs aren't symbols
            t[f"top_{direction}"] = _md_table(
                ["Gene", *(["Gene ID"] if show_id else []), "log2FC", "padj", "baseMean"],
                [[_gene_label(g), *([str(g)] if show_id else []), f"{r['log2FoldChange']:.2f}",
                  _fmt_p(r["padj"]), f"{r['baseMean']:.0f}"]
                 for g, r in top.iterrows()],
            )
    for label, e in (SESSION.enrichment_results or {}).items():
        name = {"upregulated": "enrichment_up", "downregulated": "enrichment_down"}.get(label, f"enrichment_{_key(label)}")
        rows = [[gs.replace("_", " "), term["term"], term["overlap"], _fmt_p(term["padj"])]
                for gs, terms in e.get("results", {}).items() for term in terms]
        t[name] = _md_table(["Library", "Term", "Overlap", "padj"], rows) if rows else "_No enriched terms returned._"
    versions = {**_software_versions(), **{f"nf-core: {k}": v for k, v in _pipeline_versions().items()}}
    t["versions"] = _md_table(["Software", "Version"], [[k, v] for k, v in versions.items()])
    if (paper := _paper_genes()):
        t["paper_genes"] = _md_table(
            ["Gene", "Named in (PMID)", "log2FC", "padj", "This analysis"],
            [[r["gene"], r["pmids"], r["log2FC"], r["padj"], r["result"]] for r in paper],
        ) + ("\n\n*Genes named in the fetched abstract(s), with this analysis's result "
             f"({SESSION.deseq_contrast[1]} vs {SESSION.deseq_contrast[2]}).*")
    if SESSION.counts_df is not None:
        t["methods"] = _methods()
    return t


def _required_tables() -> list[str]:
    """Tables the report must place, given what the analysis produced."""
    req = ["methods"] if SESSION.counts_df is not None else []
    if SESSION.deseq_results is not None:
        req += ["de_summary", "top_up", "top_down"]
    if _paper_genes():
        req.append("paper_genes")
    for label in SESSION.enrichment_results or {}:
        req.append({"upregulated": "enrichment_up", "downregulated": "enrichment_down"}.get(label, f"enrichment_{_key(label)}"))
    return req


def _informative_columns() -> list[str]:
    """Design columns other than condition that can carry information in a plot: more
    than one value, and not unique per sample (identifiers such as gsm or title)."""
    if SESSION.design_df is None or SESSION.counts_df is None:
        return []
    design = SESSION.design_df.loc[[s for s in SESSION.counts_df.columns if s in SESSION.design_df.index]]
    return [c for c in design.columns if c != "condition" and 1 < design[c].nunique() < len(design)]


def _draw_figures(pca_color_by: list[str]) -> list[dict]:
    """Draw every figure the current analysis state supports.

    Returns [{"path": "figures/<name>.png", "caption": <facts about the figure>}].
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable

    paths = SESSION.require_paths()
    paths.analysis_figures.mkdir(parents=True, exist_ok=True)
    figures: list[dict] = []

    # --- Consistent style + palette ---
    _PALETTE = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
                "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD"]
    _SIG_COLOR = "#C44E52"
    _NS_COLOR = "#B0B0B0"

    def _style_ax(ax, title=""):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=10)
        if title:
            ax.set_title(title, fontsize=13, fontweight="bold", pad=10)

    def _save(fig, name, caption):
        fig.savefig(paths.analysis_figures / name, dpi=150, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        figures.append({"path": f"figures/{name}", "caption": caption})

    def _condition_colors(conditions_series):
        unique = list(dict.fromkeys(conditions_series))
        return {c: _PALETTE[i % len(_PALETTE)] for i, c in enumerate(unique)}

    # --- PCA (shared with the qc.pca_* facts) ---
    pca_data = _pca() if SESSION.counts_df is not None and SESSION.counts_df.shape[1] >= 2 else None

    # --- Library size bar plot ---
    if SESSION.counts_df is not None:
        counts = SESSION.counts_df
        samples = list(counts.columns)
        lib_series, _detected = _library_stats()
        lib_sizes = [float(lib_series[s]) for s in samples]

        fig, ax = plt.subplots(figsize=(max(6, len(samples) * 0.8), 5))
        colors = [_NS_COLOR] * len(samples)
        if SESSION.design_df is not None:
            design_samples = [s for s in samples if s in SESSION.design_df.index]
            cmap = _condition_colors(SESSION.design_df.loc[design_samples, "condition"])
            colors = [cmap.get(SESSION.design_df.loc[s, "condition"], _NS_COLOR) if s in SESSION.design_df.index else _NS_COLOR for s in samples]
        bars = ax.bar(range(len(samples)), [s / 1e6 for s in lib_sizes], color=colors, edgecolor="white", linewidth=0.5)
        ax.set_xticks(range(len(samples)))
        ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=9)
        ax.set_ylabel("Library size (millions)", fontsize=11)
        _style_ax(ax, "Library sizes")
        if SESSION.design_df is not None:
            for cond, color in cmap.items():
                ax.bar([], [], color=color, label=cond)
            ax.legend(fontsize=9, frameon=False)
        _save(fig, "library_sizes.png",
              f"Total counts per sample{' before low-count filtering' if _qc_before_filtering() else ''} "
              f"({len(samples)} samples, {min(lib_sizes) / 1e6:.1f}–{max(lib_sizes) / 1e6:.1f} million), "
              f"coloured by condition.")

    # --- PCA by condition ---
    if pca_data is not None:
        coords = pca_data["coords"]
        var_exp = pca_data["var_exp"]
        samples = pca_data["samples"]

        fig, ax = plt.subplots(figsize=(7, 6))
        if SESSION.design_df is not None:
            conditions = SESSION.design_df.loc[samples, "condition"]
            cmap = _condition_colors(conditions)
            for cond in dict.fromkeys(conditions):
                mask = [i for i, c in enumerate(conditions) if c == cond]
                ax.scatter(coords[mask, 0], coords[mask, 1] if pca_data["n_pcs"] > 1 else np.zeros(len(mask)),
                           c=cmap[cond], label=cond, s=70, edgecolors="white", linewidths=0.5, zorder=3)
            ax.legend(fontsize=9, frameon=False)
        else:
            ax.scatter(coords[:, 0], coords[:, 1] if pca_data["n_pcs"] > 1 else np.zeros(len(samples)),
                       c=_PALETTE[0], s=70, edgecolors="white", linewidths=0.5, zorder=3)
        for i, s in enumerate(samples):
            ax.annotate(s, (coords[i, 0], coords[i, 1] if pca_data["n_pcs"] > 1 else 0),
                        fontsize=8, alpha=0.7, textcoords="offset points", xytext=(5, 5))
        ax.set_xlabel(f"PC1 ({var_exp[0]:.1%} variance)", fontsize=11)
        ax.set_ylabel(f"PC2 ({var_exp[1]:.1%} variance)" if pca_data["n_pcs"] > 1 else "PC2", fontsize=11)
        _style_ax(ax, "PCA — log2(counts + 1)")
        _save(fig, "pca.png",
              f"PCA of log2(counts + 1): PC1 {var_exp[0]:.1%}"
              + (f", PC2 {var_exp[1]:.1%}" if pca_data["n_pcs"] > 1 else "")
              + " of variance, coloured by condition."
              + (f" {_pca_summary(groups)}" if (groups := _pca_groups()) else ""))

    # --- PCA coloured by other metadata variables ---
    if pca_data is not None and SESSION.design_df is not None:
        for col in pca_color_by:
            vals = SESSION.design_df.loc[samples, col]
            fig, ax = plt.subplots(figsize=(7, 6))

            if vals.dtype == object or vals.nunique() <= 8:
                unique_vals = list(dict.fromkeys(vals))
                col_cmap = {v: _PALETTE[i % len(_PALETTE)] for i, v in enumerate(unique_vals)}
                for val in unique_vals:
                    mask = [i for i, v in enumerate(vals) if v == val]
                    ax.scatter(coords[mask, 0], coords[mask, 1], c=col_cmap[val],
                               label=str(val), s=70, edgecolors="white", linewidths=0.5, zorder=3)
                ax.legend(fontsize=9, frameon=False, title=col)
            else:
                numeric_vals = pd.to_numeric(vals, errors="coerce")
                sc = ax.scatter(coords[:, 0], coords[:, 1], c=numeric_vals, cmap="viridis",
                                s=70, edgecolors="white", linewidths=0.5, zorder=3)
                plt.colorbar(sc, ax=ax, label=col)

            for i, s in enumerate(samples):
                ax.annotate(s, (coords[i, 0], coords[i, 1]),
                            fontsize=8, alpha=0.7, textcoords="offset points", xytext=(5, 5))
            ax.set_xlabel(f"PC1 ({var_exp[0]:.1%} variance)", fontsize=11)
            ax.set_ylabel(f"PC2 ({var_exp[1]:.1%} variance)", fontsize=11)
            _style_ax(ax, f"PCA coloured by {col}")
            _save(fig, f"pca_{col.lower().replace(' ', '_')}.png",
                  f"PCA of log2(counts + 1) coloured by {col} ({vals.nunique()} values).")

    # --- PC–metadata association heatmap ---
    if pca_data is not None and SESSION.design_df is not None:
        design_cols = ["condition", *_informative_columns()]
        n_pcs_show = min(pca_data["n_pcs"], 10)
        if design_cols and n_pcs_show >= 2:
            from scipy import stats as sp_stats
            assoc_matrix = np.full((len(design_cols), n_pcs_show), np.nan)
            for ci, col in enumerate(design_cols):
                vals = SESSION.design_df.loc[samples, col]
                for pc_i in range(n_pcs_show):
                    pc_vals = coords[:, pc_i]
                    numeric_vals = pd.to_numeric(vals, errors="coerce")
                    if numeric_vals.notna().all():
                        r, _ = sp_stats.pearsonr(numeric_vals.values, pc_vals)
                        assoc_matrix[ci, pc_i] = r ** 2
                    elif vals.nunique() > 1:
                        groups = [pc_vals[[j for j in range(len(vals)) if vals.iloc[j] == g]]
                                  for g in vals.unique() if sum(vals == g) > 0]
                        if len(groups) >= 2 and all(len(g) > 0 for g in groups):
                            f_stat, p_val = sp_stats.f_oneway(*groups)
                            ss_between = sum(len(g) * (g.mean() - pc_vals.mean()) ** 2 for g in groups)
                            ss_total = np.sum((pc_vals - pc_vals.mean()) ** 2)
                            assoc_matrix[ci, pc_i] = ss_between / ss_total if ss_total > 0 else 0

            fig, ax = plt.subplots(figsize=(max(6, n_pcs_show * 0.8), max(3, len(design_cols) * 0.6 + 1)))
            im = ax.imshow(assoc_matrix, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)
            ax.set_xticks(range(n_pcs_show))
            ax.set_xticklabels([f"PC{i+1}\n({pca_data['var_exp'][i]:.1%})" for i in range(n_pcs_show)], fontsize=9)
            ax.set_yticks(range(len(design_cols)))
            ax.set_yticklabels(design_cols, fontsize=10)
            for ci in range(len(design_cols)):
                for pi in range(n_pcs_show):
                    val = assoc_matrix[ci, pi]
                    if not np.isnan(val):
                        ax.text(pi, ci, f"{val:.2f}", ha="center", va="center",
                                fontsize=8, color="white" if val > 0.5 else "black")
            plt.colorbar(im, ax=ax, label="R² (association)", shrink=0.8)
            _style_ax(ax, "PC–metadata associations")
            ax.spines["left"].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            _save(fig, "pc_association.png",
                  f"Association (R²) between the first {n_pcs_show} principal components and "
                  f"design variables: {', '.join(design_cols)}.")

    # --- Sample correlation heatmap ---
    if SESSION.counts_df is not None:
        corr = _sample_correlation()
        samples = list(corr.columns)
        n = len(samples)

        fig, ax = plt.subplots(figsize=(max(5, n * 0.6), max(4, n * 0.5)))
        im = ax.imshow(corr.values, cmap="RdYlBu_r", vmin=corr.values.min(), vmax=1)
        ax.set_xticks(range(n))
        ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=9)
        ax.set_yticks(range(n))
        ax.set_yticklabels(samples, fontsize=9)
        for i in range(n):
            for j in range(n):
                val = corr.values[i, j]
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7,
                        color="white" if val < 0.95 else "black")
        plt.colorbar(im, ax=ax, label="Pearson r", shrink=0.8)
        _style_ax(ax, "Sample correlation — log2(counts + 1)")
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        _save(fig, "sample_correlation.png",
              f"Pearson correlation of log2(counts + 1) between all {n} samples "
              f"(range {corr.values.min():.2f}–1.00).")

    # --- Volcano plot ---
    if SESSION.deseq_results is not None:
        results = SESSION.deseq_results.dropna(subset=["padj"])
        sig_mask = results["padj"] < 0.05

        fig, ax = plt.subplots(figsize=(7, 6))
        ax.scatter(results.loc[~sig_mask, "log2FoldChange"],
                   -np.log10(results.loc[~sig_mask, "padj"]),
                   c=_NS_COLOR, alpha=0.4, s=12, label="NS", zorder=2)
        ax.scatter(results.loc[sig_mask, "log2FoldChange"],
                   -np.log10(results.loc[sig_mask, "padj"]),
                   c=_SIG_COLOR, alpha=0.6, s=12, label="padj < 0.05", zorder=3)
        ax.axhline(-np.log10(0.05), color="#888888", linewidth=0.8, linestyle="--", alpha=0.5)
        ax.set_xlabel("log2 Fold Change", fontsize=11)
        ax.set_ylabel("-log10(padj)", fontsize=11)
        ax.legend(fontsize=9, frameon=False)
        _style_ax(ax, "Volcano plot")
        n_sig = int(sig_mask.sum())
        n_up = int((results.loc[sig_mask, "log2FoldChange"] > 0).sum())
        de_facts = (f"{len(results)} genes tested; {n_sig} with padj < 0.05 "
                    f"({n_up} up, {n_sig - n_up} down; design {SESSION.deseq_design}).")
        _save(fig, "volcano.png", f"Volcano plot: {de_facts}")

        # --- MA plot ---
        fig, ax = plt.subplots(figsize=(7, 6))
        ax.scatter(np.log10(results.loc[~sig_mask, "baseMean"] + 1),
                   results.loc[~sig_mask, "log2FoldChange"],
                   c=_NS_COLOR, alpha=0.4, s=12, label="NS", zorder=2)
        ax.scatter(np.log10(results.loc[sig_mask, "baseMean"] + 1),
                   results.loc[sig_mask, "log2FoldChange"],
                   c=_SIG_COLOR, alpha=0.6, s=12, label="padj < 0.05", zorder=3)
        ax.axhline(0, color="#333333", linewidth=0.8)
        ax.set_xlabel("log10(baseMean + 1)", fontsize=11)
        ax.set_ylabel("log2 Fold Change", fontsize=11)
        ax.legend(fontsize=9, frameon=False)
        _style_ax(ax, "MA plot")
        _save(fig, "ma_plot.png", f"MA plot: {de_facts}")

    # --- Top DE genes heatmap ---
    if SESSION.deseq_results is not None and SESSION.counts_df is not None:
        sig = SESSION.deseq_results.dropna(subset=["padj"])
        sig = sig[sig["padj"] < 0.05]
        if len(sig) > 0:
            up = _ranked("up").head(_HEATMAP_PER_DIRECTION).index
            down = _ranked("down").head(_HEATMAP_PER_DIRECTION).index
            top_genes = list(up) + list(down)
            n_top = len(top_genes)
            counts_top = SESSION.counts_df.loc[top_genes]
            log_vals = np.log2(counts_top.values.astype(float) + 1)
            row_means = log_vals.mean(axis=1, keepdims=True)
            row_stds = log_vals.std(axis=1, keepdims=True)
            row_stds[row_stds == 0] = 1
            z_scores = (log_vals - row_means) / row_stds

            gene_labels = []
            for gid in top_genes:
                if SESSION.gene_names is not None and gid in SESSION.gene_names.index:
                    gene_labels.append(str(SESSION.gene_names.loc[gid]))
                else:
                    gene_labels.append(str(gid))

            heatmap_samples = list(counts_top.columns)
            fig, ax = plt.subplots(figsize=(max(5, len(heatmap_samples) * 0.7), max(6, n_top * 0.25)))
            im = ax.imshow(z_scores, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
            ax.set_xticks(range(len(heatmap_samples)))
            ax.set_xticklabels(heatmap_samples, rotation=45, ha="right", fontsize=9)
            ax.set_yticks(range(n_top))
            ax.set_yticklabels(gene_labels, fontsize=7)
            if len(up) and len(down):
                ax.axhline(len(up) - 0.5, color="black", linewidth=1)
            plt.colorbar(im, ax=ax, label="z-score", shrink=0.6)
            _style_ax(ax, f"Top {len(up)} up / {len(down)} down DE genes (z-scored)")
            ax.spines["left"].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            _save(fig, "de_heatmap.png",
                  f"Top {len(up)} up-regulated (above the line) and {len(down)} down-regulated "
                  f"genes with padj < 0.05, ranked by padj then |log2FC|; z-scored log2(counts + 1) "
                  f"across {len(heatmap_samples)} samples.")

    # --- Enrichment bar plots (one per label × gene set) ---
    if SESSION.enrichment_results:
        palette_idx = 0
        for label, entry in SESSION.enrichment_results.items():
            results = entry.get("results", {})
            for gs, terms in results.items():
                if not terms:
                    continue
                top = terms[:10]
                fig, ax = plt.subplots(figsize=(9, max(3, len(top) * 0.4)))
                term_names = [t["term"][:65] for t in reversed(top)]
                pvals = [-np.log10(t["padj"]) if t["padj"] > 0 else 10 for t in reversed(top)]
                color = _PALETTE[palette_idx % len(_PALETTE)]
                ax.barh(term_names, pvals, color=color, edgecolor="white", linewidth=0.5)
                ax.set_xlabel("-log10(padj)", fontsize=11)
                title_label = label if label != "default" else ""
                title_gs = gs.replace("_", " ")
                plot_title = f"{title_gs} — {title_label}" if title_label else title_gs
                _style_ax(ax, plot_title)
                if label != "default":
                    fig_name = f"enrichment_{label}_{gs.lower().replace(' ', '_')}.png"
                else:
                    fig_name = f"enrichment_{gs.lower().replace(' ', '_')}.png"
                n_in = entry.get("n_input_genes")
                _save(fig, fig_name,
                      f"Enrichr {gs}, {title_label or 'input'} genes"
                      + (f" ({n_in} input genes)" if n_in else "")
                      + f": top {len(top)} terms by adjusted p-value.")
            palette_idx += 1

    return figures


_IMAGE_LINK = re.compile(r"!\[([^\]]*)\]\(([^)\s]+)\)")
_BARE_IMAGE = re.compile(r"!\[([^\]]+)\](?!\()")  # ![figures/x.png] — missing the (path) part
_DISCLAIMER = (
    "\n\n---\n\n"
    "**Disclaimer:** This report was generated by an AI system. "
    "Large language models can produce inaccurate statements (hallucinations). "
    "All biological claims, gene annotations, pathway interpretations, and "
    "cited references should be independently verified before use in "
    "publications or clinical decisions. "
    "Biological roles listed alongside gene names are AI-generated summaries "
    "and must be verified against primary databases (UniProt, NCBI Gene)."
)


def generate_figures(pca_color_by: list[str] | None = None) -> Summary:
    """Draw all figures for the current analysis state and return their paths and
    factual captions. Call before writing the report so it describes real figures.

    pca_color_by: design columns to draw extra PCA plots for. Defaults to every
    informative column (identifier and constant columns are never plotted).
    """
    if SESSION.counts_df is None:
        return {"error": "counts_not_loaded", "message": "Call load_counts before generate_figures."}
    informative = _informative_columns()
    if pca_color_by is None:
        pca_color_by = informative
    pca_color_by = [c for c in pca_color_by if c != "condition"]  # always drawn as pca.png
    bad = [c for c in pca_color_by if c not in informative]
    if bad:
        return {"error": "bad_pca_columns", "message": (
            f"Cannot colour PCA by {bad}: not a design column, or constant / unique per sample "
            f"(identifiers). Options: {informative}")}

    SESSION.figures = _draw_figures(pca_color_by)
    SESSION.report_attempts = 0
    return {
        "figures": {f["path"]: f["caption"] for f in SESSION.figures},
        "pca_color_options": informative,
        "note": "Reference figures only by these paths. write_report adds each figure's "
                "path and caption beneath it automatically.",
    }


_PLACEHOLDER = re.compile(r"\{\{\s*([^{}]+?)\s*\}\}")
_HEADING = re.compile(r"^(#{1,6})\s+(.*)$", re.MULTILINE)
_METHODS_ONLY = re.compile(r"\{\{\s*table:(methods|versions)\s*\}\}")


_INTERPRETATION_MAX_WORDS = 250


def _sections(report_markdown: str, title: str) -> list[tuple[str, str]]:
    """(heading, body) for each section whose heading matches title (a regex); the body
    runs to the next heading of the same or a higher level."""
    headings = list(_HEADING.finditer(report_markdown))
    out = []
    for i, h in enumerate(headings):
        if re.match(rf"(\d+\.?\s*)?({title})\b", h.group(2).strip(), re.IGNORECASE):
            level = len(h.group(1))
            end = next((n.start() for n in headings[i + 1:] if len(n.group(1)) <= level), len(report_markdown))
            out.append((h.group(0).strip(), report_markdown[h.end():end]))
    return out


def _interpretation_problems(report_markdown: str) -> list[str]:
    """Biological Interpretation stays short, and every paragraph or list item points at
    a result through a placeholder, so claims are tied to evidence."""
    problems = []
    for heading, body in _sections(report_markdown, "biological interpretation|interpretation|discussion"):
        units, para = [], []
        for line in body.splitlines() + [""]:
            text = line.strip()
            if re.match(r"([-*+]|\d+\.)\s", text) or not text or text.startswith(("#", "|", "![")):
                if para:
                    units.append(" ".join(para))
                para = []
                if re.match(r"([-*+]|\d+\.)\s", text):
                    para = [text]
                continue
            para.append(text)
        words = len(_PLACEHOLDER.sub("X", body).split())
        if words > _INTERPRETATION_MAX_WORDS:
            problems.append(f"'{heading}' has {words} words; keep it to {_INTERPRETATION_MAX_WORDS}.")
        bare = [u for u in units if not _PLACEHOLDER.search(u)]
        problems += [f"'{heading}': this paragraph points at no result — anchor it with a {{{{gene:}}}}, "
                     f"{{{{genes:PREFIX*}}}}, {{{{term:...}}}}, fact or {{{{cite:}}}} placeholder, or drop it: "
                     f"\"{u[:80]}...\"" for u in bare]
    return problems


def _methods_prose(report_markdown: str) -> list[str]:
    """Headings of Methods sections that contain anything besides {{table:methods}}."""
    headings = list(_HEADING.finditer(report_markdown))
    bad = []
    for i, h in enumerate(headings):
        if not re.match(r"(\d+\.?\s*)?(materials and )?methods\b", h.group(2).strip(), re.IGNORECASE):
            continue
        level = len(h.group(1))
        end = next((n.start() for n in headings[i + 1:] if len(n.group(1)) <= level), len(report_markdown))
        if _METHODS_ONLY.sub("", report_markdown[h.end():end]).strip():
            bad.append(h.group(0).strip())
    return bad
_MAX_REPORT_REJECTIONS = 2


def _render_gene(query: str) -> str | None:
    if SESSION.deseq_results is None:
        return None
    res = SESSION.deseq_results
    if query in res.index:
        gene_id = query
    else:
        names = SESSION.gene_names if SESSION.gene_names is not None else pd.Series(dtype=str)
        hits = [g for g, n in names.items() if isinstance(n, str) and n == query and g in res.index]
        if not hits:
            hits = [g for g, n in names.items()
                    if isinstance(n, str) and n.lower() == query.lower() and g in res.index]
        if len(hits) != 1:
            return None
        gene_id = hits[0]
    row = res.loc[gene_id]
    padj = "n/a (not tested)" if pd.isna(row["padj"]) else _fmt_p(row["padj"])
    return f"{_gene_label(gene_id)} (log2FC {row['log2FoldChange']:.2f}, padj {padj})"


def _render_gene_family(pattern: str) -> str | None:
    if SESSION.deseq_results is None or not pattern.endswith("*") or len(pattern) < 2:
        return None
    labels = _gene_labels()
    family = SESSION.deseq_results.loc[labels.index[labels.str.upper().str.startswith(pattern[:-1].upper())]]
    if family.empty:
        return None
    status = family.apply(_gene_status, axis=1)
    n_up, n_down = int((status == "up").sum()), int((status == "down").sum())
    return (f"{pattern} genes: {n_up + n_down} of {len(family)} significant "
            f"({n_up} up, {n_down} down)")


def _render_term(arg: str) -> str | None:
    direction, _, term = arg.partition(":")
    label = {"up": "upregulated", "down": "downregulated"}.get(direction)
    table = _enrichment_table(label) if label else None
    if table is None:
        return None
    hits = table[table["Term"].str.lower() == term.strip().lower()]
    if hits.empty:
        return None
    r = hits.sort_values("Adjusted P-value").iloc[0]
    return f"{r['Term']} ({r['Overlap']} genes, padj {_fmt_p(r['Adjusted P-value'])})"


def _render_cite(pmid: str) -> str | None:
    ref = SESSION.references.get(pmid)
    if ref is None:
        return None
    authors = ref.get("authors") or []
    surname = " ".join(authors[0].split()[:-1]) or authors[0] if authors else "Anonymous"
    who = surname if len(authors) == 1 else f"{surname} et al."
    return f"({who}, {ref.get('year') or 'n.d.'}; PMID: {pmid})"


def write_report(report_markdown: str) -> Summary:
    """Write the report. Every figure link, number, gene statistic, table and citation
    must come from code:

    - every figure from generate_figures must be placed, as ![description](path);
      a missing figures/ prefix, or ![figures/x.png] without the (path) part, is fixed
      losslessly when it names a real figure; a caption with path and facts is added;
    - {{name}} renders a fact from summarize_findings; {{gene:SYMBOL}} renders the
      gene's log2FC and padj; {{cite:PMID}} renders a fetched reference;
      {{table:name}} inserts a code-generated table; required tables must be placed.

    All problems are reported together and nothing is written. After
    _MAX_REPORT_REJECTIONS rejections the report is written anyway with visible ⚠
    markers where it's wrong. LLM text is never deleted. The standard disclaimer is
    appended.
    """
    if SESSION.figures is None:
        return {"error": "no_figures", "message": "Call generate_figures before write_report."}
    by_name = {Path(f["path"]).name: f for f in SESSION.figures}
    facts, tables = _facts(), _tables()
    forced = SESSION.report_attempts >= _MAX_REPORT_REJECTIONS
    problems: list[str] = []
    referenced: list[str] = []

    def _figure(m: re.Match) -> str:
        target = m.group(2)
        fig = by_name.get(Path(target).name)
        if fig is None:
            problems.append(f"figure not generated: {target}")
            return f"{m.group(0)}\n\n*⚠ Figure not generated: `{target}`*"
        referenced.append(fig["path"])
        return f"![{m.group(1)}]({fig['path']})\n\n*File: `{fig['path']}` — {fig['caption']}*"

    placed_tables: set[str] = set()

    def _placeholder(m: re.Match) -> str:
        key = m.group(1)
        kind, _, arg = key.partition(":")
        if kind == "table" and arg in tables:
            placed_tables.add(arg)
            return f"\n\n{tables[arg]}\n\n"
        rendered = (facts.get(key) if not arg else
                    _render_gene(arg) if kind == "gene" else
                    _render_gene_family(arg) if kind == "genes" else
                    _render_term(arg) if kind == "term" else
                    _render_cite(arg) if kind == "cite" else None)
        if rendered is None:
            problems.append(f"unknown placeholder: {{{{{key}}}}}")
            return f"⚠[unknown: {{{{{key}}}}}]"
        if key.startswith(("versions.", "pipeline.")):
            # Lossless: "PyDESeq2 {{versions.pydeseq2}}" must not print the name twice.
            name, _, version = rendered.rpartition(" ")
            if re.search(rf"{re.escape(name)}\s*\(?\s*$", m.string[:m.start()], re.IGNORECASE):
                return version
        return rendered

    def _bare_image(m: re.Match) -> str:
        # Lossless fix when the bracket names a real figure; otherwise it's a problem.
        fig = by_name.get(Path(m.group(1).strip()).name)
        if fig is not None:
            return f"![{Path(fig['path']).stem}]({fig['path']})"
        problems.append(f"malformed image reference {m.group(0)} — use ![description](figures/<name>.png)")
        return f"{m.group(0)} *⚠ malformed image reference*"

    problems += [f"'{h}' must contain only {{{{table:methods}}}} — the Methods section is generated from "
                 "what the tools did. Move other text (e.g. why a covariate was dropped) to Limitations "
                 "or the relevant results section." for h in _methods_prose(report_markdown)]
    problems += _interpretation_problems(report_markdown)
    report = _BARE_IMAGE.sub(_bare_image, report_markdown)
    report = _IMAGE_LINK.sub(_figure, report)
    report = _PLACEHOLDER.sub(_placeholder, report)
    missing = [t for t in _required_tables() if t not in placed_tables]
    problems += [f"required table not placed: {{{{table:{t}}}}}" for t in missing]
    missing_figs = [f for f in SESSION.figures if f["path"] not in referenced]
    problems += [f"figure not placed: ![...]({f['path']})" for f in missing_figs]

    if problems and not forced:
        SESSION.report_attempts += 1
        return {
            "error": "report_problems",
            "message": f"Nothing was written. Fix ALL of these and resubmit: {problems}",
            "available_facts": sorted(facts),
            "available_tables": sorted(tables),
            "required_figures": [f["path"] for f in SESSION.figures],
            "fetched_references": sorted(SESSION.references),
        }

    if missing:
        report += "\n\n## ⚠ Required tables the report did not place\n\n" + "\n\n".join(
            tables[t] for t in missing)
    if missing_figs:
        report += "\n\n## ⚠ Figures the report did not place\n\n" + "\n\n".join(
            f"![{Path(f['path']).stem}]({f['path']})\n\n*File: `{f['path']}` — {f['caption']}*"
            for f in missing_figs)
    report = re.sub(r"\n{3,}", "\n\n", report)  # blank-line runs left by table insertion

    # Strip trailing horizontal rules to avoid doubling before the disclaimer
    report = report.rstrip()
    while report.endswith("---"):
        report = report[:-3].rstrip()

    paths = SESSION.require_paths()
    paths.analysis_report.write_text(report + _DISCLAIMER + "\n")
    result: Summary = {
        "report_path": str(paths.analysis_report),
        "figures_referenced": len(referenced),
        "unreferenced_figures": [f["path"] for f in SESSION.figures if f["path"] not in referenced],
    }
    if problems:
        result["warning"] = (f"Written after {SESSION.report_attempts} rejections with ⚠ markers — "
                             f"the user must review these: {problems}")
    return result


def generate_report(report_markdown: str) -> Summary:
    """Back-compat for replay scripts from older runs: figures + report in one call."""
    result = generate_figures()
    return result if "error" in result else write_report(report_markdown)

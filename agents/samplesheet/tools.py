"""Sample sheet generation tools.

THE CONTRACT:
  Every tool is a deterministic Python function that does real work AND returns a
  structured, JSON-serializable summary dict. That dict is the *only* thing the LLM
  sees — it decides the next step entirely from it. A tool that returns None is a bug.

Tools scan a directory of FASTQ files, read metadata, match read pairs, draft a
nf-core/rnaseq sample sheet, and validate the result.
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path
from typing import Any

Summary = dict[str, Any]

_FASTQ_EXTENSIONS = {".fastq", ".fq", ".fastq.gz", ".fq.gz"}

_PAIR_PATTERN = re.compile(
    r"^(?P<prefix>.+?)[_.]R?(?P<read>[12])(?P<tail>[_.].*?)?(?P<suffix>\.fastq(?:\.gz)?|\.fq(?:\.gz)?)$",
    re.IGNORECASE,
)

_ILLUMINA_SUFFIX = re.compile(r"_S\d+(?:_L\d+)?$")


def list_directory(path: str) -> Summary:
    """List the contents of a directory — files, subdirectories, and symlinks."""
    dirpath = Path(path)
    if not dirpath.is_dir():
        return {"error": "not_a_directory", "message": f"'{path}' is not a directory."}

    entries = []
    for p in sorted(dirpath.iterdir()):
        entries.append({"name": p.name, "type": "dir" if p.is_dir() else "file"})

    return {"path": str(dirpath), "n_entries": len(entries), "entries": entries}


def read_file(filepath: str, max_lines: int = 50) -> Summary:
    """Read the first N lines of a text file."""
    path = Path(filepath)
    if not path.is_file():
        return {"error": "not_a_file", "message": f"'{filepath}' does not exist."}

    try:
        lines = []
        with path.open() as f:
            for i, line in enumerate(f):
                if i >= max_lines:
                    break
                lines.append(line.rstrip("\n"))
        truncated = False
        with path.open() as f:
            for i, _ in enumerate(f):
                if i >= max_lines:
                    truncated = True
                    break
    except UnicodeDecodeError:
        return {"error": "binary_file", "message": f"'{filepath}' is not a text file."}

    return {
        "filepath": str(path),
        "n_lines": len(lines),
        "truncated": truncated,
        "content": "\n".join(lines),
    }


def scan_fastqs(directory: str) -> Summary:
    """Recursively scan a directory for FASTQ files.

    Returns a summary of discovered files: total count, list of filenames, any files
    that don't match expected naming patterns, and the directory scanned.
    """
    dirpath = Path(directory)
    if not dirpath.is_dir():
        return {"error": "not_a_directory", "message": f"'{directory}' is not a directory."}

    fastqs = []
    for p in sorted(dirpath.rglob("*")):
        if p.is_file() and any(p.name.endswith(ext) for ext in _FASTQ_EXTENSIONS):
            fastqs.append(str(p.resolve()))

    from core.session import SESSION
    SESSION.fastq_files = fastqs
    SESSION.fastq_pairs = {}

    return {
        "directory": str(dirpath),
        "n_files": len(fastqs),
        "files": fastqs,
    }


def read_metadata(filepath: str) -> Summary:
    """Read a metadata file (CSV or TSV) and return its contents as a summary.

    The metadata typically maps sample names to conditions, strandedness, or other
    attributes that the sample sheet needs. Returns column names and the first 10 rows.
    """
    path = Path(filepath)
    if not path.is_file():
        return {"error": "not_a_file", "message": f"'{filepath}' does not exist."}

    delimiter = "\t" if path.suffix in (".tsv", ".txt") else ","
    with path.open() as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        columns = reader.fieldnames or []
        rows = []
        for i, row in enumerate(reader):
            if i >= 10:
                break
            rows.append(dict(row))

    return {
        "filepath": str(path),
        "columns": list(columns),
        "n_rows_shown": len(rows),
        "rows": rows,
    }


def match_pairs() -> Summary:
    """Match the FASTQs from the last scan_fastqs into R1/R2 pairs by filename pattern.

    Every pair, incomplete pair and unpaired file gets a pair_id; the tools keep the
    paths, and draft_samplesheet refers to them by pair_id only.
    """
    from core.session import SESSION

    if not SESSION.fastq_files:
        return {"error": "no_fastqs", "message": "No FASTQs scanned. Call scan_fastqs first."}
    result = _match_pairs(SESSION.fastq_files)
    SESSION.fastq_pairs = {
        e["pair_id"]: {k: e[k] for k in ("fastq_1", "fastq_2") if k in e}
        for e in result["matched"] + result["incomplete"] + result["unpaired"]
    }
    return result


def _match_pairs(fastq_list: list[str]) -> Summary:
    """Pair FASTQ paths by filename: _R1/_R2 (or _1/_2) with a shared prefix.

    Pairs are keyed by the filename prefix (e.g. SRR26539594, or sampleA_S1_L001 for one
    lane); 'sample' is the prefix with Illumina _S<n>_L<nnn> removed, a suggested name.
    Unpaired files keep their filename as pair_id and are treated as single-end.
    """
    pairs: dict[str, dict[str, str]] = {}
    unpaired: list[dict[str, str]] = []

    for fq in fastq_list:
        name = Path(fq).name
        m = _PAIR_PATTERN.match(name)
        if not m:
            unpaired.append({"pair_id": name, "fastq_1": fq})
            continue

        prefix = m.group("prefix")
        sample = _ILLUMINA_SUFFIX.sub("", prefix)
        read_num = m.group("read")

        if prefix not in pairs:
            pairs[prefix] = {"sample": sample}
        key = f"fastq_{read_num}"
        if key in pairs[prefix]:
            unpaired.append({"pair_id": name, "fastq_1": fq})
        else:
            pairs[prefix][key] = fq

    matched = []
    incomplete = []
    for prefix, reads in sorted(pairs.items()):
        entry = {"pair_id": prefix, **reads}
        if "fastq_1" in reads and "fastq_2" in reads:
            matched.append(entry)
        else:
            incomplete.append(entry)

    return {
        "n_matched_pairs": len(matched),
        "n_incomplete": len(incomplete),
        "n_unpaired": len(unpaired),
        "matched": matched,
        "incomplete": incomplete,
        "unpaired": unpaired,
    }


def draft_samplesheet(samples: list[dict]) -> Summary:
    """Draft an nf-core/rnaseq sample sheet from pair IDs chosen by the agent.

    Each entry is {"pair_id", "sample", "strandedness"?}; the FASTQ paths are looked up
    from match_pairs, never taken from the LLM. Several pair_ids may share a sample name
    (lanes of one sample). Unknown or repeated pair_ids reject the whole draft. The CSV
    is stored as the session draft (resetting validation).
    """
    from core.session import SESSION

    if not SESSION.fastq_pairs:
        return {"error": "no_pairs", "message": "No matched FASTQs. Call scan_fastqs then match_pairs first."}

    errors: list[str] = []
    used: set[str] = set()
    rows = []
    for i, entry in enumerate(samples):
        pair_id = entry.get("pair_id", "")
        pair = SESSION.fastq_pairs.get(pair_id)
        if pair is None:
            errors.append(f"Entry {i}: unknown pair_id '{pair_id}'.")
            continue
        if pair_id in used:
            errors.append(f"Entry {i}: pair_id '{pair_id}' is used more than once.")
            continue
        used.add(pair_id)
        rows.append({
            "sample": entry.get("sample", ""),
            "pair_id": pair_id,
            "fastq_1": pair.get("fastq_1", ""),
            "fastq_2": pair.get("fastq_2", ""),
            "strandedness": entry.get("strandedness") or "auto",
        })
    if errors:
        return {"error": "invalid_samples", "errors": errors, "known_pair_ids": sorted(SESSION.fastq_pairs)}

    import io
    output = io.StringIO()
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(["sample", "fastq_1", "fastq_2", "strandedness"])
    for row in rows:
        writer.writerow([row["sample"], row["fastq_1"], row["fastq_2"], row["strandedness"]])

    SESSION.samplesheet_draft = output.getvalue()
    SESSION.samplesheet_draft_valid = False

    return {
        "n_rows": len(rows),
        "n_samples": len({r["sample"] for r in rows}),
        "rows": [{k: r[k] for k in ("sample", "pair_id", "strandedness")} for r in rows],
        "unused_pair_ids": sorted(set(SESSION.fastq_pairs) - used),
    }


def validate_samplesheet() -> Summary:
    """Validate the current draft (from draft_samplesheet) and record whether it passed."""
    from core.session import SESSION

    if SESSION.samplesheet_draft is None:
        return {"error": "no_draft", "message": "No draft sample sheet. Call draft_samplesheet first."}
    result = _check_samplesheet(SESSION.samplesheet_draft, check_files=True)
    SESSION.samplesheet_draft_valid = result["valid"]
    return result


def _check_samplesheet(sheet: str, check_files: bool = False) -> Summary:
    """Check a sample sheet CSV string.

    Checks: required columns present, sample names non-empty without whitespace, repeated
    sample names flagged, FASTQ paths non-empty (and existing, with check_files),
    strandedness is a valid value. Returns a list of warnings
    and errors, or confirms the sheet is valid.
    """
    lines = sheet.strip().splitlines()
    if not lines:
        return {"valid": False, "errors": ["Empty sample sheet."]}

    header = lines[0].split(",")
    required = {"sample", "fastq_1", "fastq_2", "strandedness"}
    missing = required - set(header)
    if missing:
        return {"valid": False, "errors": [f"Missing columns: {', '.join(sorted(missing))}"]}

    errors: list[str] = []
    warnings: list[str] = []
    samples_seen: set[str] = set()
    valid_strandedness = {"auto", "forward", "reverse", "unstranded"}

    reader = csv.DictReader(lines)
    for i, row in enumerate(reader, start=2):
        sample = row.get("sample", "").strip()
        if not sample:
            errors.append(f"Row {i}: empty sample name.")
        elif re.search(r"\s", sample):
            errors.append(f"Row {i}: sample name '{sample}' contains whitespace.")
        elif sample in samples_seen:
            warnings.append(f"Row {i}: repeated sample name '{sample}' (fine only if these are lanes of one sample).")
        samples_seen.add(sample)

        if not row.get("fastq_1", "").strip():
            errors.append(f"Row {i}: empty fastq_1.")
        if not row.get("fastq_2", "").strip():
            warnings.append(f"Row {i}: empty fastq_2 (single-end?).")
        if check_files:
            for key in ("fastq_1", "fastq_2"):
                path = row.get(key, "").strip()
                if path and not Path(path).is_file():
                    errors.append(f"Row {i}: {key} does not exist: {path}")

        strand = row.get("strandedness", "").strip()
        if strand and strand not in valid_strandedness:
            errors.append(f"Row {i}: invalid strandedness '{strand}'.")

    return {
        "valid": len(errors) == 0,
        "n_samples": len(samples_seen),
        "errors": errors,
        "warnings": warnings,
    }


def stage_fastqs(samples: list[dict]) -> Summary:
    """Symlink pairs to clean names ({sample}_R1.fastq.gz) in the run's staged_fastqs/.

    Only needed when original filenames aren't suitable. Entries are {"pair_id",
    "sample"}; the stored pairs are updated to the links so draft_samplesheet uses them.
    Writes rename_manifest.json recording every mapping.
    """
    from core.session import SESSION

    errors: list[str] = []
    planned: list[tuple[str, str, str, str]] = []  # pair_id, key, original, clean name
    for i, entry in enumerate(samples):
        pair_id, sample = entry.get("pair_id", ""), entry.get("sample", "")
        pair = SESSION.fastq_pairs.get(pair_id)
        if pair is None:
            errors.append(f"Entry {i}: unknown pair_id '{pair_id}'.")
            continue
        if not sample or re.search(r"[\s/]", sample):
            errors.append(f"Entry {i}: invalid sample name '{sample}'.")
            continue
        for key, original in pair.items():
            ext = next(e for e in (".fastq.gz", ".fq.gz", ".fastq", ".fq") if original.endswith(e))
            read = "R1" if key == "fastq_1" else "R2"
            planned.append((pair_id, key, original, f"{sample}_{read}{ext}"))
    names = [p[3] for p in planned]
    clashes = sorted({n for n in names if names.count(n) > 1})
    if clashes:
        errors.append(f"Clean names would collide: {', '.join(clashes)}.")
    if errors:
        return {"error": "invalid_samples", "errors": errors}

    stage = SESSION.require_paths().dir / "staged_fastqs"
    stage.mkdir(parents=True, exist_ok=True)
    manifest = []
    for pair_id, key, original, clean_name in planned:
        link = stage / clean_name
        if link.is_symlink():
            link.unlink()
        link.symlink_to(Path(original).resolve())
        SESSION.fastq_pairs[pair_id][key] = str(link)
        manifest.append({"pair_id": pair_id, "original": original, "staged": str(link), "clean_name": clean_name})

    manifest_path = stage / "rename_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    return {
        "n_links": len(manifest),
        "staging_dir": str(stage),
        "manifest_path": str(manifest_path),
        "staged": [{"pair_id": m["pair_id"], "clean_name": m["clean_name"]} for m in manifest],
    }


def write_report(report_markdown: str) -> Summary:
    """Write the agent's analysis report to the run directory.

    The agent supplies the narrative: what it found, decisions it made, warnings,
    and the final sample sheet summary. Called after save_samplesheet. If no sample
    sheet was saved, a warning banner is prepended so the report can't claim otherwise.
    """
    from core.session import SESSION

    paths = SESSION.require_paths()
    report_path = paths.dir / "report.md"
    result: Summary = {"report_path": str(report_path), "report_chars": len(report_markdown)}
    if not paths.samplesheet.is_file():
        warning = "No sample sheet was saved in this run; this report does not describe a usable sample sheet."
        report_markdown = f"> ⚠ {warning}\n\n{report_markdown}"
        result["warning"] = warning
    report_path.write_text(report_markdown.rstrip() + "\n")
    return result


def save_design(rows: list[dict]) -> Summary:
    """Write a design CSV (sample -> condition + optional covariates) to the run directory.

    Each row must have at least 'sample' and 'condition' keys. Additional columns
    are preserved. The downstream analysis agent reads this to know the experimental
    structure without the user re-specifying it.
    """
    from core.session import SESSION

    if not rows:
        return {"error": "empty_design", "message": "Design must have at least one row."}
    for i, row in enumerate(rows):
        if "sample" not in row:
            return {"error": "missing_sample", "message": f"Row {i}: missing 'sample' key."}
        if "condition" not in row:
            return {"error": "missing_condition", "message": f"Row {i}: missing 'condition' key."}

    import io
    paths = SESSION.require_paths()
    columns = list(rows[0].keys())
    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow(columns)
    for row in rows:
        writer.writerow([str(row.get(c, "")) for c in columns])
    paths.design.write_text(output.getvalue())

    conditions = sorted(set(row["condition"] for row in rows))
    return {
        "design_path": str(paths.design),
        "n_samples": len(rows),
        "conditions": conditions,
        "columns": columns,
    }


def save_samplesheet() -> Summary:
    """Write the validated draft to the run directory's canonical samplesheet path.

    Takes no arguments: the CSV is the exact draft that passed validate_samplesheet,
    so the LLM never retypes it.
    """
    from core.session import SESSION

    if SESSION.samplesheet_draft is None:
        return {"error": "no_draft", "message": "No draft sample sheet. Call draft_samplesheet first."}
    if not SESSION.samplesheet_draft_valid:
        return {"error": "not_validated",
                "message": "The current draft has not passed validate_samplesheet. Validate it (and fix any errors) first."}
    path = SESSION.require_paths().samplesheet
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(SESSION.samplesheet_draft)
    n_samples = len(SESSION.samplesheet_draft.strip().splitlines()) - 1
    return {"samplesheet_path": str(path), "n_samples": n_samples}


_RUN_ACCESSION = re.compile(r"[SED]RR\d+")


def sample_label_table(paths) -> str | None:
    """Code-built table for the approval screen: sample -> condition -> run -> GEO sample -> title.

    Runs are read from each row's fastq_1 filename (symlinks resolved) and labels from
    sample_metadata.csv (GEO/ENA), so the human checks names against GEO titles rather
    than file paths. Not a tool. None when there is no sample metadata.
    """
    if not paths.sample_metadata.is_file() or not paths.samplesheet.is_file():
        return None
    with paths.sample_metadata.open() as f:
        meta = {r["run_accession"]: r for r in csv.DictReader(f)}
    conditions: dict[str, str] = {}
    if paths.design.is_file():
        with paths.design.open() as f:
            conditions = {r.get("sample", ""): r.get("condition", "") for r in csv.DictReader(f)}

    header = ["sample", "condition", "run", "GEO sample", "GEO title"] if conditions else \
        ["sample", "run", "GEO sample", "GEO title"]
    rows, used, warnings = [], set(), []
    with paths.samplesheet.open() as f:
        for r in csv.DictReader(f):
            sample = r.get("sample", "")
            m = _RUN_ACCESSION.search(Path(r.get("fastq_1", "")).resolve().name)
            run = m.group(0) if m else ""
            info = meta.get(run)
            if info is None:
                warnings.append(f"⚠ {sample}: FASTQ does not match any downloaded run")
            else:
                used.add(run)
            row = [sample, run or "?", info["sample"] if info else "?", info["title"] if info else "?"]
            if conditions:
                row.insert(1, conditions.get(sample, "?"))
            rows.append(row)
    missing = sorted(set(meta) - used)
    if missing:
        warnings.append(f"⚠ Downloaded runs not in the sample sheet: {', '.join(missing)}")

    widths = [max(len(str(x)) for x in col) for col in zip(header, *rows)]
    lines = ["  ".join(str(x).ljust(w) for x, w in zip(line, widths)).rstrip() for line in [header, *rows]]
    lines.insert(1, "  ".join("-" * w for w in widths))
    return "\n".join(lines + warnings)

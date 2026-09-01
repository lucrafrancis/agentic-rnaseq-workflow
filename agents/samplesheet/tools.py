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
import re
from pathlib import Path
from typing import Any

Summary = dict[str, Any]

_FASTQ_EXTENSIONS = {".fastq", ".fq", ".fastq.gz", ".fq.gz"}

_PAIR_PATTERN = re.compile(
    r"^(?P<prefix>.+?)(?P<read>[_.]R?[12])(?P<suffix>\.fastq(?:\.gz)?|\.fq(?:\.gz)?)$",
    re.IGNORECASE,
)


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
            fastqs.append(str(p.relative_to(dirpath)))

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


def match_pairs(fastq_list: list[str]) -> Summary:
    """Match FASTQ files into R1/R2 pairs by filename pattern.

    Groups files by sample prefix, pairing _R1/_R2 (or _1/_2) suffixes. Reports
    matched pairs, unpaired files, and any ambiguous matches.
    """
    pairs: dict[str, dict[str, str]] = {}
    unpaired: list[str] = []

    for fq in fastq_list:
        name = Path(fq).name
        m = _PAIR_PATTERN.match(name)
        if not m:
            unpaired.append(fq)
            continue

        prefix = m.group("prefix")
        read_id = m.group("read")
        read_num = "1" if "1" in read_id else "2"

        if prefix not in pairs:
            pairs[prefix] = {}
        key = f"fastq_{read_num}"
        if key in pairs[prefix]:
            unpaired.append(fq)
        else:
            pairs[prefix][key] = fq

    matched = []
    incomplete = []
    for prefix, reads in sorted(pairs.items()):
        if "fastq_1" in reads and "fastq_2" in reads:
            matched.append({"sample": prefix, **reads})
        else:
            incomplete.append({"sample": prefix, **reads})

    return {
        "n_matched_pairs": len(matched),
        "n_incomplete": len(incomplete),
        "n_unpaired": len(unpaired),
        "matched": matched,
        "incomplete": incomplete,
        "unpaired": unpaired,
    }


def draft_samplesheet(matches: list[dict], metadata: dict | None = None) -> Summary:
    """Draft an nf-core/rnaseq sample sheet from matched pairs and optional metadata.

    Produces the CSV content with columns: sample, fastq_1, fastq_2, strandedness.
    Strandedness defaults to 'auto' unless metadata provides it. Returns a preview
    of the first 10 rows, the full CSV content, and summary stats.
    """
    rows = []
    for match in matches:
        sample = match.get("sample", "unknown")
        strandedness = "auto"
        if metadata and sample in metadata:
            strandedness = metadata[sample].get("strandedness", "auto")

        rows.append({
            "sample": sample,
            "fastq_1": match.get("fastq_1", ""),
            "fastq_2": match.get("fastq_2", ""),
            "strandedness": strandedness,
        })

    header = "sample,fastq_1,fastq_2,strandedness"
    lines = [header]
    for row in rows:
        lines.append(f"{row['sample']},{row['fastq_1']},{row['fastq_2']},{row['strandedness']}")
    csv_content = "\n".join(lines) + "\n"

    preview_lines = lines[:11]  # header + first 10 data rows

    return {
        "n_samples": len(rows),
        "columns": ["sample", "fastq_1", "fastq_2", "strandedness"],
        "preview": "\n".join(preview_lines),
        "csv_content": csv_content,
    }


def validate_samplesheet(sheet: str) -> Summary:
    """Validate a draft sample sheet CSV string.

    Checks: required columns present, no empty sample names, no duplicate sample names,
    FASTQ paths are non-empty, strandedness is a valid value. Returns a list of warnings
    and errors, or confirms the sheet is valid.
    """
    lines = sheet.strip().split("\n")
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
        elif sample in samples_seen:
            warnings.append(f"Row {i}: duplicate sample name '{sample}'.")
        samples_seen.add(sample)

        if not row.get("fastq_1", "").strip():
            errors.append(f"Row {i}: empty fastq_1.")
        if not row.get("fastq_2", "").strip():
            warnings.append(f"Row {i}: empty fastq_2 (single-end?).")

        strand = row.get("strandedness", "").strip()
        if strand and strand not in valid_strandedness:
            errors.append(f"Row {i}: invalid strandedness '{strand}'.")

    return {
        "valid": len(errors) == 0,
        "n_samples": len(samples_seen),
        "errors": errors,
        "warnings": warnings,
    }

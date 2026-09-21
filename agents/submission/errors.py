"""Parse nextflow logs and diagnose common nf-core/rnaseq failures.

Scans a nextflow log or stderr output for known error patterns and returns a structured
diagnosis the agent (or the user) can act on.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

_ERROR_PATTERNS: list[tuple[str, re.Pattern, str]] = [
    (
        "samplesheet_validation",
        re.compile(r"ERROR.*samplesheet.*validation", re.IGNORECASE),
        "The sample sheet failed nf-core validation. Check column names and file paths.",
    ),
    (
        "genome_not_found",
        re.compile(r"ERROR.*genome.*not found|Unknown genome", re.IGNORECASE),
        "The specified genome is not available. Check --genome value against igenomes.config.",
    ),
    (
        "missing_files",
        re.compile(r"No such file|FileNotFoundException|does not exist", re.IGNORECASE),
        "One or more input files could not be found. Verify FASTQ paths are absolute.",
    ),
    (
        "out_of_memory",
        re.compile(r"Process.*exceed(?:s|ed).*memory|OutOfMemoryError|oom-kill", re.IGNORECASE),
        "A process ran out of memory. Consider increasing memory limits or using a retry strategy.",
    ),
    (
        "container_error",
        re.compile(r"Container.*error|docker.*pull.*failed|singularity.*error", re.IGNORECASE),
        "Container pull or execution failed. Check Docker/Singularity is available and has network access.",
    ),
    (
        "pipeline_timeout",
        re.compile(r"Time limit exceeded|wall.?time", re.IGNORECASE),
        "A process hit its time limit. Consider increasing the time allocation.",
    ),
]


@dataclass
class Diagnosis:
    """A diagnosed error from a nextflow log."""

    category: str
    line: str
    suggestion: str


def diagnose_log(log_path: str) -> list[Diagnosis]:
    """Scan a nextflow log file for known error patterns.

    Returns a list of Diagnosis objects, one per matched pattern. Returns an empty list
    if no known patterns are found (the log may still contain errors — just not ones
    we recognise).
    """
    path = Path(log_path)
    if not path.is_file():
        return []

    text = path.read_text()
    diagnoses = []
    for category, pattern, suggestion in _ERROR_PATTERNS:
        for match in pattern.finditer(text):
            line_start = text.rfind("\n", 0, match.start()) + 1
            line_end = text.find("\n", match.end())
            line = text[line_start:line_end].strip()
            diagnoses.append(Diagnosis(category=category, line=line, suggestion=suggestion))

    return diagnoses


def format_diagnoses(diagnoses: list[Diagnosis]) -> str:
    """Format a list of diagnoses into a human-readable summary."""
    if not diagnoses:
        return "No known error patterns found in the log."

    lines = [f"Found {len(diagnoses)} issue(s):\n"]
    for i, d in enumerate(diagnoses, 1):
        lines.append(f"  {i}. [{d.category}] {d.line}")
        lines.append(f"     Suggestion: {d.suggestion}")
    return "\n".join(lines)

"""JSON tool definitions sent to Claude, plus the name -> function dispatch table.

The schemas here MUST stay in sync with the signatures in tools.py: the `input_schema`
describes the arguments Claude is allowed to pass, and TOOL_FUNCTIONS maps the tool name
back to the Python function the loop actually calls.

Kept deliberately in its own module so the loop stays generic — it never hard-codes tool
names, it just iterates over whatever is registered here.
"""

from __future__ import annotations

from typing import Any, Callable

from agents.samplesheet import tools

TOOL_SCHEMAS: list[dict[str, Any]] = [
    {
        "name": "scan_fastqs",
        "description": "Recursively scan a directory for FASTQ files. Returns a list of "
        "discovered filenames, the total count, and the directory scanned. Call this first "
        "to discover the input data.",
        "input_schema": {
            "type": "object",
            "properties": {
                "directory": {
                    "type": "string",
                    "description": "Path to the directory containing FASTQ files.",
                },
            },
            "required": ["directory"],
        },
    },
    {
        "name": "read_metadata",
        "description": "Read a metadata file (CSV or TSV) and return its column names and "
        "first 10 rows. Use this to discover sample attributes like strandedness or "
        "condition that belong in the sample sheet.",
        "input_schema": {
            "type": "object",
            "properties": {
                "filepath": {
                    "type": "string",
                    "description": "Path to the metadata CSV/TSV file.",
                },
            },
            "required": ["filepath"],
        },
    },
    {
        "name": "match_pairs",
        "description": "Match FASTQ files into R1/R2 pairs by filename pattern. Returns "
        "matched pairs, incomplete pairs, and unpaired files. Run after scan_fastqs.",
        "input_schema": {
            "type": "object",
            "properties": {
                "fastq_list": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "List of FASTQ file paths from scan_fastqs.",
                },
            },
            "required": ["fastq_list"],
        },
    },
    {
        "name": "draft_samplesheet",
        "description": "Draft an nf-core/rnaseq sample sheet CSV from matched pairs and "
        "optional metadata. Returns a preview (first 10 rows), full CSV content, and "
        "summary stats. Strandedness defaults to 'auto' unless metadata overrides it.",
        "input_schema": {
            "type": "object",
            "properties": {
                "matches": {
                    "type": "array",
                    "items": {
                        "type": "object",
                        "properties": {
                            "sample": {"type": "string"},
                            "fastq_1": {"type": "string"},
                            "fastq_2": {"type": "string"},
                        },
                    },
                    "description": "Matched pairs from match_pairs.",
                },
                "metadata": {
                    "type": "object",
                    "description": "Optional: sample name -> attributes dict with strandedness overrides.",
                },
            },
            "required": ["matches"],
        },
    },
    {
        "name": "validate_samplesheet",
        "description": "Validate a draft sample sheet CSV string. Checks required columns, "
        "empty fields, duplicates, and strandedness values. Returns errors and warnings.",
        "input_schema": {
            "type": "object",
            "properties": {
                "sheet": {
                    "type": "string",
                    "description": "The sample sheet CSV content to validate.",
                },
            },
            "required": ["sheet"],
        },
    },
]

TOOL_FUNCTIONS: dict[str, Callable[..., dict[str, Any]]] = {
    "scan_fastqs": tools.scan_fastqs,
    "read_metadata": tools.read_metadata,
    "match_pairs": tools.match_pairs,
    "draft_samplesheet": tools.draft_samplesheet,
    "validate_samplesheet": tools.validate_samplesheet,
}

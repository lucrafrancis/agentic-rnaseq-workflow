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
        "name": "list_directory",
        "description": "List the contents of a directory — files and subdirectories. "
        "Use to explore the data folder structure before scanning for FASTQs.",
        "input_schema": {
            "type": "object",
            "properties": {
                "path": {"type": "string", "description": "Directory path to list."},
            },
            "required": ["path"],
        },
    },
    {
        "name": "read_file",
        "description": "Read the first N lines of a text file. Use to inspect READMEs, "
        "config files, FASTQ headers, or anything that helps understand the data.",
        "input_schema": {
            "type": "object",
            "properties": {
                "filepath": {"type": "string", "description": "Path to the file."},
                "max_lines": {"type": "integer", "description": "Lines to read (default 50)."},
            },
            "required": ["filepath"],
        },
    },
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
    {
        "name": "stage_fastqs",
        "description": "Create symlinks with clean filenames in a staging directory. Only "
        "needed when original names aren't suitable (e.g. Illumina index/lane segments). "
        "Writes a rename_manifest.json recording every mapping. Returns staged pairs with "
        "updated paths for the sample sheet.",
        "input_schema": {
            "type": "object",
            "properties": {
                "pairs": {
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
                "source_dir": {
                    "type": "string",
                    "description": "Directory the original FASTQ paths are relative to.",
                },
                "staging_dir": {
                    "type": "string",
                    "description": "Directory to create symlinks in.",
                },
            },
            "required": ["pairs", "source_dir", "staging_dir"],
        },
    },
    {
        "name": "save_samplesheet",
        "description": "Write the validated sample sheet CSV to disk. Call this after "
        "validate_samplesheet confirms no errors, then call write_report.",
        "input_schema": {
            "type": "object",
            "properties": {
                "csv_content": {
                    "type": "string",
                    "description": "The validated sample sheet CSV content.",
                },
                "output_path": {
                    "type": "string",
                    "description": "Where to write the CSV file.",
                },
            },
            "required": ["csv_content", "output_path"],
        },
    },
    {
        "name": "save_design",
        "description": "Write a design CSV mapping samples to experimental conditions. "
        "Each row needs at least 'sample' and 'condition'. Call after save_samplesheet "
        "if you can infer the experimental design from sample names or metadata.",
        "input_schema": {
            "type": "object",
            "properties": {
                "rows": {
                    "type": "array",
                    "items": {
                        "type": "object",
                        "properties": {
                            "sample": {"type": "string"},
                            "condition": {"type": "string"},
                        },
                        "required": ["sample", "condition"],
                    },
                    "description": "List of sample-to-condition mappings.",
                },
            },
            "required": ["rows"],
        },
    },
    {
        "name": "write_report",
        "description": "Write a Markdown report summarising what was found, decisions made, "
        "and the final sample sheet. Call this last, after save_samplesheet.",
        "input_schema": {
            "type": "object",
            "properties": {
                "report_markdown": {
                    "type": "string",
                    "description": "The full report as Markdown.",
                },
            },
            "required": ["report_markdown"],
        },
    },
]

TOOL_FUNCTIONS: dict[str, Callable[..., dict[str, Any]]] = {
    "list_directory": tools.list_directory,
    "read_file": tools.read_file,
    "scan_fastqs": tools.scan_fastqs,
    "read_metadata": tools.read_metadata,
    "match_pairs": tools.match_pairs,
    "draft_samplesheet": tools.draft_samplesheet,
    "validate_samplesheet": tools.validate_samplesheet,
    "stage_fastqs": tools.stage_fastqs,
    "save_samplesheet": tools.save_samplesheet,
    "save_design": tools.save_design,
    "write_report": tools.write_report,
}

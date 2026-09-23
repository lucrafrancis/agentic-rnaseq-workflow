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
        "description": "Match the FASTQs from the last scan_fastqs into R1/R2 pairs by filename "
        "pattern. Takes no arguments. Every matched pair, incomplete pair and unpaired file gets "
        "a pair_id (the filename prefix, e.g. SRR26539594, or the filename for unpaired files) "
        "with a suggested sample name. Refer to FASTQs by pair_id from here on.",
        "input_schema": {"type": "object", "properties": {}},
    },
    {
        "name": "draft_samplesheet",
        "description": "Draft the nf-core/rnaseq sample sheet. Give one entry per pair_id to "
        "include, with the sample name you chose; the tools fill in the FASTQ paths. Lanes of "
        "one sample: several pair_ids with the same sample name. Unknown or repeated pair_ids "
        "reject the draft. Returns the rows and any unused pair_ids. The draft is kept by the "
        "tools; validate_samplesheet and save_samplesheet act on it, and drafting again "
        "replaces it and resets validation.",
        "input_schema": {
            "type": "object",
            "properties": {
                "samples": {
                    "type": "array",
                    "items": {
                        "type": "object",
                        "properties": {
                            "pair_id": {"type": "string", "description": "A pair_id from match_pairs."},
                            "sample": {"type": "string", "description": "Sample name (no whitespace)."},
                            "strandedness": {
                                "type": "string",
                                "enum": ["auto", "forward", "reverse", "unstranded"],
                                "description": "Defaults to 'auto'.",
                            },
                        },
                        "required": ["pair_id", "sample"],
                    },
                },
            },
            "required": ["samples"],
        },
    },
    {
        "name": "validate_samplesheet",
        "description": "Validate the current draft from draft_samplesheet. Checks required "
        "columns, sample names, repeated names, that every FASTQ exists, and strandedness "
        "values. Returns errors and warnings.",
        "input_schema": {"type": "object", "properties": {}},
    },
    {
        "name": "stage_fastqs",
        "description": "Symlink pairs to clean filenames ({sample}_R1.fastq.gz) in the run "
        "directory. Only needed when original names aren't suitable (e.g. Illumina index/lane "
        "segments); nf-core uses the sample column for naming, so this is rarely necessary. "
        "The stored pairs then point at the links. Writes rename_manifest.json.",
        "input_schema": {
            "type": "object",
            "properties": {
                "samples": {
                    "type": "array",
                    "items": {
                        "type": "object",
                        "properties": {
                            "pair_id": {"type": "string"},
                            "sample": {"type": "string"},
                        },
                        "required": ["pair_id", "sample"],
                    },
                },
            },
            "required": ["samples"],
        },
    },
    {
        "name": "save_samplesheet",
        "description": "Write the current draft to the run directory. Takes no arguments; "
        "refuses unless validate_samplesheet passed on the current draft. Then call write_report.",
        "input_schema": {"type": "object", "properties": {}},
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

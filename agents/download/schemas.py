"""JSON tool definitions sent to Claude, plus the name -> function dispatch table.

Kept in its own module so the loop stays generic.
"""

from __future__ import annotations

from typing import Any, Callable

from agents.download import tools

TOOL_SCHEMAS: list[dict[str, Any]] = [
    {
        "name": "resolve_accession",
        "description": (
            "Resolve a GEO/SRA accession to download metadata. Supports GSE (GEO series), "
            "SRP/ERP/DRP (SRA study), PRJNA (BioProject), and SRR (individual run). "
            "Queries NCBI and ENA APIs. Saves full metadata (URLs, MD5 checksums) to disk — "
            "downstream tools read it directly so checksums never pass through you. "
            "Call this first."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "accession": {
                    "type": "string",
                    "description": "The GEO or SRA accession (e.g. GSE157852, SRP282091, SRR12626034).",
                },
            },
            "required": ["accession"],
        },
    },
    {
        "name": "check_existing_files",
        "description": (
            "Check if expected FASTQ files already exist in a directory and validate their "
            "MD5 checksums against ENA metadata. Call after resolve_accession to see what's "
            "already downloaded. Returns status for each file: valid, missing, or corrupted."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "search_dir": {
                    "type": "string",
                    "description": "Directory to check for existing FASTQ files.",
                },
            },
            "required": ["search_dir"],
        },
    },
    {
        "name": "generate_download_script",
        "description": (
            "Generate a shell script to download FASTQ files from ENA. Automatically skips "
            "files that already exist and pass MD5 validation in the output directory. "
            "Writes download.sh to the run directory. The script will be presented to the "
            "user for approval before execution. Call after check_existing_files shows "
            "missing or corrupted files."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "output_dir": {
                    "type": "string",
                    "description": "Directory to download FASTQ files into.",
                },
            },
            "required": ["output_dir"],
        },
    },
]

TOOL_FUNCTIONS: dict[str, Callable[..., dict[str, Any]]] = {
    "resolve_accession": tools.resolve_accession,
    "check_existing_files": tools.check_existing_files,
    "generate_download_script": tools.generate_download_script,
}

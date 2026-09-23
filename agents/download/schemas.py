"""JSON tool definitions sent to Claude, plus the name -> function dispatch table.

Kept in its own module so the loop stays generic.
"""

from __future__ import annotations

from typing import Any, Callable

from agents.download import counts, tools

TOOL_SCHEMAS: list[dict[str, Any]] = [
    {
        "name": "resolve_accession",
        "description": (
            "Resolve a GEO/SRA accession to download metadata. Supports GSE (GEO series), "
            "SRP/ERP/DRP (SRA study), PRJNA (BioProject), and SRR (individual run). "
            "Queries NCBI and ENA APIs. Returns each run with its GEO sample (GSM) and title. "
            "Saves full metadata (URLs, MD5 checksums) to disk — "
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
                "runs": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": (
                        "Optional run accessions (SRR...) to restrict to. Omit for all runs. "
                        "Use when the prompt asks for a subset of samples."
                    ),
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
                "runs": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": (
                        "Optional run accessions (SRR...) to restrict to. Omit for all runs. "
                        "Use when the prompt asks for a subset of samples."
                    ),
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


# --- Counts mode: fetch a processed count matrix from GEO instead of FASTQs ---

COUNTS_TOOL_SCHEMAS: list[dict[str, Any]] = [
    {
        "name": "list_geo_count_sources",
        "description": (
            "List candidate count files for a GEO series (author supplementary files and "
            "NCBI-generated raw counts) plus every sample's title and characteristics. "
            "URLs are saved to disk — refer to files by name. Call this first."
        ),
        "input_schema": {
            "type": "object",
            "properties": {"accession": {"type": "string", "description": "GEO series, e.g. GSE164073."}},
            "required": ["accession"],
        },
    },
    {
        "name": "preview_geo_file",
        "description": (
            "Download (cached) and preview a listed file: detected format, every column with "
            "numeric/integer flags and any exact GSM/title match, and the first rows. Use it "
            "to choose the gene ID column and map sample columns to GSMs."
        ),
        "input_schema": {
            "type": "object",
            "properties": {"filename": {"type": "string", "description": "A file name from list_geo_count_sources."}},
            "required": ["filename"],
        },
    },
    {
        "name": "fetch_geo_counts",
        "description": (
            "Parse a listed file into the run's count matrix. Map every sample column to its "
            "GSM (one column per GSM); unmapped columns are dropped (annotation columns, or "
            "samples you intend to exclude). Validates numeric, non-negative values, sums "
            "duplicate gene IDs, renames samples to their GEO titles and adds gene symbols."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "filename": {"type": "string"},
                "gene_id_column": {"type": "string", "description": "Column holding gene IDs (exact name from the preview)."},
                "sample_map": {
                    "type": "object",
                    "additionalProperties": {"type": "string"},
                    "description": "{file column: GSM ID}, e.g. {\"MW1_cornea_mock_1\": \"GSM4996084\"}.",
                },
            },
            "required": ["filename", "gene_id_column", "sample_map"],
        },
    },
    {
        "name": "save_geo_design",
        "description": (
            "Write design.csv for the fetched samples. Give exactly one of condition_field "
            "(a characteristics key such as \"treatment\") or conditions ({GSM: label}) when "
            "the condition is only encoded in sample titles. Other characteristics that vary "
            "across samples are kept automatically as covariate columns."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "condition_field": {"type": "string"},
                "conditions": {"type": "object", "additionalProperties": {"type": "string"}},
            },
        },
    },
]

COUNTS_TOOL_FUNCTIONS: dict[str, Callable[..., dict[str, Any]]] = {
    "list_geo_count_sources": counts.list_geo_count_sources,
    "preview_geo_file": counts.preview_geo_file,
    "fetch_geo_counts": counts.fetch_geo_counts,
    "save_geo_design": counts.save_geo_design,
}

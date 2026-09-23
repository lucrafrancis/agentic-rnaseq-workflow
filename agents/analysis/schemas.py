"""JSON tool definitions sent to Claude, plus the name -> function dispatch table.

The schemas here MUST stay in sync with the signatures in tools.py: the `input_schema`
describes the arguments Claude is allowed to pass, and TOOL_FUNCTIONS maps the tool name
back to the Python function the loop actually calls.
"""

from __future__ import annotations

from typing import Any, Callable

from agents.analysis import tools

TOOL_SCHEMAS: list[dict[str, Any]] = [
    {
        "name": "fetch_geo_metadata",
        "description": "Fetch GEO series metadata via NCBI E-utilities. Returns title, "
        "summary, organism, cell type, sample descriptions, and PubMed IDs. "
        "Call this early if the prompt mentions a GEO accession — it provides "
        "essential context (cell type, experimental design, associated paper).",
        "input_schema": {
            "type": "object",
            "properties": {
                "accession": {
                    "type": "string",
                    "description": "GEO series accession (e.g. GSE245856).",
                },
            },
            "required": ["accession"],
        },
    },
    {
        "name": "fetch_abstract",
        "description": "Fetch a PubMed abstract via NCBI E-utilities. Returns title, "
        "authors, journal, year, and abstract text. Use this to get citable context "
        "for the report — cite as (Author et al., Year; PMID: <id>).",
        "input_schema": {
            "type": "object",
            "properties": {
                "pmid": {
                    "type": "string",
                    "description": "PubMed ID (numeric string).",
                },
            },
            "required": ["pmid"],
        },
    },
    {
        "name": "scan_results",
        "description": "Scan a directory for analysis inputs: count matrices, design files, "
        "and MultiQC data. Detects both nf-core/rnaseq output structure and loose "
        "user-provided files. Lists all tabular files found. Call this first.",
        "input_schema": {
            "type": "object",
            "properties": {
                "results_dir": {
                    "type": "string",
                    "description": "Path to the directory containing analysis inputs.",
                },
            },
            "required": ["results_dir"],
        },
    },
    {
        "name": "load_counts",
        "description": "Load a gene count matrix and optional design file. Auto-detects "
        "format: nf-core TSV (gene_id/gene_name), featureCounts (Geneid/Chr/Start/...), "
        "or generic CSV/TSV (first column = gene IDs, rest = samples). Also auto-detects "
        "separator. Reports sample names, library sizes, and conditions if design loaded.",
        "input_schema": {
            "type": "object",
            "properties": {
                "counts_path": {
                    "type": "string",
                    "description": "Path to the count matrix file (CSV or TSV).",
                },
                "design_path": {
                    "type": "string",
                    "description": "Optional path to a design file (CSV/TSV) with 'sample' and 'condition' columns.",
                },
            },
            "required": ["counts_path"],
        },
    },
    {
        "name": "inspect_counts",
        "description": "Summarise the loaded count matrix: value range, fraction of "
        "non-integer values, fraction of zeros, per-sample stats. Use this to assess "
        "whether data is raw counts, TPM/FPKM, or log-transformed before running DE. "
        "DESeq2 requires raw counts — do not run it on normalised data. "
        "Requires load_counts first.",
        "input_schema": {"type": "object", "properties": {}},
    },
    {
        "name": "set_design",
        "description": "Define the sample-to-condition mapping when no design.csv exists. "
        "Each row needs 'sample' and 'condition'. Validates sample names match the loaded "
        "count matrix. Call after load_counts if no design was loaded.",
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
        "name": "compute_qc",
        "description": "Compute QC metrics: library sizes, gene detection rates, PCA on "
        "log2(counts+1) for outlier detection. Parses MultiQC stats if available. "
        "Requires load_counts first.",
        "input_schema": {"type": "object", "properties": {}, "required": []},
    },
    {
        "name": "filter_low_counts",
        "description": "Drop genes where fewer than min_samples have at least min_count "
        "reads. Run before DE so multiple-testing correction is not diluted. "
        "Requires load_counts first.",
        "input_schema": {
            "type": "object",
            "properties": {
                "min_count": {
                    "type": "integer",
                    "description": "Minimum read count threshold (default 10).",
                },
                "min_samples": {
                    "type": "integer",
                    "description": "Minimum number of samples meeting the threshold (default 2).",
                },
            },
            "required": [],
        },
    },
    {
        "name": "run_deseq2",
        "description": "Run differential expression with PyDESeq2 for a given contrast. "
        "The contrast is [factor, test, reference], e.g. ['condition', 'treated', 'control']. "
        "Optional covariates (other design columns such as batch, donor, tissue) are adjusted "
        "for: design = ~ covariates + factor. Refuses confounded designs. "
        "Requires load_counts and a design (set_design or loaded via load_counts).",
        "input_schema": {
            "type": "object",
            "properties": {
                "contrast": {
                    "type": "array",
                    "items": {"type": "string"},
                    "minItems": 3,
                    "maxItems": 3,
                    "description": "DE contrast as [factor, test_level, reference_level].",
                },
                "covariates": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Design columns to adjust for (categorical), e.g. ['donor'] or ['batch'].",
                },
            },
            "required": ["contrast"],
        },
    },
    {
        "name": "get_top_genes",
        "description": "Return top DE genes by adjusted p-value for inspection. "
        "Requires run_deseq2 first.",
        "input_schema": {
            "type": "object",
            "properties": {
                "n": {
                    "type": "integer",
                    "description": "Number of top genes to return (default 20).",
                },
                "direction": {
                    "type": "string",
                    "enum": ["up", "down", "both"],
                    "description": "Filter by direction: 'up', 'down', or 'both' (default).",
                },
            },
            "required": [],
        },
    },
    {
        "name": "run_enrichment",
        "description": "Run gene set enrichment via Enrichr (GO and KEGG). Pass gene "
        "symbols (gene_name), not Ensembl IDs. Requires a non-empty gene list.",
        "input_schema": {
            "type": "object",
            "properties": {
                "gene_list": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Gene symbols to test for enrichment.",
                },
                "organism": {
                    "type": "string",
                    "description": "Organism: 'human', 'mouse', or 'yeast' (default 'human').",
                },
                "label": {
                    "type": "string",
                    "description": "Label for this enrichment run, e.g. 'upregulated' or "
                    "'downregulated'. Results accumulate across calls — each label gets "
                    "its own figures (enrichment_<label>_<gene_set>.png).",
                },
            },
            "required": ["gene_list"],
        },
    },
    {
        "name": "summarize_findings",
        "description": "Consolidate the analysis state (QC, DE, enrichment) into one "
        "factual summary including data source, software versions, and all results. "
        "Use the returned software_versions for the Methods section — do not guess. "
        "Non-mutating.",
        "input_schema": {"type": "object", "properties": {}, "required": []},
    },
    {
        "name": "generate_figures",
        "description": "Draw all figures for the current analysis and return each figure's "
        "path with a factual caption (e.g. how many genes the heatmap shows). Call after "
        "summarize_findings and BEFORE writing the report, so the report describes the real "
        "figures. Optional pca_color_by chooses which design columns get their own PCA plot "
        "(default: all informative columns; identifier and constant columns are never allowed).",
        "input_schema": {
            "type": "object",
            "properties": {
                "pca_color_by": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Design columns to colour extra PCA plots by, e.g. ['tissue'].",
                },
            },
        },
    },
    {
        "name": "write_report",
        "description": "Write the final Markdown report. Link figures only by the paths "
        "returned by generate_figures — links to other figures are rejected and nothing is "
        "written. The tool inserts each figure's file path and factual caption beneath it and "
        "appends the standard disclaimer. Call last.",
        "input_schema": {
            "type": "object",
            "properties": {
                "report_markdown": {
                    "type": "string",
                    "description": "The full analysis report as Markdown with inline figure references.",
                },
            },
            "required": ["report_markdown"],
        },
    },
]

TOOL_FUNCTIONS: dict[str, Callable[..., dict[str, Any]]] = {
    "fetch_geo_metadata": tools.fetch_geo_metadata,
    "fetch_abstract": tools.fetch_abstract,
    "scan_results": tools.scan_results,
    "load_counts": tools.load_counts,
    "inspect_counts": tools.inspect_counts,
    "set_design": tools.set_design,
    "compute_qc": tools.compute_qc,
    "filter_low_counts": tools.filter_low_counts,
    "run_deseq2": tools.run_deseq2,
    "get_top_genes": tools.get_top_genes,
    "run_enrichment": tools.run_enrichment,
    "summarize_findings": tools.summarize_findings,
    "generate_figures": tools.generate_figures,
    "write_report": tools.write_report,
}

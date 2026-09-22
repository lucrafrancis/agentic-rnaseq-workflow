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
            },
            "required": ["gene_list"],
        },
    },
    {
        "name": "summarize_findings",
        "description": "Consolidate the analysis state (QC, DE, enrichment) into one "
        "factual summary to write the report from. Non-mutating.",
        "input_schema": {"type": "object", "properties": {}, "required": []},
    },
    {
        "name": "generate_report",
        "description": "Generate all analysis figures and write the Markdown report. "
        "Figures are generated at known paths under figures/: pca.png, library_sizes.png, "
        "sample_correlation.png, volcano.png, ma_plot.png, de_heatmap.png, "
        "enrichment_*.png, pca_<variable>.png, pc_association.png. "
        "Place figure references inline in the relevant report sections using "
        "![caption](figures/<name>.png). The tool writes the report as-is. Call this last.",
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
    "generate_report": tools.generate_report,
}

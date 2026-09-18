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
        "description": "Scan an nf-core/rnaseq output directory for count matrices, "
        "TPM values, and MultiQC data. Call this first to discover what the pipeline "
        "produced. Also checks for a design.csv from the samplesheet agent.",
        "input_schema": {
            "type": "object",
            "properties": {
                "results_dir": {
                    "type": "string",
                    "description": "Path to the nf-core/rnaseq results directory.",
                },
            },
            "required": ["results_dir"],
        },
    },
    {
        "name": "load_counts",
        "description": "Load the gene count matrix TSV and optional design CSV. Reports "
        "sample names, gene counts, library sizes, and conditions if a design is loaded. "
        "Also loads TPM values if available in the same directory.",
        "input_schema": {
            "type": "object",
            "properties": {
                "counts_path": {
                    "type": "string",
                    "description": "Path to the salmon.merged.gene_counts.tsv file.",
                },
                "design_path": {
                    "type": "string",
                    "description": "Optional path to a design CSV with sample and condition columns.",
                },
            },
            "required": ["counts_path"],
        },
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
        "description": "Render figures (volcano, MA, PCA, enrichment) and assemble the "
        "final Markdown report. Provide the report narrative as report_markdown; the "
        "tool adds the figures. Call this last.",
        "input_schema": {
            "type": "object",
            "properties": {
                "report_markdown": {
                    "type": "string",
                    "description": "The full analysis report as Markdown, written by you from the findings.",
                },
            },
            "required": ["report_markdown"],
        },
    },
]

TOOL_FUNCTIONS: dict[str, Callable[..., dict[str, Any]]] = {
    "scan_results": tools.scan_results,
    "load_counts": tools.load_counts,
    "set_design": tools.set_design,
    "compute_qc": tools.compute_qc,
    "filter_low_counts": tools.filter_low_counts,
    "run_deseq2": tools.run_deseq2,
    "get_top_genes": tools.get_top_genes,
    "run_enrichment": tools.run_enrichment,
    "summarize_findings": tools.summarize_findings,
    "generate_report": tools.generate_report,
}

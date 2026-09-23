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
        "log2(counts+1) for outlier detection. Requires load_counts first.",
        "input_schema": {"type": "object", "properties": {}, "required": []},
    },
    {
        "name": "read_multiqc",
        "description": "Read MultiQC general statistics from the nf-core results (e.g. mapping "
        "rate, duplication, GC, trimming) — whichever columns this pipeline version reports. "
        "Use when scan_results lists multiqc_files. Requires scan_results; call after "
        "load_counts so rows are matched to the analysis samples.",
        "input_schema": {
            "type": "object",
            "properties": {
                "path": {
                    "type": "string",
                    "description": "One of scan_results' multiqc_files. Optional when there is only one.",
                },
            },
            "required": [],
        },
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
        "description": "Over-representation analysis via Enrichr (GO Biological Process + KEGG) "
        "for one direction. The tool selects the genes itself from the DESeq2 results — you never "
        "pass gene names: padj < padj_max AND |log2FC| >= lfc_min, split by direction, ranked by "
        "padj then |log2FC|, capped at max_genes, converted to gene symbols. Call once with "
        "direction='up' and once with direction='down'. Saves the exact input genes and "
        "Enrichr's raw results to analysis/.",
        "input_schema": {
            "type": "object",
            "properties": {
                "direction": {"type": "string", "enum": ["up", "down"]},
                "padj_max": {"type": "number", "description": "Adjusted p-value cutoff (default 0.05)."},
                "lfc_min": {
                    "type": "number",
                    "description": "Minimum |log2FC| (default 1.0 = 2-fold). Lower (e.g. 0.585 = "
                    "1.5-fold) only if too few genes pass; say so in the report.",
                },
                "max_genes": {
                    "type": "integer",
                    "description": "Cap on genes per direction (default 500), taken in rank order.",
                },
                "organism": {"type": "string", "enum": ["human", "mouse", "yeast"]},
            },
            "required": ["direction"],
        },
    },
    {
        "name": "summarize_findings",
        "description": "Consolidate the analysis state into one factual summary. Returns "
        "'facts' (placeholder names with their current values, for {{name}} in the report), "
        "'tables_available' and 'tables_required' (for {{table:name}}), 'references' "
        "(for {{cite:PMID}}), plus provenance, versions and results. Non-mutating.",
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
        "description": "Write the final Markdown report. Never type a number, statistic, "
        "parameter, version or citation — use placeholders: {{fact.name}} (from "
        "summarize_findings facts), {{gene:SYMBOL}} (renders log2FC and padj), {{cite:PMID}} "
        "(fetched references only), {{table:name}} (code-generated tables; every table in "
        "tables_required must be placed). Place EVERY figure from generate_figures as "
        "![description](figures/<name>.png) — no other figures. "
        "Any problem rejects the whole report with a list of all issues — fix them all and "
        "resubmit. The tool adds figure captions and the standard disclaimer. Call last.",
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
    "read_multiqc": tools.read_multiqc,
    "filter_low_counts": tools.filter_low_counts,
    "run_deseq2": tools.run_deseq2,
    "get_top_genes": tools.get_top_genes,
    "run_enrichment": tools.run_enrichment,
    "summarize_findings": tools.summarize_findings,
    "generate_figures": tools.generate_figures,
    "write_report": tools.write_report,
}

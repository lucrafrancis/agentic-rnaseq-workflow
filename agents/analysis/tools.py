"""Downstream analysis tools for bulk RNA-seq.

THE CONTRACT:
  Every tool is a deterministic Python function that does real work AND returns a
  structured, JSON-serializable summary dict. That dict is the *only* thing the LLM
  sees — it decides the next step entirely from it. A tool that returns None is a bug.

Tools operate on DataFrames held by the running session (see core/session.py). Each
tool validates the state it expects and returns a clear error summary the agent can
react to.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from core.session import SESSION

Summary = dict[str, Any]

_QUANTILES = {"min": 0.0, "p25": 0.25, "median": 0.5, "p75": 0.75, "p95": 0.95, "max": 1.0}


def _quantile_summary(values) -> dict[str, float]:
    arr = np.asarray(values, dtype=float)
    return {name: round(float(np.quantile(arr, q)), 4) for name, q in _QUANTILES.items()}


_COUNT_KEYWORDS = {"count", "counts", "gene_counts", "raw_counts", "featurecounts"}
_DESIGN_KEYWORDS = {"design", "metadata", "coldata", "sample_info", "conditions", "phenotype"}
_TABULAR_EXTENSIONS = {".csv", ".tsv", ".txt"}


def scan_results(results_dir: str) -> Summary:
    """Scan a directory for analysis inputs: count matrices, design files, MultiQC.

    Detects both nf-core/rnaseq output structure and loose user-provided files.
    Non-mutating; the agent calls this first to orient itself.
    """
    rdir = Path(results_dir)
    if not rdir.is_dir():
        return {"error": "not_a_directory", "message": f"'{results_dir}' is not a directory."}

    counts_path = None
    tpm_path = None
    mqc_path = None
    design_path = None
    source = "unknown"

    # --- Try nf-core structure first ---
    star_salmon = rdir / "star_salmon"
    salmon_dir = star_salmon if star_salmon.is_dir() else rdir / "salmon"

    if salmon_dir.is_dir():
        nf_counts = salmon_dir / "salmon.merged.gene_counts.tsv"
        nf_tpm = salmon_dir / "salmon.merged.gene_tpm.tsv"
        if nf_counts.is_file():
            counts_path = nf_counts
            source = "nf-core"
        if nf_tpm.is_file():
            tpm_path = nf_tpm

        for mqc_subdir in ("star_salmon", "salmon"):
            candidate = rdir / "multiqc" / mqc_subdir / "multiqc_report_data" / "multiqc_general_stats.txt"
            if candidate.is_file():
                mqc_path = candidate
                break

    # --- Fall back to scanning for loose tabular files ---
    if counts_path is None:
        candidates = []
        for p in sorted(rdir.iterdir()):
            if not p.is_file():
                continue
            stem_lower = p.stem.lower()
            suffixes = "".join(p.suffixes).lower()
            if suffixes not in _TABULAR_EXTENSIONS:
                continue
            if any(kw in stem_lower for kw in _COUNT_KEYWORDS):
                candidates.append(p)
        if len(candidates) == 1:
            counts_path = candidates[0]
            source = "user-provided"
        elif len(candidates) > 1:
            counts_path = candidates[0]
            source = "user-provided"

    # --- Look for design files ---
    for search_dir in (rdir, rdir.parent):
        for p in sorted(search_dir.iterdir()):
            if not p.is_file():
                continue
            stem_lower = p.stem.lower()
            suffixes = "".join(p.suffixes).lower()
            if suffixes in _TABULAR_EXTENSIONS and any(kw in stem_lower for kw in _DESIGN_KEYWORDS):
                design_path = p
                break
        if design_path:
            break
    if design_path is None:
        candidate = rdir.parent / "design.csv"
        if candidate.is_file():
            design_path = candidate

    # --- List all tabular files for the agent to inspect if needed ---
    tabular_files = []
    for p in sorted(rdir.iterdir()):
        if p.is_file() and "".join(p.suffixes).lower() in _TABULAR_EXTENSIONS:
            tabular_files.append(p.name)

    SESSION.results_dir = rdir

    return {
        "results_dir": str(rdir),
        "source": source,
        "counts_found": counts_path is not None,
        "counts_path": str(counts_path) if counts_path else None,
        "tpm_found": tpm_path is not None,
        "tpm_path": str(tpm_path) if tpm_path else None,
        "multiqc_found": mqc_path is not None,
        "multiqc_path": str(mqc_path) if mqc_path else None,
        "design_found": design_path is not None,
        "design_path": str(design_path) if design_path else None,
        "tabular_files": tabular_files,
    }


_FEATURECOUNTS_META = {"Geneid", "Chr", "Start", "End", "Strand", "Length"}


def _detect_separator(path: Path) -> str:
    """Guess CSV vs TSV from the first line."""
    first_line = path.open().readline()
    if "\t" in first_line:
        return "\t"
    return ","


def _load_design(design_path: str) -> tuple[pd.DataFrame | None, str | None]:
    """Load a design CSV/TSV. Returns (design_df, error_message)."""
    dp = Path(design_path)
    if not dp.is_file():
        return None, f"Design file '{design_path}' does not exist."
    sep = _detect_separator(dp)
    design = pd.read_csv(dp, sep=sep)
    if "sample" not in design.columns or "condition" not in design.columns:
        return None, "Design file must have 'sample' and 'condition' columns."
    design = design.set_index("sample")
    return design, None


def load_counts(counts_path: str, design_path: str | None = None) -> Summary:
    """Load a gene count matrix and optional design file.

    Handles multiple formats:
    - nf-core/rnaseq: TSV with gene_id, gene_name, then sample columns
    - featureCounts: TSV with Geneid, Chr, Start, End, Strand, Length, then samples
    - Generic CSV/TSV: first column as gene index, remaining as samples
    Auto-detects separator (tab vs comma).
    """
    path = Path(counts_path)
    if not path.is_file():
        return {"error": "not_a_file", "message": f"'{counts_path}' does not exist."}

    sep = _detect_separator(path)
    df = pd.read_csv(path, sep=sep)

    if df.shape[1] < 2:
        return {"error": "bad_format", "message": "Count matrix has fewer than 2 columns."}

    # Identify format and extract the count matrix
    gene_names = None
    gene_id_col = None

    if "gene_id" in df.columns:
        # nf-core format
        gene_id_col = "gene_id"
        if "gene_name" in df.columns:
            gene_names = df["gene_name"]
        meta_cols = {"gene_id", "gene_name"} & set(df.columns)

    elif "Geneid" in df.columns:
        # featureCounts format
        gene_id_col = "Geneid"
        meta_cols = _FEATURECOUNTS_META & set(df.columns)

    else:
        # Generic: first column is gene IDs, rest are samples
        gene_id_col = df.columns[0]
        meta_cols = {gene_id_col}

    gene_ids = df[gene_id_col]
    counts = df.drop(columns=list(meta_cols))

    # Check that remaining columns look numeric
    non_numeric = [c for c in counts.columns if not pd.api.types.is_numeric_dtype(counts[c])]
    if non_numeric:
        return {
            "error": "non_numeric_columns",
            "message": f"Columns {non_numeric[:5]} are not numeric. Check the file format — "
            "the count matrix should have gene IDs in one column and numeric values in the rest.",
        }

    counts = counts.set_index(gene_ids)

    SESSION.counts_df = counts
    SESSION.gene_names = gene_names.set_axis(gene_ids) if gene_names is not None else None

    # Load TPM if available (nf-core convention: same directory)
    nf_tpm = path.parent / "salmon.merged.gene_tpm.tsv"
    if nf_tpm.is_file():
        tpm_df = pd.read_csv(nf_tpm, sep="\t")
        tpm_meta = {"gene_id", "gene_name"} & set(tpm_df.columns)
        tpm_id_col = "gene_id" if "gene_id" in tpm_df.columns else tpm_df.columns[0]
        tpm = tpm_df.drop(columns=list(tpm_meta)).set_index(tpm_df[tpm_id_col])
        SESSION.tpm_df = tpm

    samples = list(counts.columns)
    lib_sizes = {s: round(float(counts[s].sum()), 2) for s in samples}

    # Load design if provided
    conditions = None
    if design_path:
        design, err = _load_design(design_path)
        if err:
            return {"error": "bad_design", "message": err}
        SESSION.design_df = design
        conditions = sorted(design["condition"].unique().tolist())

    detected_format = "nf-core" if "gene_id" in df.columns else (
        "featureCounts" if "Geneid" in df.columns else "generic"
    )

    return {
        "n_genes": int(counts.shape[0]),
        "n_samples": len(samples),
        "samples": samples,
        "library_sizes": lib_sizes,
        "has_gene_names": gene_names is not None,
        "tpm_loaded": SESSION.tpm_df is not None,
        "design_loaded": SESSION.design_df is not None,
        "conditions": conditions,
        "detected_format": detected_format,
    }


def inspect_counts() -> Summary:
    """Summarise the loaded count matrix so the agent can assess normalisation status.

    Returns per-sample and global statistics: value range, fraction of non-integer
    values, fraction of zeros, and distribution quantiles. These let the agent
    determine whether the data is raw counts, TPM/FPKM, or log-transformed.
    """
    if SESSION.counts_df is None:
        return {"error": "counts_not_loaded", "message": "Call load_counts before inspect_counts."}

    counts = SESSION.counts_df
    values = counts.values.ravel()
    n_values = values.size
    n_zero = int((values == 0).sum())
    n_non_integer = int(np.sum(values != np.floor(values)))

    global_stats = {
        "min": round(float(np.nanmin(values)), 4),
        "max": round(float(np.nanmax(values)), 4),
        "mean": round(float(np.nanmean(values)), 4),
        "median": round(float(np.nanmedian(values)), 4),
        "fraction_zero": round(n_zero / n_values, 4) if n_values else 0,
        "fraction_non_integer": round(n_non_integer / n_values, 4) if n_values else 0,
        "n_genes": int(counts.shape[0]),
        "n_samples": int(counts.shape[1]),
    }

    per_sample = {}
    for col in counts.columns:
        col_vals = counts[col].values
        per_sample[col] = {
            "mean": round(float(np.nanmean(col_vals)), 2),
            "median": round(float(np.nanmedian(col_vals)), 2),
            "max": round(float(np.nanmax(col_vals)), 2),
            "fraction_zero": round(float((col_vals == 0).sum() / len(col_vals)), 4),
        }

    hints = []
    if global_stats["fraction_non_integer"] < 0.01 and global_stats["max"] > 100:
        hints.append("Values are almost entirely integers with a wide range — likely raw counts.")
    elif global_stats["fraction_non_integer"] > 0.5 and global_stats["max"] < 25:
        hints.append("Mostly non-integer values in a narrow range — likely log-transformed. "
                      "DESeq2 requires raw counts; do not run it on this data.")
    elif global_stats["fraction_non_integer"] > 0.5 and global_stats["max"] > 100:
        hints.append("Non-integer values with a wide range — likely TPM/FPKM normalised. "
                      "DESeq2 requires raw counts; do not run it on this data.")

    return {
        "global": global_stats,
        "per_sample": per_sample,
        "hints": hints,
    }


def set_design(rows: list[dict]) -> Summary:
    """Define the sample-to-condition mapping when no design.csv exists.

    The agent infers conditions from sample names and the user's prompt, then calls
    this to set the mapping. Validates that sample names match the loaded count matrix.
    """
    if SESSION.counts_df is None:
        return {"error": "counts_not_loaded", "message": "Call load_counts before set_design."}
    if not rows:
        return {"error": "empty_design", "message": "Design must have at least one row."}

    for i, row in enumerate(rows):
        if "sample" not in row:
            return {"error": "missing_sample", "message": f"Row {i}: missing 'sample' key."}
        if "condition" not in row:
            return {"error": "missing_condition", "message": f"Row {i}: missing 'condition' key."}

    design = pd.DataFrame(rows).set_index("sample")
    count_samples = set(SESSION.counts_df.columns)
    design_samples = set(design.index)

    missing = count_samples - design_samples
    if missing:
        return {
            "error": "samples_missing_from_design",
            "message": f"Samples in counts but not in design: {sorted(missing)}",
        }

    SESSION.design_df = design

    # Write to disk for traceability
    paths = SESSION.require_paths()
    design.reset_index().to_csv(paths.design, index=False)

    conditions = sorted(design["condition"].unique().tolist())
    return {
        "n_samples": len(design),
        "conditions": conditions,
        "design_path": str(paths.design),
    }


def compute_qc() -> Summary:
    """Library sizes, gene detection rates, and PCA for outlier detection.

    Performs PCA on log2(counts + 1) so the agent can see whether samples cluster
    by condition or whether there are outliers. Parses MultiQC general stats if the
    results directory was scanned.
    """
    if SESSION.counts_df is None:
        return {"error": "counts_not_loaded", "message": "Call load_counts before compute_qc."}

    counts = SESSION.counts_df
    samples = list(counts.columns)

    lib_sizes = {s: int(counts[s].sum()) for s in samples}
    genes_detected = {s: int((counts[s] > 0).sum()) for s in samples}

    # PCA on log2(counts + 1)
    log_counts = np.log2(counts.values.astype(float).T + 1)  # samples x genes
    centered = log_counts - log_counts.mean(axis=0)
    U, S, _Vt = np.linalg.svd(centered, full_matrices=False)
    n_comps = min(2, len(S))
    pca_coords = (U[:, :n_comps] * S[:n_comps]).tolist()
    var_explained = ((S ** 2) / (S ** 2).sum())[:n_comps].tolist()
    pca_summary = {
        samples[i]: {"PC1": round(pca_coords[i][0], 4), "PC2": round(pca_coords[i][1], 4) if n_comps > 1 else 0.0}
        for i in range(len(samples))
    }

    # Parse MultiQC if available
    multiqc_summary = None
    if SESSION.results_dir:
        for mqc_subdir in ("star_salmon", "salmon"):
            mqc_path = SESSION.results_dir / "multiqc" / mqc_subdir / "multiqc_report_data" / "multiqc_general_stats.txt"
            if mqc_path.is_file():
                mqc = pd.read_csv(mqc_path, sep="\t")
                # Filter out per-read rows (e.g., "WT_REP1 Read 1")
                mqc = mqc[~mqc["Sample"].str.contains(" Read ", na=False)]
                key_cols = [c for c in mqc.columns if any(k in c for k in ("mapped_percent", "PERCENT_DUPLICATION", "total_sequences"))]
                if key_cols:
                    mqc_subset = mqc.set_index("Sample")[key_cols]
                    multiqc_summary = {
                        sample: {col.split("-")[-1]: round(float(val), 2) for col, val in row.items() if pd.notna(val)}
                        for sample, row in mqc_subset.iterrows()
                    }
                SESSION.multiqc_stats = mqc
                break

    return {
        "library_sizes": lib_sizes,
        "library_size_distribution": _quantile_summary(list(lib_sizes.values())),
        "genes_detected": genes_detected,
        "pca": pca_summary,
        "variance_explained": [round(v, 4) for v in var_explained],
        "multiqc_summary": multiqc_summary,
    }


def filter_low_counts(min_count: int = 10, min_samples: int = 2) -> Summary:
    """Drop genes where fewer than min_samples have at least min_count reads.

    Must run before DE so the multiple-testing correction isn't diluted by
    uninformative genes.
    """
    if SESSION.counts_df is None:
        return {"error": "counts_not_loaded", "message": "Call load_counts before filter_low_counts."}

    counts = SESSION.counts_df
    before = int(counts.shape[0])
    keep = (counts >= min_count).sum(axis=1) >= min_samples
    SESSION.counts_df = counts.loc[keep]
    after = int(SESSION.counts_df.shape[0])

    return {
        "min_count": min_count,
        "min_samples": min_samples,
        "genes_before": before,
        "genes_after": after,
        "genes_removed": before - after,
    }


def run_deseq2(contrast: list[str]) -> Summary:
    """Run PyDESeq2 for a given contrast.

    contrast is [factor, test, reference], e.g. ["condition", "treated", "control"].
    PyDESeq2 expects raw integer counts as samples (rows) x genes (cols).
    """
    if SESSION.counts_df is None:
        return {"error": "counts_not_loaded", "message": "Call load_counts before run_deseq2."}
    if SESSION.design_df is None:
        return {"error": "design_not_set", "message": "Call set_design or load_counts with a design before run_deseq2."}
    if len(contrast) != 3:
        return {"error": "bad_contrast", "message": "Contrast must be [factor, test, reference]."}

    factor, test, ref = contrast
    if factor not in SESSION.design_df.columns:
        return {"error": "bad_factor", "message": f"'{factor}' is not a column in the design."}

    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    # PyDESeq2 expects samples x genes
    counts_t = SESSION.counts_df.T
    # Align design to count matrix samples
    design = SESSION.design_df.loc[counts_t.index]

    dds = DeseqDataSet(counts=counts_t, metadata=design, design=f"~{factor}")
    dds.deseq2()

    stat_res = DeseqStats(dds, contrast=[factor, test, ref])
    stat_res.summary()
    results = stat_res.results_df.copy()

    SESSION.deseq_results = results

    # Write to disk
    paths = SESSION.require_paths()
    paths.analysis_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(paths.de_results)

    sig = results[results["padj"] < 0.05].dropna(subset=["padj"])
    n_up = int((sig["log2FoldChange"] > 0).sum())
    n_down = int((sig["log2FoldChange"] < 0).sum())

    # Top genes by padj
    top_up = sig[sig["log2FoldChange"] > 0].nsmallest(5, "padj")
    top_down = sig[sig["log2FoldChange"] < 0].nsmallest(5, "padj")

    def _gene_row(gene_id, row):
        out = {"gene_id": str(gene_id), "log2FC": round(float(row["log2FoldChange"]), 4), "padj": float(row["padj"])}
        if SESSION.gene_names is not None and gene_id in SESSION.gene_names.index:
            out["gene_name"] = str(SESSION.gene_names.loc[gene_id])
        return out

    return {
        "contrast": contrast,
        "n_tested": int(results["padj"].notna().sum()),
        "n_significant": int(len(sig)),
        "n_up": n_up,
        "n_down": n_down,
        "top_up": [_gene_row(gid, row) for gid, row in top_up.iterrows()],
        "top_down": [_gene_row(gid, row) for gid, row in top_down.iterrows()],
        "de_results_path": str(paths.de_results),
    }


def get_top_genes(n: int = 20, direction: str = "both") -> Summary:
    """Return top DE genes by adjusted p-value for the agent to inspect.

    direction: "up" (log2FC > 0), "down" (log2FC < 0), or "both".
    """
    if SESSION.deseq_results is None:
        return {"error": "no_deseq_results", "message": "Call run_deseq2 before get_top_genes."}

    results = SESSION.deseq_results.dropna(subset=["padj"])
    sig = results[results["padj"] < 0.05]

    if direction == "up":
        sig = sig[sig["log2FoldChange"] > 0]
    elif direction == "down":
        sig = sig[sig["log2FoldChange"] < 0]

    top = sig.nsmallest(n, "padj")
    genes = []
    for gene_id, row in top.iterrows():
        entry = {
            "gene_id": str(gene_id),
            "log2FC": round(float(row["log2FoldChange"]), 4),
            "padj": float(row["padj"]),
            "baseMean": round(float(row["baseMean"]), 2),
        }
        if SESSION.gene_names is not None and gene_id in SESSION.gene_names.index:
            entry["gene_name"] = str(SESSION.gene_names.loc[gene_id])
        genes.append(entry)

    return {
        "direction": direction,
        "n_returned": len(genes),
        "n_significant_total": int(len(sig)),
        "genes": genes,
    }


def run_enrichment(gene_list: list[str], organism: str = "human") -> Summary:
    """Run gene set enrichment via Enrichr (GO and KEGG).

    gene_list should contain gene symbols (not Ensembl IDs). Enrichr expects HGNC
    symbols for human, MGI symbols for mouse.
    """
    if not gene_list:
        return {"error": "empty_gene_list", "message": "gene_list is empty."}

    import gseapy

    gene_sets = ["GO_Biological_Process_2023", "KEGG_2021_Human"]
    if organism.lower() == "mouse":
        gene_sets = ["GO_Biological_Process_2023", "KEGG_2019_Mouse"]
    elif organism.lower() == "yeast":
        gene_sets = ["GO_Biological_Process_2023", "KEGG_2019"]

    try:
        enr = gseapy.enrichr(
            gene_list=gene_list,
            gene_sets=gene_sets,
            organism=organism,
            outdir=None,
            no_plot=True,
        )
    except Exception as exc:
        return {"error": "enrichr_failed", "message": str(exc)}

    results_df = enr.results
    SESSION.enrichment_results = {"organism": organism, "gene_sets": gene_sets}

    enrichment = {}
    for gs in gene_sets:
        subset = results_df[results_df["Gene_set"] == gs].nsmallest(10, "Adjusted P-value")
        enrichment[gs] = [
            {
                "term": str(row["Term"]),
                "pval": float(row["P-value"]),
                "padj": float(row["Adjusted P-value"]),
                "overlap": str(row["Overlap"]),
                "genes": str(row["Genes"]),
            }
            for _, row in subset.iterrows()
        ]

    SESSION.enrichment_results["results"] = enrichment
    return {
        "organism": organism,
        "gene_sets_queried": gene_sets,
        "n_input_genes": len(gene_list),
        "enrichment": enrichment,
    }


def summarize_findings() -> Summary:
    """Consolidate the final analysis state into one factual summary for the report."""
    out: Summary = {}

    if SESSION.counts_df is not None:
        out["n_genes"] = int(SESSION.counts_df.shape[0])
        out["n_samples"] = int(SESSION.counts_df.shape[1])

    if SESSION.design_df is not None:
        out["conditions"] = sorted(SESSION.design_df["condition"].unique().tolist())
        out["samples_per_condition"] = SESSION.design_df["condition"].value_counts().to_dict()

    if SESSION.deseq_results is not None:
        sig = SESSION.deseq_results.dropna(subset=["padj"])
        sig = sig[sig["padj"] < 0.05]
        out["de_summary"] = {
            "n_tested": int(SESSION.deseq_results["padj"].notna().sum()),
            "n_significant": int(len(sig)),
            "n_up": int((sig["log2FoldChange"] > 0).sum()),
            "n_down": int((sig["log2FoldChange"] < 0).sum()),
        }

    if SESSION.enrichment_results and "results" in SESSION.enrichment_results:
        top_terms = {}
        for gs, terms in SESSION.enrichment_results["results"].items():
            top_terms[gs] = [t["term"] for t in terms[:5]]
        out["top_enrichment_terms"] = top_terms

    return out


def generate_report(report_markdown: str) -> Summary:
    """Generate all analysis figures and write the report.

    The agent supplies the full report as Markdown with inline figure references
    (e.g. ![PCA](figures/pca.png)). This tool generates the figures at known paths
    and writes the report as-is.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable

    paths = SESSION.require_paths()
    paths.analysis_figures.mkdir(parents=True, exist_ok=True)
    figures: list[Path] = []

    # --- Consistent style + palette ---
    _PALETTE = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
                "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD"]
    _SIG_COLOR = "#C44E52"
    _NS_COLOR = "#B0B0B0"

    def _style_ax(ax, title=""):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=10)
        if title:
            ax.set_title(title, fontsize=13, fontweight="bold", pad=10)

    def _save(fig, name):
        p = paths.analysis_figures / name
        fig.savefig(p, dpi=150, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        figures.append(p)

    def _condition_colors(conditions_series):
        unique = list(dict.fromkeys(conditions_series))
        return {c: _PALETTE[i % len(_PALETTE)] for i, c in enumerate(unique)}

    # --- PCA (reusable computation) ---
    pca_data = None
    if SESSION.counts_df is not None:
        counts = SESSION.counts_df
        log_counts = np.log2(counts.values.astype(float).T + 1)
        centered = log_counts - log_counts.mean(axis=0)
        U, S, Vt = np.linalg.svd(centered, full_matrices=False)
        n_pcs = min(10, len(S))
        coords = U[:, :n_pcs] * S[:n_pcs]
        var_exp = (S ** 2) / (S ** 2).sum()
        samples = list(counts.columns)
        pca_data = {"coords": coords, "var_exp": var_exp, "samples": samples, "n_pcs": n_pcs}

    # --- Library size bar plot ---
    if SESSION.counts_df is not None:
        counts = SESSION.counts_df
        samples = list(counts.columns)
        lib_sizes = [float(counts[s].sum()) for s in samples]

        fig, ax = plt.subplots(figsize=(max(6, len(samples) * 0.8), 5))
        colors = [_NS_COLOR] * len(samples)
        if SESSION.design_df is not None:
            cmap = _condition_colors(SESSION.design_df.loc[samples, "condition"])
            colors = [cmap[SESSION.design_df.loc[s, "condition"]] for s in samples]
        bars = ax.bar(range(len(samples)), [s / 1e6 for s in lib_sizes], color=colors, edgecolor="white", linewidth=0.5)
        ax.set_xticks(range(len(samples)))
        ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=9)
        ax.set_ylabel("Library size (millions)", fontsize=11)
        _style_ax(ax, "Library sizes")
        if SESSION.design_df is not None:
            for cond, color in cmap.items():
                ax.bar([], [], color=color, label=cond)
            ax.legend(fontsize=9, frameon=False)
        _save(fig, "library_sizes.png")

    # --- PCA by condition ---
    if pca_data is not None:
        coords = pca_data["coords"]
        var_exp = pca_data["var_exp"]
        samples = pca_data["samples"]

        fig, ax = plt.subplots(figsize=(7, 6))
        if SESSION.design_df is not None:
            conditions = SESSION.design_df.loc[samples, "condition"]
            cmap = _condition_colors(conditions)
            for cond in dict.fromkeys(conditions):
                mask = [i for i, c in enumerate(conditions) if c == cond]
                ax.scatter(coords[mask, 0], coords[mask, 1] if pca_data["n_pcs"] > 1 else np.zeros(len(mask)),
                           c=cmap[cond], label=cond, s=70, edgecolors="white", linewidths=0.5, zorder=3)
            ax.legend(fontsize=9, frameon=False)
        else:
            ax.scatter(coords[:, 0], coords[:, 1] if pca_data["n_pcs"] > 1 else np.zeros(len(samples)),
                       c=_PALETTE[0], s=70, edgecolors="white", linewidths=0.5, zorder=3)
        for i, s in enumerate(samples):
            ax.annotate(s, (coords[i, 0], coords[i, 1] if pca_data["n_pcs"] > 1 else 0),
                        fontsize=8, alpha=0.7, textcoords="offset points", xytext=(5, 5))
        ax.set_xlabel(f"PC1 ({var_exp[0]:.1%} variance)", fontsize=11)
        ax.set_ylabel(f"PC2 ({var_exp[1]:.1%} variance)" if pca_data["n_pcs"] > 1 else "PC2", fontsize=11)
        _style_ax(ax, "PCA — log2(counts + 1)")
        _save(fig, "pca.png")

    # --- PCA coloured by other metadata variables ---
    if pca_data is not None and SESSION.design_df is not None:
        extra_cols = [c for c in SESSION.design_df.columns if c != "condition"]
        for col in extra_cols:
            vals = SESSION.design_df.loc[samples, col]
            fig, ax = plt.subplots(figsize=(7, 6))

            if vals.dtype == object or vals.nunique() <= 8:
                unique_vals = list(dict.fromkeys(vals))
                col_cmap = {v: _PALETTE[i % len(_PALETTE)] for i, v in enumerate(unique_vals)}
                for val in unique_vals:
                    mask = [i for i, v in enumerate(vals) if v == val]
                    ax.scatter(coords[mask, 0], coords[mask, 1], c=col_cmap[val],
                               label=str(val), s=70, edgecolors="white", linewidths=0.5, zorder=3)
                ax.legend(fontsize=9, frameon=False, title=col)
            else:
                numeric_vals = pd.to_numeric(vals, errors="coerce")
                sc = ax.scatter(coords[:, 0], coords[:, 1], c=numeric_vals, cmap="viridis",
                                s=70, edgecolors="white", linewidths=0.5, zorder=3)
                plt.colorbar(sc, ax=ax, label=col)

            for i, s in enumerate(samples):
                ax.annotate(s, (coords[i, 0], coords[i, 1]),
                            fontsize=8, alpha=0.7, textcoords="offset points", xytext=(5, 5))
            ax.set_xlabel(f"PC1 ({var_exp[0]:.1%} variance)", fontsize=11)
            ax.set_ylabel(f"PC2 ({var_exp[1]:.1%} variance)", fontsize=11)
            _style_ax(ax, f"PCA coloured by {col}")
            _save(fig, f"pca_{col.lower().replace(' ', '_')}.png")

    # --- PC–metadata association heatmap ---
    if pca_data is not None and SESSION.design_df is not None:
        design_cols = list(SESSION.design_df.columns)
        n_pcs_show = min(pca_data["n_pcs"], 10)
        if design_cols and n_pcs_show >= 2:
            from scipy import stats as sp_stats
            assoc_matrix = np.full((len(design_cols), n_pcs_show), np.nan)
            for ci, col in enumerate(design_cols):
                vals = SESSION.design_df.loc[samples, col]
                for pc_i in range(n_pcs_show):
                    pc_vals = coords[:, pc_i]
                    numeric_vals = pd.to_numeric(vals, errors="coerce")
                    if numeric_vals.notna().all():
                        r, _ = sp_stats.pearsonr(numeric_vals.values, pc_vals)
                        assoc_matrix[ci, pc_i] = r ** 2
                    elif vals.nunique() > 1:
                        groups = [pc_vals[[j for j in range(len(vals)) if vals.iloc[j] == g]]
                                  for g in vals.unique() if sum(vals == g) > 0]
                        if len(groups) >= 2 and all(len(g) > 0 for g in groups):
                            f_stat, p_val = sp_stats.f_oneway(*groups)
                            ss_between = sum(len(g) * (g.mean() - pc_vals.mean()) ** 2 for g in groups)
                            ss_total = np.sum((pc_vals - pc_vals.mean()) ** 2)
                            assoc_matrix[ci, pc_i] = ss_between / ss_total if ss_total > 0 else 0

            fig, ax = plt.subplots(figsize=(max(6, n_pcs_show * 0.8), max(3, len(design_cols) * 0.6 + 1)))
            im = ax.imshow(assoc_matrix, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)
            ax.set_xticks(range(n_pcs_show))
            ax.set_xticklabels([f"PC{i+1}\n({pca_data['var_exp'][i]:.1%})" for i in range(n_pcs_show)], fontsize=9)
            ax.set_yticks(range(len(design_cols)))
            ax.set_yticklabels(design_cols, fontsize=10)
            for ci in range(len(design_cols)):
                for pi in range(n_pcs_show):
                    val = assoc_matrix[ci, pi]
                    if not np.isnan(val):
                        ax.text(pi, ci, f"{val:.2f}", ha="center", va="center",
                                fontsize=8, color="white" if val > 0.5 else "black")
            plt.colorbar(im, ax=ax, label="R² (association)", shrink=0.8)
            _style_ax(ax, "PC–metadata associations")
            ax.spines["left"].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            _save(fig, "pc_association.png")

    # --- Sample correlation heatmap ---
    if SESSION.counts_df is not None:
        counts = SESSION.counts_df
        log_counts_df = np.log2(counts.astype(float) + 1)
        corr = log_counts_df.corr(method="pearson")
        samples = list(corr.columns)
        n = len(samples)

        fig, ax = plt.subplots(figsize=(max(5, n * 0.6), max(4, n * 0.5)))
        im = ax.imshow(corr.values, cmap="RdYlBu_r", vmin=corr.values.min(), vmax=1)
        ax.set_xticks(range(n))
        ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=9)
        ax.set_yticks(range(n))
        ax.set_yticklabels(samples, fontsize=9)
        for i in range(n):
            for j in range(n):
                val = corr.values[i, j]
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7,
                        color="white" if val < 0.95 else "black")
        plt.colorbar(im, ax=ax, label="Pearson r", shrink=0.8)
        _style_ax(ax, "Sample correlation — log2(counts + 1)")
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)
        _save(fig, "sample_correlation.png")

    # --- Volcano plot ---
    if SESSION.deseq_results is not None:
        results = SESSION.deseq_results.dropna(subset=["padj"])
        sig_mask = results["padj"] < 0.05

        fig, ax = plt.subplots(figsize=(7, 6))
        ax.scatter(results.loc[~sig_mask, "log2FoldChange"],
                   -np.log10(results.loc[~sig_mask, "padj"]),
                   c=_NS_COLOR, alpha=0.4, s=12, label="NS", zorder=2)
        ax.scatter(results.loc[sig_mask, "log2FoldChange"],
                   -np.log10(results.loc[sig_mask, "padj"]),
                   c=_SIG_COLOR, alpha=0.6, s=12, label="padj < 0.05", zorder=3)
        ax.axhline(-np.log10(0.05), color="#888888", linewidth=0.8, linestyle="--", alpha=0.5)
        ax.set_xlabel("log2 Fold Change", fontsize=11)
        ax.set_ylabel("-log10(padj)", fontsize=11)
        ax.legend(fontsize=9, frameon=False)
        _style_ax(ax, "Volcano plot")
        _save(fig, "volcano.png")

        # --- MA plot ---
        fig, ax = plt.subplots(figsize=(7, 6))
        ax.scatter(np.log10(results.loc[~sig_mask, "baseMean"] + 1),
                   results.loc[~sig_mask, "log2FoldChange"],
                   c=_NS_COLOR, alpha=0.4, s=12, label="NS", zorder=2)
        ax.scatter(np.log10(results.loc[sig_mask, "baseMean"] + 1),
                   results.loc[sig_mask, "log2FoldChange"],
                   c=_SIG_COLOR, alpha=0.6, s=12, label="padj < 0.05", zorder=3)
        ax.axhline(0, color="#333333", linewidth=0.8)
        ax.set_xlabel("log10(baseMean + 1)", fontsize=11)
        ax.set_ylabel("log2 Fold Change", fontsize=11)
        ax.legend(fontsize=9, frameon=False)
        _style_ax(ax, "MA plot")
        _save(fig, "ma_plot.png")

    # --- Top DE genes heatmap ---
    if SESSION.deseq_results is not None and SESSION.counts_df is not None:
        sig = SESSION.deseq_results.dropna(subset=["padj"])
        sig = sig[sig["padj"] < 0.05]
        if len(sig) > 0:
            n_top = min(40, len(sig))
            top_genes = sig.nsmallest(n_top, "padj").index
            counts_top = SESSION.counts_df.loc[top_genes]
            log_vals = np.log2(counts_top.values.astype(float) + 1)
            row_means = log_vals.mean(axis=1, keepdims=True)
            row_stds = log_vals.std(axis=1, keepdims=True)
            row_stds[row_stds == 0] = 1
            z_scores = (log_vals - row_means) / row_stds

            gene_labels = []
            for gid in top_genes:
                if SESSION.gene_names is not None and gid in SESSION.gene_names.index:
                    gene_labels.append(str(SESSION.gene_names.loc[gid]))
                else:
                    gene_labels.append(str(gid))

            fig, ax = plt.subplots(figsize=(max(5, len(samples) * 0.7), max(6, n_top * 0.25)))
            im = ax.imshow(z_scores, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
            ax.set_xticks(range(len(samples)))
            ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=9)
            ax.set_yticks(range(n_top))
            ax.set_yticklabels(gene_labels, fontsize=7)
            plt.colorbar(im, ax=ax, label="z-score", shrink=0.6)
            _style_ax(ax, f"Top {n_top} DE genes (z-scored)")
            ax.spines["left"].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            _save(fig, "de_heatmap.png")

    # --- Enrichment bar plots ---
    if SESSION.enrichment_results and "results" in SESSION.enrichment_results:
        for gs, terms in SESSION.enrichment_results["results"].items():
            if not terms:
                continue
            top = terms[:10]
            fig, ax = plt.subplots(figsize=(9, max(3, len(top) * 0.4)))
            term_names = [t["term"][:65] for t in reversed(top)]
            pvals = [-np.log10(t["padj"]) if t["padj"] > 0 else 10 for t in reversed(top)]
            ax.barh(term_names, pvals, color=_PALETTE[0], edgecolor="white", linewidth=0.5)
            ax.set_xlabel("-log10(padj)", fontsize=11)
            _style_ax(ax, gs.replace("_", " "))
            fig_name = f"enrichment_{gs.lower().replace(' ', '_')}.png"
            _save(fig, fig_name)

    # --- Write report ---
    paths.analysis_report.write_text(report_markdown.rstrip() + "\n")

    return {
        "report_path": str(paths.analysis_report),
        "figures": {p.stem: f"figures/{p.name}" for p in figures},
        "n_figures": len(figures),
    }

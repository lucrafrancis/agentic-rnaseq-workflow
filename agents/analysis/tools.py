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


def scan_results(results_dir: str) -> Summary:
    """Discover nf-core/rnaseq outputs: count matrices, TPM, MultiQC data.

    Non-mutating; the agent calls this first to orient itself. Also checks the parent
    run directory for a design.csv from the samplesheet agent.
    """
    rdir = Path(results_dir)
    if not rdir.is_dir():
        return {"error": "not_a_directory", "message": f"'{results_dir}' is not a directory."}

    # nf-core uses star_salmon/ by default; also check salmon/ for salmon-only runs
    star_salmon = rdir / "star_salmon"
    salmon_dir = star_salmon if star_salmon.is_dir() else rdir / "salmon"

    counts_path = salmon_dir / "salmon.merged.gene_counts.tsv" if salmon_dir.is_dir() else None
    tpm_path = salmon_dir / "salmon.merged.gene_tpm.tsv" if salmon_dir.is_dir() else None

    counts_found = counts_path is not None and counts_path.is_file()
    tpm_found = tpm_path is not None and tpm_path.is_file()

    # MultiQC: look for the general stats file
    mqc_path = None
    for mqc_subdir in ("star_salmon", "salmon"):
        candidate = rdir / "multiqc" / mqc_subdir / "multiqc_report_data" / "multiqc_general_stats.txt"
        if candidate.is_file():
            mqc_path = candidate
            break

    # Check parent directory for design.csv (from samplesheet agent)
    design_path = rdir.parent / "design.csv"
    design_found = design_path.is_file()

    SESSION.results_dir = rdir

    return {
        "results_dir": str(rdir),
        "counts_found": counts_found,
        "counts_path": str(counts_path) if counts_found else None,
        "tpm_found": tpm_found,
        "tpm_path": str(tpm_path) if tpm_found else None,
        "multiqc_found": mqc_path is not None,
        "multiqc_path": str(mqc_path) if mqc_path else None,
        "design_found": design_found,
        "design_path": str(design_path) if design_found else None,
    }


def load_counts(counts_path: str, design_path: str | None = None) -> Summary:
    """Load the gene count matrix and optional design CSV.

    The nf-core count matrix has columns: gene_id, gene_name, then one column per
    sample with integer counts. Sets gene_id as index, stores gene_name separately.
    """
    path = Path(counts_path)
    if not path.is_file():
        return {"error": "not_a_file", "message": f"'{counts_path}' does not exist."}

    df = pd.read_csv(path, sep="\t")
    if "gene_id" not in df.columns:
        return {"error": "bad_format", "message": "Count matrix missing 'gene_id' column."}

    gene_names = df["gene_name"] if "gene_name" in df.columns else None
    meta_cols = {"gene_id", "gene_name"} & set(df.columns)
    counts = df.drop(columns=list(meta_cols)).set_index(df["gene_id"])
    counts = counts.astype(int)

    SESSION.counts_df = counts
    SESSION.gene_names = gene_names.set_axis(df["gene_id"]) if gene_names is not None else None

    # Load TPM if available (same directory, known filename)
    tpm_path = path.parent / "salmon.merged.gene_tpm.tsv"
    if tpm_path.is_file():
        tpm_df = pd.read_csv(tpm_path, sep="\t")
        tpm = tpm_df.drop(columns=list(meta_cols & set(tpm_df.columns))).set_index(tpm_df["gene_id"])
        SESSION.tpm_df = tpm

    samples = list(counts.columns)
    lib_sizes = {s: int(counts[s].sum()) for s in samples}

    # Load design if provided
    conditions = None
    if design_path:
        dp = Path(design_path)
        if not dp.is_file():
            return {"error": "design_not_found", "message": f"Design file '{design_path}' does not exist."}
        design = pd.read_csv(dp)
        if "sample" not in design.columns or "condition" not in design.columns:
            return {"error": "bad_design", "message": "Design CSV must have 'sample' and 'condition' columns."}
        design = design.set_index("sample")
        SESSION.design_df = design
        conditions = sorted(design["condition"].unique().tolist())

    return {
        "n_genes": int(counts.shape[0]),
        "n_samples": len(samples),
        "samples": samples,
        "library_sizes": lib_sizes,
        "has_gene_names": gene_names is not None,
        "tpm_loaded": SESSION.tpm_df is not None,
        "design_loaded": SESSION.design_df is not None,
        "conditions": conditions,
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
    """Render figures, write DE results, and assemble the analysis report.

    The agent supplies the narrative (report_markdown); this tool handles the
    deterministic artifacts: figures and stitching them into the report.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    paths = SESSION.require_paths()
    paths.analysis_figures.mkdir(parents=True, exist_ok=True)
    figures: list[Path] = []

    # Volcano plot
    if SESSION.deseq_results is not None:
        results = SESSION.deseq_results.dropna(subset=["padj"])
        fig, ax = plt.subplots(figsize=(8, 6))
        sig_mask = results["padj"] < 0.05
        ax.scatter(
            results.loc[~sig_mask, "log2FoldChange"],
            -np.log10(results.loc[~sig_mask, "padj"]),
            c="grey", alpha=0.5, s=10, label="NS",
        )
        ax.scatter(
            results.loc[sig_mask, "log2FoldChange"],
            -np.log10(results.loc[sig_mask, "padj"]),
            c="red", alpha=0.6, s=10, label="padj < 0.05",
        )
        ax.set_xlabel("log2 Fold Change")
        ax.set_ylabel("-log10(padj)")
        ax.legend()
        ax.set_title("Volcano Plot")
        path = paths.analysis_figures / "volcano.png"
        fig.savefig(path, dpi=120, bbox_inches="tight")
        plt.close(fig)
        figures.append(path)

        # MA plot
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.scatter(
            np.log10(results.loc[~sig_mask, "baseMean"] + 1),
            results.loc[~sig_mask, "log2FoldChange"],
            c="grey", alpha=0.5, s=10, label="NS",
        )
        ax.scatter(
            np.log10(results.loc[sig_mask, "baseMean"] + 1),
            results.loc[sig_mask, "log2FoldChange"],
            c="red", alpha=0.6, s=10, label="padj < 0.05",
        )
        ax.set_xlabel("log10(baseMean + 1)")
        ax.set_ylabel("log2 Fold Change")
        ax.axhline(0, color="black", linewidth=0.5)
        ax.legend()
        ax.set_title("MA Plot")
        path = paths.analysis_figures / "ma_plot.png"
        fig.savefig(path, dpi=120, bbox_inches="tight")
        plt.close(fig)
        figures.append(path)

    # PCA plot
    if SESSION.counts_df is not None:
        counts = SESSION.counts_df
        log_counts = np.log2(counts.values.astype(float).T + 1)
        centered = log_counts - log_counts.mean(axis=0)
        U, S, _Vt = np.linalg.svd(centered, full_matrices=False)
        n_comps = min(2, len(S))
        coords = U[:, :n_comps] * S[:n_comps]
        var_exp = (S ** 2) / (S ** 2).sum()

        fig, ax = plt.subplots(figsize=(8, 6))
        samples = list(counts.columns)

        if SESSION.design_df is not None:
            conditions = SESSION.design_df.loc[samples, "condition"]
            for cond in conditions.unique():
                mask = conditions == cond
                idx = [i for i, m in enumerate(mask) if m]
                ax.scatter(
                    coords[idx, 0], coords[idx, 1] if n_comps > 1 else np.zeros(len(idx)),
                    label=cond, s=60,
                )
        else:
            ax.scatter(coords[:, 0], coords[:, 1] if n_comps > 1 else np.zeros(len(samples)), s=60)

        for i, s in enumerate(samples):
            ax.annotate(s, (coords[i, 0], coords[i, 1] if n_comps > 1 else 0), fontsize=7, alpha=0.7)

        ax.set_xlabel(f"PC1 ({var_exp[0]:.1%} variance)")
        ax.set_ylabel(f"PC2 ({var_exp[1]:.1%} variance)" if n_comps > 1 else "PC2")
        ax.legend()
        ax.set_title("PCA — log2(counts + 1)")
        path = paths.analysis_figures / "pca.png"
        fig.savefig(path, dpi=120, bbox_inches="tight")
        plt.close(fig)
        figures.append(path)

    # Enrichment bar plot
    if SESSION.enrichment_results and "results" in SESSION.enrichment_results:
        for gs, terms in SESSION.enrichment_results["results"].items():
            if not terms:
                continue
            top = terms[:10]
            fig, ax = plt.subplots(figsize=(10, 6))
            term_names = [t["term"][:60] for t in reversed(top)]
            pvals = [-np.log10(t["padj"]) if t["padj"] > 0 else 10 for t in reversed(top)]
            ax.barh(term_names, pvals)
            ax.set_xlabel("-log10(padj)")
            ax.set_title(gs.replace("_", " "))
            path = paths.analysis_figures / f"enrichment_{gs.lower().replace(' ', '_')}.png"
            fig.savefig(path, dpi=120, bbox_inches="tight")
            plt.close(fig)
            figures.append(path)

    # Assemble report
    figures_md = ""
    if figures:
        figures_md = "\n\n## Figures\n\n" + "\n\n".join(
            f"![{p.stem}](figures/{p.name})" for p in figures
        )
    paths.analysis_report.write_text(report_markdown.rstrip() + figures_md + "\n")

    return {
        "report_path": str(paths.analysis_report),
        "figures": [str(p) for p in figures],
        "report_chars": len(report_markdown),
    }

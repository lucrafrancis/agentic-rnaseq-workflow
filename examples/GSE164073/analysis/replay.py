"""Replay of the analysis agent's tool calls for run 20260930_GSE164073_counts_4.

Generated: 2026-09-30T10:58:35
Git commit: dbce4ddf4fd2dad866a3009cd89a5a0a8791fc3d-dirty
Source log: tool_calls.jsonl (10 of 19 calls; failed and read-only calls omitted)
Copied to examples/GSE164073/ for the repository: data paths were changed from runs/<run>/ to
examples/GSE164073/ (the only edit; the tool logs are unchanged).

Re-runs the same tool functions with the same arguments, without the LLM. Output goes
to a new run directory (runs/<date>_20260930_GSE164073_counts_4_replay/). Package versions are pinned by
uv.lock at the commit above — check it out first if the code has changed since.

Caveat: run_enrichment queries Enrichr live and fetch_abstract queries PubMed; Enrichr's
gene-set libraries change over time, so enrichment results may differ from the original
run (the original inputs and results are saved in analysis/enrichment_*.{txt,csv}).

  uv run python examples/GSE164073/analysis/replay.py
"""

import os
import sys
from pathlib import Path

# The repo this script sits in (runs/<run>/analysis/, or a copy under examples/)
ROOT = next(p for p in Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)  # paths below are relative to the repo root

from agents.analysis import tools
from core.session import SESSION

NETWORK_TOOLS = ['fetch_abstract', 'run_enrichment']


def _run(name, **kwargs):
    print(f"-> {name}")
    result = getattr(tools, name)(**kwargs)
    if "error" in result:
        if name in NETWORK_TOOLS:
            print(f"   WARNING: {result['error']}: {result.get('message', '')}")
            return result
        sys.exit(f"   FAILED: {result['error']}: {result.get('message', '')}")
    return result


SESSION.begin_run('20260930_GSE164073_counts_4_replay')
print(f"Replay output: {SESSION.paths.dir}")

_run(
    "scan_results",
    results_dir='examples/GSE164073',
)

_run(
    "fetch_abstract",
    pmid='34022129',
)

_run(
    "load_counts",
    counts_path='examples/GSE164073/counts.tsv',
    design_path='examples/GSE164073/design.csv',
)

_run("compute_qc")

_run("filter_low_counts")

_run(
    "run_deseq2",
    contrast=['condition', 'SARS-CoV-2, MOI = 1.0', 'mock'],
    covariates=['tissue'],
)

_run(
    "run_enrichment",
    direction='up',
)

_run(
    "run_enrichment",
    direction='down',
)

_run("generate_figures")

_run(
    "write_report",
    report_markdown=('## 1. Executive Summary\n'
     '\n'
     'Raw gene counts from {{fact.provenance.accession}} ({{fact.provenance.source}}, file '
     '{{fact.provenance.file}}) were analysed to compare SARS-CoV-2-infected human ocular '
     'surface tissue against mock-infected controls. Using a design that adjusts for '
     'tissue of origin ({{fact.de.design}}), differential expression identified '
     '{{de.n_significant}} genes changed (padj < {{fact.de.padj_threshold}}: {{de.n_up}} '
     'up, {{de.n_down}} down out of {{fact.de.n_tested}} tested). Upregulated genes are '
     'strongly enriched for inflammatory/NF-κB and chemokine signalling '
     '({{term:up:Inflammatory Response (GO:0006954)}}, {{term:up:NF-kappa B signaling '
     "pathway}}), consistent with the associated study's report of an NF-κB-driven "
     'response to infection ({{cite:34022129}}), while downregulated genes are enriched '
     'for epithelial differentiation programs ({{term:down:Epithelium Development '
     '(GO:0060429)}}).\n'
     '\n'
     '## 2. Experimental Design\n'
     '\n'
     'This dataset profiles the transcriptional response of human ocular surface tissue '
     '(cornea, limbus, sclera; donor cadaver explants) to SARS-CoV-2 infection versus mock '
     'treatment, as described in {{cite:34022129}}. The loaded design covers {{samples.n}} '
     'samples across the conditions {{fact.samples.conditions}} ({{fact.samples.n_mock}} '
     'mock, {{fact.samples.n_sars_cov_2_moi_1_0}} infected), with tissue levels '
     '{{design.tissue_values}} and a time point of {{design.time_point}}. Because tissue '
     'is a major source of variation shared across both conditions (three tissues × two '
     'conditions × triplicate donors), tissue was included as a covariate in the DESeq2 '
     'model (design {{fact.de.design}}) so that the SARS-CoV-2-vs-mock contrast is '
     'estimated within each tissue and pooled across tissues.\n'
     '\n'
     '## 3. Quality Control\n'
     '\n'
     '{{table:qc}}\n'
     '\n'
     'Library sizes ranged from {{fact.qc.library_size_min_millions}} to '
     '{{fact.qc.library_size_max_millions}} million reads (median '
     '{{fact.qc.library_size_median_millions}} million), and the minimum number of genes '
     'detected in any sample was {{fact.qc.genes_detected_min}} (maximum '
     '{{fact.qc.genes_detected_max}}). Pairwise sample correlations ranged from '
     '{{fact.qc.sample_correlation_min}} to {{fact.qc.sample_correlation_max}}.\n'
     '\n'
     '![Library sizes per sample, coloured by condition](figures/library_sizes.png)\n'
     '\n'
     '![PCA of samples coloured by condition](figures/pca.png)\n'
     '\n'
     '![PCA of samples coloured by tissue](figures/pca_tissue.png)\n'
     '\n'
     '![PC–design variable associations](figures/pc_association.png)\n'
     '\n'
     '![Sample-to-sample correlation heatmap](figures/sample_correlation.png)\n'
     '\n'
     '{{fact.qc.pca_summary}} As the PCA plots show, PC1 and PC2 are dominated by tissue '
     'identity ({{fact.qc.pca_r2_pc1.tissue}} and {{fact.qc.pca_r2_pc2.tissue}} of '
     'variance respectively) rather than infection status, which is why tissue was '
     'modelled as a covariate rather than left unadjusted. Within-tissue replicate spread '
     'is similar between the two conditions ({{fact.qc.pca_spread.sars_cov_2_moi_1_0}} vs '
     '{{fact.qc.pca_spread.mock}}), and only two samples ({{fact.qc.pca_misclustered}}) '
     "fall closer to the other condition's centroid within their tissue group. No samples "
     'were excluded as outliers.\n'
     '\n'
     '## 4. Differential Expression\n'
     '\n'
     '{{table:de_summary}}\n'
     '\n'
     '![Volcano plot of differential expression](figures/volcano.png)\n'
     '\n'
     '![MA plot of differential expression](figures/ma_plot.png)\n'
     '\n'
     'Top upregulated genes include {{gene:SOD2}}, {{gene:TNFAIP3}} and {{gene:C3}}; top '
     'downregulated genes include {{gene:ACAN}} and {{gene:ACTC1}}.\n'
     '\n'
     '{{table:top_up}}\n'
     '\n'
     '{{table:top_down}}\n'
     '\n'
     '![Heatmap of top differentially expressed genes](figures/de_heatmap.png)\n'
     '\n'
     '## 5. Gene Set Enrichment\n'
     '\n'
     '{{table:enrichment_up}}\n'
     '\n'
     '![GO Biological Process enrichment, upregulated '
     'genes](figures/enrichment_upregulated_go_biological_process_2023.png)\n'
     '\n'
     '![KEGG enrichment, upregulated '
     'genes](figures/enrichment_upregulated_kegg_2021_human.png)\n'
     '\n'
     '{{table:enrichment_down}}\n'
     '\n'
     '![GO Biological Process enrichment, downregulated '
     'genes](figures/enrichment_downregulated_go_biological_process_2023.png)\n'
     '\n'
     '![KEGG enrichment, downregulated '
     'genes](figures/enrichment_downregulated_kegg_2021_human.png)\n'
     '\n'
     'Upregulated genes are dominated by inflammatory and cytokine/chemokine terms such as '
     '{{term:up:Inflammatory Response (GO:0006954)}}, {{term:up:Chemokine-Mediated '
     'Signaling Pathway (GO:0070098)}} and {{term:up:NF-kappa B signaling pathway}}, '
     'alongside {{term:up:Coronavirus disease}}. Downregulated genes are enriched for '
     '{{term:down:Epithelium Development (GO:0060429)}} and related '
     'epithelial/keratinocyte differentiation terms, indicating a loss of epithelial '
     'identity programs alongside the inflammatory response.\n'
     '\n'
     '## 6. Comparison with the Published Study\n'
     '\n'
     'The associated study ({{cite:34022129}}) reports that SARS-CoV-2-infected ocular '
     'surface tissue shows "robust induction of NF-κB in infected cells as well as '
     'diminished type I/III interferon signaling." This analysis agrees with the '
     'NF-κB/inflammatory induction: {{gene:RELB}}, {{gene:NFKBIA}} and {{gene:NFKBIZ}} are '
     'all significantly upregulated, and the enrichment results independently recover '
     '{{term:up:NF-kappa B signaling pathway}} among upregulated genes. Regarding '
     'interferon signalling, the picture in this contrast is more mixed than "diminished": '
     'canonical ISGs such as {{gene:ISG15}}, {{gene:IFIT1}}, {{gene:IFIT3}}, '
     '{{gene:OAS1}}, {{gene:OAS2}} and {{gene:OAS3}} are not significant in this analysis, '
     'while a few others ({{gene:MX1}}, {{gene:IRF7}}, {{gene:STAT1}}) are significantly, '
     'though modestly, upregulated. Type I/III interferon genes themselves '
     '({{fact.samples.conditions}} aside) were not testable — IFNB1, IFNL1, IFNL2 and '
     'IFNL3 were filtered out before DE due to low counts, so this dataset cannot directly '
     'confirm or refute the reported interferon suppression at the ligand level.\n'
     '\n'
     '## 7. Biological Interpretation\n'
     '\n'
     '- The dominant transcriptional signature of SARS-CoV-2 infection across ocular '
     'tissues is an innate-immune/inflammatory chemokine response: {{term:up:Cellular '
     'Response To Chemokine (GO:1990869)}} and {{term:up:Cytokine-cytokine receptor '
     'interaction}} are both strongly enriched among upregulated genes, consistent with '
     'recruitment of innate immune cells to infected tissue.\n'
     '- NF-κB pathway activation is evident at the gene level ({{gene:RELB}}, '
     '{{gene:NFKBIA}}, {{gene:TNFAIP3}} all significantly upregulated) and is '
     'independently supported by {{term:up:NF-kappa B signaling pathway}} enrichment, '
     'matching the qualitative direction reported in {{cite:34022129}}.\n'
     '- Classical interferon-stimulated genes ({{gene:ISG15}}, {{gene:IFIT1}}, '
     '{{gene:OAS1}}) are not significantly changed in this contrast, which does not '
     'straightforwardly support a strong interferon response either way in the '
     'pooled-tissue analysis; this should be interpreted cautiously given the '
     'covariate-adjusted, tissue-pooled design rather than tissue-specific comparisons.\n'
     '- Downregulated genes point to a loss of epithelial/keratinocyte differentiation '
     'programs ({{term:down:Epithelium Development (GO:0060429)}}), a pattern not '
     'discussed in the abstract but visible directly in this dataset.\n'
     '- Because tissue explains the large majority of overall expression variance '
     '({{fact.qc.pca_r2_pc1.tissue}} of PC1), condition-driven effects were only '
     'detectable after adjusting for tissue in the model; the infection signature '
     'described above is the shared, tissue-adjusted component of the response.\n'
     '\n'
     '## 8. Limitations\n'
     '\n'
     '- Each tissue × condition group has a modest number of donor replicates, and PCA '
     'shows replicate spread within condition groups is not negligible '
     '({{fact.qc.pca_spread.sars_cov_2_moi_1_0}} vs {{fact.qc.pca_spread.mock}}), with a '
     'small number of samples ({{fact.qc.pca_misclustered}}) landing closer to the '
     "opposite condition's centroid within their tissue.\n"
     '- The DESeq2 model pools cornea, limbus and sclera and adjusts for tissue as a '
     'covariate; it estimates a shared infection effect across tissues rather than '
     'tissue-specific responses, even though the source study highlights the limbus '
     'specifically as a site of productive infection.\n'
     '- Type I/III interferon ligand genes (IFNB1, IFNL1/2/3) were filtered out for low '
     'counts before DE and could not be evaluated, limiting direct comparison with the '
     "paper's interferon-suppression claim.\n"
     "- Enrichment analysis used Enrichr's whole-genome background and the default "
     'significance/fold-change thresholds ({{fact.enrichment.up.padj_max}} padj, '
     '{{fact.enrichment.up.lfc_min}} log2FC) rather than a tissue-matched expressed-gene '
     'background.\n'
     '- Data originate from a processed GEO count matrix rather than raw FASTQ '
     'reprocessing, so upstream alignment/quantification choices made by the original '
     'authors could not be independently verified.\n'
     '\n'
     '## 9. Methods\n'
     '\n'
     '{{table:methods}}\n'
     '\n'
     '## 10. References\n'
     '\n'
     '- {{cite:34022129}}\n'),
)

print(f"Done. Report: {SESSION.paths.analysis_report}")

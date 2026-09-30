"""Replay of the analysis agent's tool calls for run 20260930_GSE246386_full.

Generated: 2026-09-30T13:59:11
Git commit: b2ebf7493e98afd85c8e0709719dd15aec4e29f2
Source log: tool_calls.jsonl (11 of 23 calls; failed and read-only calls omitted)
Copied to examples/GSE246386/ for the repository: data paths were changed from runs/<run>/ to
examples/GSE246386/ (the only edit to this script). results/ is trimmed to the files the analysis reads.

Re-runs the same tool functions with the same arguments, without the LLM. Output goes
to a new run directory (runs/<date>_20260930_GSE246386_full_replay/). Package versions are pinned by
uv.lock at the commit above — check it out first if the code has changed since.

Caveat: run_enrichment queries Enrichr live and fetch_abstract queries PubMed; Enrichr's
gene-set libraries change over time, so enrichment results may differ from the original
run (the original inputs and results are saved in analysis/enrichment_*.{txt,csv}).

  uv run python examples/GSE246386/analysis/replay.py
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


SESSION.begin_run('20260930_GSE246386_full_replay')
print(f"Replay output: {SESSION.paths.dir}")

_run(
    "scan_results",
    results_dir='examples/GSE246386/results',
)

_run(
    "fetch_abstract",
    pmid='38961746',
)

_run(
    "load_counts",
    counts_path='examples/GSE246386/results/salmon/salmon.merged.gene_counts.tsv',
    design_path='examples/GSE246386/design.csv',
)

_run("compute_qc")

_run("read_multiqc")

_run("filter_low_counts")

_run(
    "run_deseq2",
    contrast=['condition', 'GFI1B', 'EV'],
)

_run(
    "run_enrichment",
    direction='up',
    organism='human',
)

_run(
    "run_enrichment",
    direction='down',
    organism='human',
)

_run("generate_figures")

_run(
    "write_report",
    report_markdown=('# GFI1B ectopic expression in iPSC-derived hemogenic endothelium (GSE246386)\n'
     '\n'
     '## 1. Executive Summary\n'
     'This analysis compared GFI1B-overexpressing CD34+ iPSC-derived hemogenic endothelium '
     'with empty-vector (EV) controls using Salmon gene-level counts. {{de.n_significant}} '
     'of {{de.n_tested}} tested genes were significant ({{de.pct_significant}}). The '
     'hematopoietic regulators named in the associated study ({{cite:38961746}}) rose, and '
     'endothelial and extracellular-matrix programmes fell. Upregulated genes were also '
     'enriched for cell-cycle and interferon-response terms.\n'
     '\n'
     '## 2. Experimental Design\n'
     'The dataset is GSE246386: human (GRCh38) CD34+CD43-CD73- iPSC-derived hemogenic '
     'endothelium, sampled at the day-four stage according to the GEO sample titles. The '
     'design has {{samples.n}} samples, with conditions {{samples.conditions}} '
     '({{samples.n_ev}} EV and {{samples.n_gfi1b}} GFI1B replicates). The contrast is '
     '{{de.contrast}}, with EV as the reference. The design file has no other varying '
     'factor, so no covariates were used. The counts come from {{provenance.data_source}} '
     '({{provenance.quantification}}).\n'
     '\n'
     '## 3. Quality Control\n'
     '{{table:qc}}\n'
     '\n'
     'Library sizes ranged from {{qc.library_size_min_millions}} to '
     '{{qc.library_size_max_millions}} million (median {{qc.library_size_median_millions}} '
     'million). Between {{qc.genes_detected_min}} and {{qc.genes_detected_max}} genes were '
     'detected per sample. Salmon mapping rates were '
     '{{multiqc.salmon_percent_mapped.min}}–{{multiqc.salmon_percent_mapped.max}}%. Raw '
     'FastQC duplication was '
     '{{multiqc.fastqc_raw_percent_duplicates.min}}–{{multiqc.fastqc_raw_percent_duplicates.max}}%.\n'
     '\n'
     '![Library sizes](figures/library_sizes.png)\n'
     '\n'
     '{{qc.pca_summary}} PC1 explains {{qc.pca_pc1_pct}} of variance and PC2 '
     '{{qc.pca_pc2_pct}}. The EV replicates are much more dispersed than the GFI1B '
     'replicates (spread {{qc.pca_spread.ev}} vs {{qc.pca_spread.gfi1b}}). No sample '
     'clusters with the other group.\n'
     '\n'
     '![PCA](figures/pca.png)\n'
     '\n'
     '![PC-metadata association](figures/pc_association.png)\n'
     '\n'
     'Pairwise sample correlations ranged from {{qc.sample_correlation_min}} to '
     '{{qc.sample_correlation_max}}.\n'
     '\n'
     '![Sample correlation](figures/sample_correlation.png)\n'
     '\n'
     "No sample was flagged as an outlier by the PCA facts. The EV group's uneven spread "
     'is noted under Limitations.\n'
     '\n'
     '## 4. Differential Expression\n'
     '{{table:de_summary}}\n'
     '\n'
     '![Volcano plot](figures/volcano.png)\n'
     '\n'
     '![MA plot](figures/ma_plot.png)\n'
     '\n'
     'Top upregulated genes:\n'
     '\n'
     '{{table:top_up}}\n'
     '\n'
     'Top downregulated genes:\n'
     '\n'
     '{{table:top_down}}\n'
     '\n'
     '![DE heatmap](figures/de_heatmap.png)\n'
     '\n'
     '{{gene:GFI1B}} is the most strongly induced gene, as expected for ectopic '
     'expression. {{gene:ALOX15}} is also strongly up.\n'
     '\n'
     '## 5. Gene Set Enrichment\n'
     'Upregulated genes:\n'
     '\n'
     '{{table:enrichment_up}}\n'
     '\n'
     '![GO BP up](figures/enrichment_upregulated_go_biological_process_2023.png)\n'
     '\n'
     '![KEGG up](figures/enrichment_upregulated_kegg_2021_human.png)\n'
     '\n'
     'Downregulated genes:\n'
     '\n'
     '{{table:enrichment_down}}\n'
     '\n'
     '![GO BP down](figures/enrichment_downregulated_go_biological_process_2023.png)\n'
     '\n'
     '![KEGG down](figures/enrichment_downregulated_kegg_2021_human.png)\n'
     '\n'
     '## 6. Comparison with the published study\n'
     '{{table:paper_genes}}\n'
     '\n'
     'The abstract ({{cite:38961746}}) reports that ectopic GFI1B expression downregulated '
     'endothelial genes and upregulated hematopoietic genes, including GATA2, KIT, RUNX1 '
     'and SPI1. In this analysis {{gene:GATA2}}, {{gene:KIT}}, {{gene:RUNX1}} and '
     '{{gene:SPI1}} are all significantly up, which agrees with the abstract. '
     '{{gene:KDM1A}} is not significantly changed. The abstract describes the '
     'hematopoietic specification as "partial". Endothelial markers {{gene:CDH5}} and '
     '{{gene:PECAM1}} fall only modestly, and {{gene:KDR}} is not significant.\n'
     '\n'
     '## 7. Biological Interpretation\n'
     '- The hematopoietic regulators {{gene:GATA2}}, {{gene:RUNX1}} and {{gene:KIT}} are '
     'induced, consistent with {{cite:38961746}}.\n'
     '- Downregulated genes are enriched for {{term:down:Sprouting Angiogenesis '
     '(GO:0002040)}} and {{term:down:Endothelial Cell Proliferation (GO:0001935)}}. They '
     'are also enriched for {{term:down:Extracellular Matrix Organization (GO:0030198)}}. '
     'The endothelial genes {{gene:DLL4}}, {{gene:HEY1}} and {{gene:TEK}} are all lower. '
     'This fits reduced endothelial identity.\n'
     '- Upregulated genes show strong enrichment for {{term:up:Cell cycle}} and '
     '{{term:up:DNA replication}}. The replication-licensing family is also shifted: '
     '{{genes:MCM*}}. GFI1B overexpression therefore coincides with a proliferative '
     'signature.\n'
     '- An interferon-response signature is also enriched: {{term:up:Defense Response To '
     'Virus (GO:0051607)}}. {{gene:IFIT1}} is up. This design does not show whether the '
     'signature reflects GFI1B, the vector or transduction.\n'
     '- Collagen genes do not move uniformly: {{genes:COL*}}. For example {{gene:COL5A2}} '
     'is down while {{gene:COL1A1}} is up, so the ECM signal is not a simple global '
     'decrease.\n'
     '\n'
     '## 8. Limitations\n'
     '- There are only three replicates per group. EV replicates are much more dispersed '
     'in PCA than GFI1B replicates, which may affect variance estimates.\n'
     '- Only a single stage is sampled. The analysis does not separate direct GFI1B '
     'targets from secondary effects.\n'
     '- Enrichment used the top-ranked genes per direction, capped at the maximum given in '
     'Methods. More genes passed the thresholds than were tested, so lower-ranked genes '
     'were not included. The enrichment background is the Enrichr default.\n'
     '- Salmon counts are fractional and were rounded before DE. Alignment was skipped, so '
     'there are no alignment-based QC metrics.\n'
     '- The data come from a public GEO series, and no details of the expression vector '
     'beyond the GEO metadata are recorded here.\n'
     '\n'
     '## 9. Methods\n'
     '{{table:methods}}\n'
     '\n'
     '## 10. References\n'
     '{{cite:38961746}}\n'),
)

print(f"Done. Report: {SESSION.paths.analysis_report}")

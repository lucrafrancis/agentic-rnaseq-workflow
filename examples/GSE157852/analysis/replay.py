"""Replay of the analysis agent's tool calls for run 20260930_GSE157852_counts_2.

Generated: 2026-09-30T11:30:56
Git commit: 3e100049331194e5b3085014b342647d00cf49d1-dirty
Source log: tool_calls.jsonl (10 of 19 calls; failed and read-only calls omitted)
Copied to examples/GSE157852/ for the repository: data paths were changed from runs/<run>/ to
examples/GSE157852/ (the only edit; the tool logs are unchanged).

Re-runs the same tool functions with the same arguments, without the LLM. Output goes
to a new run directory (runs/<date>_20260930_GSE157852_counts_2_replay/). Package versions are pinned by
uv.lock at the commit above — check it out first if the code has changed since.

Caveat: run_enrichment queries Enrichr live and fetch_abstract queries PubMed; Enrichr's
gene-set libraries change over time, so enrichment results may differ from the original
run (the original inputs and results are saved in analysis/enrichment_*.{txt,csv}).

  uv run python examples/GSE157852/analysis/replay.py
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


SESSION.begin_run('20260930_GSE157852_counts_2_replay')
print(f"Replay output: {SESSION.paths.dir}")

_run(
    "scan_results",
    results_dir='examples/GSE157852',
)

_run(
    "fetch_abstract",
    pmid='33010822',
)

_run(
    "load_counts",
    counts_path='examples/GSE157852/counts.tsv',
    design_path='examples/GSE157852/design.csv',
)

_run("compute_qc")

_run("filter_low_counts")

_run(
    "run_deseq2",
    contrast=['condition', 'SARS-CoV-2 72 hpi', 'Mock 72 hpi'],
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
    report_markdown=('# SARS-CoV-2 infection of human choroid plexus organoids: {{de.contrast}} '
     '(GSE157852)\n'
     '\n'
     '## 1. Executive Summary\n'
     'This analysis compared SARS-CoV-2-infected and mock-infected hiPSC-derived choroid '
     "plexus organoids using the authors' {{provenance.value_type}} from GEO "
     '({{provenance.accession}}). {{de.n_significant}} of {{de.n_tested}} tested genes '
     'were differentially expressed. Viral transcripts were strongly detected in infected '
     'samples. The host response was up-regulation of migration, extracellular matrix, '
     'focal adhesion and cytokine/chemokine genes, and down-regulation of ion-transport '
     'and secretory genes. The design has few replicates per group, so the results should '
     'be read with that in mind.\n'
     '\n'
     '## 2. Experimental Design\n'
     '- Dataset: {{provenance.accession}}, {{provenance.source}} ({{provenance.file}}); '
     'data source: {{provenance.data_source}}.\n'
     '- Model system: {{design.cell_type}}, Homo sapiens. The associated study is '
     '{{cite:33010822}}.\n'
     '- Conditions: {{samples.conditions}}, with {{samples.n_mock_72_hpi}} Mock and '
     '{{samples.n_sars_cov_2_72_hpi}} SARS-CoV-2 replicates ({{samples.n}} samples in '
     'total).\n'
     '- The contrast was {{de.contrast}}, with Mock as the reference. The GEO series also '
     'contains earlier-time-point infected samples, which are not in the supplied matrix '
     'and were not analysed.\n'
     '\n'
     '## 3. Quality Control\n'
     '{{table:qc}}\n'
     '\n'
     'Library sizes ranged from {{qc.library_size_min_millions}} to '
     '{{qc.library_size_max_millions}} million reads (median '
     '{{qc.library_size_median_millions}} million). Detected genes per sample ranged from '
     '{{qc.genes_detected_min}} to {{qc.genes_detected_max}}.\n'
     '\n'
     '![Library sizes](figures/library_sizes.png)\n'
     '\n'
     'PC1 and PC2 explain {{qc.pca_pc1_pct}} and {{qc.pca_pc2_pct}} of variance. '
     '{{qc.pca_summary}} The infected group is therefore noticeably more spread than the '
     'mock group.\n'
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
     'No sample was flagged as misclustered ({{qc.pca_misclustered}}), and no sample was '
     'removed.\n'
     '\n'
     '## 4. Differential Expression\n'
     '{{table:de_summary}}\n'
     '\n'
     'The volcano and MA plots show {{de.pct_up}} of tested genes up and {{de.pct_down}} '
     'down at the significance threshold.\n'
     '\n'
     '![Volcano plot](figures/volcano.png)\n'
     '\n'
     '![MA plot](figures/ma_plot.png)\n'
     '\n'
     '### Top up-regulated genes\n'
     '{{table:top_up}}\n'
     '\n'
     'The top of the list is dominated by SARS-CoV-2 genomic and subgenomic transcripts, '
     'for example {{gene:N}}, {{gene:S}} and {{gene:ORF3a}}. These reflect viral reads in '
     'infected samples, since the Mock samples have essentially none. Among host genes, '
     '{{gene:TM4SF1}}, {{gene:CTGF}} and {{gene:CCL2}} are strongly up.\n'
     '\n'
     '### Top down-regulated genes\n'
     '{{table:top_down}}\n'
     '\n'
     '![DE heatmap](figures/de_heatmap.png)\n'
     '\n'
     'Other results of interest:\n'
     '- The choroid plexus marker {{gene:TTR}} is modestly decreased.\n'
     '- {{gene:ACE2}} and {{gene:TMPRSS2}} are not significantly changed.\n'
     '- {{gene:IFIT1}} is down, so there is no sign of a canonical interferon-stimulated '
     'gene induction. It was the only interferon-stimulated gene I queried that passed '
     'filtering.\n'
     '- {{genes:CLDN*}} claudins are mixed: {{gene:CLDN5}}, {{gene:CLDN2}} and '
     '{{gene:CLDN16}} are down, while {{gene:CLDN4}} and {{gene:CLDN6}} are up.\n'
     '\n'
     '## 5. Gene Set Enrichment\n'
     '### Up-regulated genes\n'
     '{{table:enrichment_up}}\n'
     '\n'
     '![GO BP up](figures/enrichment_upregulated_go_biological_process_2023.png)\n'
     '\n'
     '![KEGG up](figures/enrichment_upregulated_kegg_2021_human.png)\n'
     '\n'
     '### Down-regulated genes\n'
     '{{table:enrichment_down}}\n'
     '\n'
     '![GO BP down](figures/enrichment_downregulated_go_biological_process_2023.png)\n'
     '\n'
     '![KEGG down](figures/enrichment_downregulated_kegg_2021_human.png)\n'
     '\n'
     'The up-regulated set is enriched for cell migration, wound healing, extracellular '
     'matrix and adhesion terms, for example {{term:up:Regulation Of Cell Migration '
     '(GO:0030334)}}, {{term:up:Extracellular Matrix Organization (GO:0030198)}} and '
     '{{term:up:Focal adhesion}}. The apoptosis term {{term:up:Regulation Of Apoptotic '
     'Process (GO:0042981)}} is also enriched. Among the down-regulated genes, '
     '{{term:down:Metal Ion Transport (GO:0030001)}} and {{term:down:Bile secretion}} are '
     'the terms that pass the significance threshold. Most other down-regulated terms do '
     'not reach it.\n'
     '\n'
     '## 6. Comparison with the published study\n'
     'The abstract of {{cite:33010822}} names no specific genes, so no paper-gene table is '
     'shown. The abstract reports increased cell death and transcriptional dysregulation '
     'indicative of an inflammatory response and cellular function deficits after '
     'infection. Here the apoptosis-related term {{term:up:Regulation Of Apoptotic Process '
     '(GO:0042981)}} and the inflammatory term {{term:up:Inflammatory Response '
     '(GO:0006954)}} are enriched among up-regulated genes, which is compatible with that '
     'description. The inflammatory term is only modestly significant. The down-regulation '
     'of ion transport and secretory genes could correspond to the reported function '
     'deficits, but the abstract does not say so.\n'
     '\n'
     '## 7. Biological Interpretation\n'
     '- Infection is clearly present, with very large increases in viral transcripts such '
     'as {{gene:N}} and {{gene:ORF8}}.\n'
     '- The host response is dominated by extracellular matrix, adhesion and migration '
     'programmes ({{term:up:ECM-receptor interaction}}). The chemokine {{gene:CCL2}} and '
     'the adhesion molecule {{gene:ICAM1}} are also up.\n'
     '- Apoptosis-related genes go up, for example {{gene:CASP8}}. This is consistent with '
     'the cell death described in {{cite:33010822}}.\n'
     '- Down-regulated genes point to loss of transport and secretory functions '
     '({{term:down:Bile secretion}}). {{gene:AQP1}}, {{gene:ATP1A2}} and {{gene:SLC4A4}} '
     'are all down.\n'
     '- An interferon response is not evident: {{gene:IFIT1}} falls, and the '
     'interferon-related terms searched were not significantly enriched '
     '({{term:up:Response To Type II Interferon (GO:0034341)}}).\n'
     '\n'
     '## 8. Limitations\n'
     '- There are only {{samples.n_sars_cov_2_72_hpi}} replicates per group, and the '
     'infected group is more spread in PCA than the mock group '
     '({{qc.pca_spread.sars_cov_2_72_hpi}} versus {{qc.pca_spread.mock_72_hpi}}).\n'
     '- The data are bulk RNA-seq of whole organoids, so changes in cell composition (for '
     'example cell loss) cannot be separated from changes in expression within cells.\n'
     '- Enrichment used Enrichr over-representation. The up-regulated list was capped at '
     'the maximum number of genes, so {{enrichment.up.n_passing}} genes passed the '
     'thresholds but {{enrichment.up.n_input_genes}} were used.\n'
     '- Several down-regulated terms have padj above the significance threshold, so they '
     'are suggestive only.\n'
     '- Viral genes appear in the count matrix and were included in the DE analysis and '
     'the top-gene lists.\n'
     '- Only the single contrast requested was performed. Other time points in the GEO '
     'series were not analysed.\n'
     '\n'
     '## 9. Methods\n'
     '{{table:methods}}\n'
     '\n'
     '## 10. References\n'
     '{{cite:33010822}}\n'),
)

print(f"Done. Report: {SESSION.paths.analysis_report}")

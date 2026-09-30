"""Replay of the analysis agent's tool calls for run 20260930_GSE157852_counts.

Generated: 2026-09-30T11:00:30
Git commit: dbce4ddf4fd2dad866a3009cd89a5a0a8791fc3d-dirty
Source log: tool_calls.jsonl (10 of 18 calls; failed and read-only calls omitted)
Copied to examples/GSE157852/ for the repository: data paths were changed from runs/<run>/ to
examples/GSE157852/ (the only edit; the tool logs are unchanged).

Re-runs the same tool functions with the same arguments, without the LLM. Output goes
to a new run directory (runs/<date>_20260930_GSE157852_counts_replay/). Package versions are pinned by
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


SESSION.begin_run('20260930_GSE157852_counts_replay')
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
     'This analysis performed differential expression and gene set enrichment on a '
     'processed GEO count matrix ({{provenance.accession}}) comparing SARS-CoV-2-infected '
     'human choroid plexus organoids (CPOs) to mock-infected controls at the matched '
     'post-infection time point ({{de.contrast}}). Of {{de.n_tested}} genes tested, '
     '{{de.n_significant}} were significant ({{de.pct_significant}} of tested genes; '
     '{{de.n_up}} up, {{de.n_down}} down), with viral transcripts among the most strongly '
     'enriched genes, confirming productive infection. Enrichment of the upregulated gene '
     'set highlighted cell migration/adhesion, extracellular matrix, coagulation, and '
     'inflammatory response and apoptosis pathways, broadly consistent with the '
     'transcriptional dysregulation and cell-death phenotype reported by '
     '{{cite:33010822}}.\n'
     '\n'
     '## 2. Experimental Design\n'
     '\n'
     'Samples are {{design.cell_type}} ({{provenance.data_source}}, {{provenance.file}}), '
     'organism Homo sapiens. This analysis compares {{samples.n_sars_cov_2_72_hpi}} '
     'SARS-CoV-2-infected replicates against {{samples.n_mock_72_hpi}} mock replicates at '
     'the same time point ({{de.contrast}}). Data were provided as '
     '{{provenance.value_type}} with gene identifiers as {{provenance.gene_id_type}}.\n'
     '\n'
     '## 3. Quality Control\n'
     '\n'
     '{{table:qc}}\n'
     '\n'
     'Library sizes ranged from {{qc.library_size_min_millions}} to '
     '{{qc.library_size_max_millions}} million reads (median '
     '{{qc.library_size_median_millions}} million); gene detection ranged from '
     '{{qc.genes_detected_min}} to {{qc.genes_detected_max}} genes per sample. Pairwise '
     'sample correlations ranged from {{qc.sample_correlation_min}} to '
     '{{qc.sample_correlation_max}}.\n'
     '\n'
     '![Library sizes per sample before filtering](figures/library_sizes.png)\n'
     '\n'
     '![PCA of samples coloured by condition](figures/pca.png)\n'
     '\n'
     '{{qc.pca_summary}} Condition explains {{qc.pca_condition_r2_pc1}} of PC1 variance '
     'and {{qc.pca_condition_r2_pc2}} of PC2 variance, and no sample was misclustered '
     '({{qc.pca_misclustered}}). The SARS-CoV-2 group shows notably more within-group '
     'spread than the Mock group ({{qc.pca_spread.sars_cov_2_72_hpi}} vs '
     '{{qc.pca_spread.mock_72_hpi}}), consistent with variable infection/response '
     'magnitude across replicates.\n'
     '\n'
     '![PC association with design variables](figures/pc_association.png)\n'
     '\n'
     '![Sample-sample correlation heatmap](figures/sample_correlation.png)\n'
     '\n'
     'No samples were excluded as outliers; all six were retained for differential '
     'expression.\n'
     '\n'
     '## 4. Differential Expression\n'
     '\n'
     'Low-count genes were filtered prior to DE ({{filter.pct_removed}} of genes removed, '
     'leaving {{filter.genes_after}} of {{filter.genes_before}}). DESeq2 was run with '
     'design {{de.design}}, contrast {{de.contrast}}.\n'
     '\n'
     '{{table:de_summary}}\n'
     '\n'
     '![Volcano plot of differential expression](figures/volcano.png)\n'
     '\n'
     '![MA plot of differential expression](figures/ma_plot.png)\n'
     '\n'
     'The top upregulated genes include several SARS-CoV-2 viral transcripts, e.g. '
     '{{gene:N}}, {{gene:ORF6}} and {{gene:S}}, reflecting direct detection of viral RNA '
     'in infected samples rather than host transcriptional change. Among host genes, '
     '{{gene:CTGF}} and {{gene:TM4SF1}} are strongly upregulated. The top downregulated '
     'gene is {{gene:UGT2A1}}.\n'
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
     'Upregulated genes are enriched for {{term:up:Regulation Of Cell Migration '
     '(GO:0030334)}}, {{term:up:Extracellular Matrix Organization (GO:0030198)}}, '
     '{{term:up:Focal adhesion}} and {{term:up:Complement and coagulation cascades}}, '
     'together with {{term:up:Inflammatory Response (GO:0006954)}} and '
     '{{term:up:Apoptosis}}.\n'
     '\n'
     '{{table:enrichment_down}}\n'
     '\n'
     '![GO Biological Process enrichment, downregulated '
     'genes](figures/enrichment_downregulated_go_biological_process_2023.png)\n'
     '\n'
     '![KEGG enrichment, downregulated '
     'genes](figures/enrichment_downregulated_kegg_2021_human.png)\n'
     '\n'
     'Downregulated genes are enriched for {{term:down:Metal Ion Transport (GO:0030001)}} '
     'and {{term:down:Bile secretion}}/{{term:down:Pancreatic secretion}} KEGG terms, '
     'which largely reflect ion-transport and secretory-epithelium functions rather than '
     'classic immune pathways.\n'
     '\n'
     '## 6. Comparison with the Published Study\n'
     '\n'
     '{{cite:33010822}} reports that SARS-CoV-2 infection of hiPSC-derived choroid plexus '
     'organoids caused increased cell death and "transcriptional dysregulation indicative '
     'of an inflammatory response and cellular function deficits." This analysis is '
     'consistent with that description at a pathway level: the upregulated set shows '
     'significant enrichment of {{term:up:Inflammatory Response (GO:0006954)}} and '
     '{{term:up:Apoptosis}}, and the downregulated set shows loss of '
     'ion-transport/secretory terms such as {{term:down:Metal Ion Transport '
     '(GO:0030001)}}, consistent with a "cellular function deficit" in a secretory '
     'epithelium. The abstract does not name specific genes, so no gene-level comparison '
     'table is included; the agreement here is at the level of the general pathway themes '
     'it describes.\n'
     '\n'
     '## 7. Biological Interpretation\n'
     '\n'
     '- Viral transcripts ({{gene:N}}, {{gene:ORF6}}, {{gene:S}}, among others) dominate '
     'the largest log2 fold-changes, confirming that the SARS-CoV-2 samples were '
     'productively infected, as expected from the study design.\n'
     '- The upregulated host response is dominated by {{term:up:Focal adhesion}}, '
     '{{term:up:Extracellular Matrix Organization (GO:0030198)}} and {{term:up:Regulation '
     'Of Cell Migration (GO:0030334)}} terms, alongside a smaller but significant '
     '{{term:up:Inflammatory Response (GO:0006954)}} signature, indicating that ECM '
     'remodeling/adhesion changes are the largest-magnitude transcriptional theme, with '
     'inflammation present but less dominant by gene count.\n'
     '- {{term:up:Apoptosis}} enrichment among upregulated genes is consistent with the '
     'increased cell death reported in the source study, though this analysis only shows '
     'transcriptional correlates, not direct cell-death measurements.\n'
     '- Downregulated genes are enriched for {{term:down:Metal Ion Transport '
     '(GO:0030001)}}, pointing to reduced expression of ion-channel/transporter genes that '
     'are not part of the inflammatory narrative and represent a result outside the '
     "paper's stated focus.\n"
     '- {{genes:ITG*}}, several of which appear among the focal-adhesion/ECM upregulated '
     'genes, illustrate a broad integrin-pathway shift accompanying infection.\n'
     '\n'
     '## 8. Limitations\n'
     '\n'
     '- Each condition has only {{samples.n_mock_72_hpi}} (Mock) and '
     '{{samples.n_sars_cov_2_72_hpi}} (SARS-CoV-2) replicates, limiting statistical power '
     'and the ability to detect subtler effects.\n'
     '- The SARS-CoV-2 group shows higher PCA spread than the Mock group '
     '({{qc.pca_spread.sars_cov_2_72_hpi}} vs {{qc.pca_spread.mock_72_hpi}}), indicating '
     'more variable infection response across replicates that should be considered when '
     'interpreting individual gene calls.\n'
     '- Enrichment analysis used the Enrichr default background rather than a '
     'study-specific background, and the lowest-magnitude viral-gene fold changes are not '
     'directly comparable to host gene changes because viral transcripts are absent in '
     'Mock samples by design.\n'
     "- This is a re-analysis of the authors' processed count matrix ({{provenance.file}}) "
     'rather than a re-alignment from FASTQ, so any upstream processing choices made by '
     'the original authors are inherited as-is.\n'
     '- No paired/batch covariate was available beyond condition; the design used is '
     '{{de.design}}.\n'
     '\n'
     '## Methods\n'
     '\n'
     '{{table:methods}}\n'
     '\n'
     '## References\n'
     '\n'
     '{{cite:33010822}}\n'),
)

print(f"Done. Report: {SESSION.paths.analysis_report}")

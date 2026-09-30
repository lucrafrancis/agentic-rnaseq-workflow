"""Replay of the analysis agent's tool calls for run 20260930_GSE164073_counts_5.

Generated: 2026-09-30T11:28:54
Git commit: 3e100049331194e5b3085014b342647d00cf49d1-dirty
Source log: tool_calls.jsonl (10 of 19 calls; failed and read-only calls omitted)
Copied to examples/GSE164073/ for the repository: data paths were changed from runs/<run>/ to
examples/GSE164073/ (the only edit; the tool logs are unchanged).

Re-runs the same tool functions with the same arguments, without the LLM. Output goes
to a new run directory (runs/<date>_20260930_GSE164073_counts_5_replay/). Package versions are pinned by
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


SESSION.begin_run('20260930_GSE164073_counts_5_replay')
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
    report_markdown=('# SARS-CoV-2 infection of human ocular surface tissue (GSE164073)\n'
     '\n'
     '## Executive Summary\n'
     'We compared SARS-CoV-2-infected against mock-treated ocular surface tissue from '
     "{{provenance.accession}}, adjusting for tissue of origin, using the authors' raw "
     'count matrix. Infection was associated with {{de.n_significant}} differentially '
     'expressed genes ({{de.n_up}} up, {{de.n_down}} down). The up-regulated set was '
     'dominated by NF-κB/TNF-linked chemokine and inflammatory signatures. The '
     'down-regulated set was enriched for epithelial differentiation and keratin genes.\n'
     '\n'
     '## Experimental Design\n'
     'The dataset is {{provenance.data_source}} ({{provenance.source}}, '
     '{{provenance.file}}), {{provenance.value_type}} for human (Homo sapiens) cornea, '
     'limbus and sclera from cadaver donors ({{cite:34022129}}). Tissues: '
     '{{design.tissue_values}}. Conditions: {{samples.conditions}}; replicates per '
     'condition: mock {{samples.n_mock}}, infected {{samples.n_sars_cov_2_moi_1_0}} '
     '({{samples.n}} samples in total). Time point: {{design.time_point}}. The contrast is '
     '{{de.contrast}}, with the design {{de.design}}.\n'
     '\n'
     '## Quality Control\n'
     '{{table:qc}}\n'
     '\n'
     'Library sizes are shown below (median {{qc.library_size_median_millions}} million '
     'reads; genes detected per sample ranged from {{qc.genes_detected_min}} to '
     '{{qc.genes_detected_max}}).\n'
     '\n'
     '![Library sizes per sample](figures/library_sizes.png)\n'
     '\n'
     'PC1 and PC2 explain {{qc.pca_pc1_pct}} and {{qc.pca_pc2_pct}} of variance. '
     '{{qc.pca_summary}}\n'
     '\n'
     '![PCA coloured by condition](figures/pca.png)\n'
     '\n'
     '![PCA coloured by tissue](figures/pca_tissue.png)\n'
     '\n'
     'Tissue, not infection, dominates the leading principal components, which is why '
     'tissue is included as a covariate in the model. Between-sample correlations ranged '
     'from {{qc.sample_correlation_min}} to {{qc.sample_correlation_max}}.\n'
     '\n'
     '![Sample correlation](figures/sample_correlation.png)\n'
     '\n'
     '![PC–metadata associations](figures/pc_association.png)\n'
     '\n'
     "The PCA flags {{qc.pca_misclustered}} as closer to the other condition's centroid "
     'within their tissue. No samples were removed.\n'
     '\n'
     '## Differential Expression\n'
     'After filtering, {{filter.genes_after}} genes remained ({{filter.pct_removed}} '
     'removed). {{de.n_significant}} genes ({{de.pct_significant}}) had padj below '
     '{{de.padj_threshold}}: {{de.n_up}} up ({{de.pct_up}}) and {{de.n_down}} down '
     '({{de.pct_down}}).\n'
     '\n'
     '{{table:de_summary}}\n'
     '\n'
     '![Volcano plot](figures/volcano.png)\n'
     '\n'
     '![MA plot](figures/ma_plot.png)\n'
     '\n'
     '{{table:top_up}}\n'
     '\n'
     '{{table:top_down}}\n'
     '\n'
     'The most significant up-regulated genes include {{gene:SOD2}}, {{gene:TNFAIP3}} and '
     '{{gene:C3}}. The top down-regulated genes include {{gene:ACAN}}, {{gene:ACTC1}} and '
     '{{gene:CTGF}}.\n'
     '\n'
     '![DE heatmap](figures/de_heatmap.png)\n'
     '\n'
     '## Gene Set Enrichment\n'
     'Enrichment used the genes passing the thresholds described in Methods '
     '({{enrichment.up.n_input_genes}} up, {{enrichment.down.n_input_genes}} down).\n'
     '\n'
     '{{table:enrichment_up}}\n'
     '\n'
     '![Up GO BP](figures/enrichment_upregulated_go_biological_process_2023.png)\n'
     '\n'
     '![Up KEGG](figures/enrichment_upregulated_kegg_2021_human.png)\n'
     '\n'
     'Up-regulated genes were enriched for {{term:up:Inflammatory Response (GO:0006954)}} '
     'and {{term:up:Chemokine-Mediated Signaling Pathway (GO:0070098)}}. KEGG terms '
     'include {{term:up:TNF signaling pathway}} and {{term:up:NF-kappa B signaling '
     'pathway}}.\n'
     '\n'
     '{{table:enrichment_down}}\n'
     '\n'
     '![Down GO BP](figures/enrichment_downregulated_go_biological_process_2023.png)\n'
     '\n'
     '![Down KEGG](figures/enrichment_downregulated_kegg_2021_human.png)\n'
     '\n'
     'Down-regulated genes were enriched for {{term:down:Epithelium Development '
     '(GO:0060429)}} and {{term:down:Keratinocyte Differentiation (GO:0030216)}}. The '
     'down-regulated KEGG terms are weaker; see the table for which reach the significance '
     'threshold.\n'
     '\n'
     '## Comparison with the published study\n'
     "No table of the paper's named genes was generated for this run, so the comparison is "
     'qualitative and based on the abstract alone. The abstract ({{cite:34022129}}) '
     'reports robust NF-κB induction in infected cells and diminished type I/III '
     'interferon signaling. NF-κB induction agrees with our results: {{term:up:NF-kappa B '
     'signaling pathway}} is enriched. Separately, {{gene:RELB}}, {{gene:NFKBIA}} and '
     '{{gene:NFKB2}} are up. The interferon part is not clearly supported here. '
     '{{gene:MX1}} and {{gene:STAT1}} are up, while {{gene:IFIT1}}, {{gene:ISG15}} and '
     '{{gene:OAS1}} are not significant. The interferon genes IFNB1 and IFNL1 were removed '
     'by low-count filtering, so this analysis cannot assess their expression directly. '
     'Our pooled, tissue-adjusted model does not test whether interferon signaling is '
     'diminished.\n'
     '\n'
     '## Biological Interpretation\n'
     '- Infection is accompanied by a chemokine and inflammatory response: '
     '{{genes:CXCL*}}, plus {{gene:IL6}}, which is also up. {{term:up:Neutrophil '
     'Chemotaxis (GO:0030593)}} is enriched among up-regulated genes.\n'
     '- NF-κB pathway components are induced ({{gene:RELB}}, {{gene:TNFAIP3}}), consistent '
     'with the NF-κB signature described in {{cite:34022129}}.\n'
     '- Down-regulated genes point to reduced epithelial differentiation and keratin '
     'expression ({{term:down:Epidermal Cell Differentiation (GO:0009913)}}). The data '
     'cannot distinguish epithelial damage from an altered cell state.\n'
     '- {{gene:ACE2}} is not significantly changed, with very low expression in this '
     'dataset.\n'
     '\n'
     '## Limitations\n'
     '- Replication is limited to a few replicates per tissue and condition. '
     'Tissue-specific infection responses, which the abstract suggests differ by region, '
     'are not tested here because the model estimates one pooled infection effect adjusted '
     'for tissue.\n'
     '- The data are from cadaver donor tissue. Donor identity is not in the design, so '
     'pairing could not be modelled.\n'
     '- Enrichment used only strongly changed genes (see Methods). Enrichr uses its '
     'default library background rather than the expressed-gene universe. '
     'Interferon-related genes with low counts were filtered out.\n'
     "- The data are the authors' processed count matrix, so alignment and quantification "
     'were not assessed here.\n'
     '- Viral read content was not examined.\n'
     '\n'
     '## Methods\n'
     '{{table:methods}}\n'
     '\n'
     '## References\n'
     '- {{cite:34022129}}\n'),
)

print(f"Done. Report: {SESSION.paths.analysis_report}")

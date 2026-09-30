## 1. Executive Summary

This analysis performed differential expression and gene set enrichment on a processed GEO count matrix (GSE157852) comparing SARS-CoV-2-infected human choroid plexus organoids (CPOs) to mock-infected controls at the matched post-infection time point (SARS-CoV-2 72 hpi vs Mock 72 hpi). Of 15,308 genes tested, 3,212 were significant (21.0% of tested genes; 1,688 up, 1,524 down), with viral transcripts among the most strongly enriched genes, confirming productive infection. Enrichment of the upregulated gene set highlighted cell migration/adhesion, extracellular matrix, coagulation, and inflammatory response and apoptosis pathways, broadly consistent with the transcriptional dysregulation and cell-death phenotype reported by (Jacob et al., 2020; PMID: 33010822).

## 2. Experimental Design

Samples are hiPSC-derived choroid plexus organoids (CPOs) (GEO count matrix, GSE157852_CPO_RawCounts.txt.gz), organism Homo sapiens. This analysis compares 3 SARS-CoV-2-infected replicates against 3 mock replicates at the same time point (SARS-CoV-2 72 hpi vs Mock 72 hpi). Data were provided as raw integer counts with gene identifiers as symbol.

## 3. Quality Control

| Sample | Condition | Total counts | Genes detected |
|---|---|---|---|
| Mock_72_hpi_rep1 | Mock 72 hpi | 23,507,394 | 20,894 |
| Mock_72_hpi_rep2 | Mock 72 hpi | 22,884,961 | 20,393 |
| Mock_72_hpi_rep3 | Mock 72 hpi | 20,121,047 | 19,643 |
| SARS_CoV_2_72_hpi_rep1 | SARS-CoV-2 72 hpi | 31,305,936 | 21,330 |
| SARS_CoV_2_72_hpi_rep2 | SARS-CoV-2 72 hpi | 16,889,301 | 20,072 |
| SARS_CoV_2_72_hpi_rep3 | SARS-CoV-2 72 hpi | 23,628,250 | 19,946 |

*Genes detected: genes with at least one count, before low-count filtering.*

Library sizes ranged from 16.9 to 31.3 million reads (median 23.2 million); gene detection ranged from 19,643 to 21,330 genes per sample. Pairwise sample correlations ranged from 0.95 to 0.99.

![Library sizes per sample before filtering](figures/library_sizes.png)

*File: `figures/library_sizes.png` — Total counts per sample before low-count filtering (6 samples, 16.9–31.3 million), coloured by condition.*

![PCA of samples coloured by condition](figures/pca.png)

*File: `figures/pca.png` — PCA of log2(counts + 1): PC1 44.5%, PC2 31.2% of variance, coloured by condition. On PC1–PC2, mean distance of replicates to their group centroid: Mock 72 hpi 20.4, SARS-CoV-2 72 hpi 45.8; group centroids are 87.5 apart. The most spread group is 2.2× the least spread. Condition explains 89% of PC1 and 8% of PC2 variance. Samples closer to another group's centroid: none.*

On PC1–PC2, mean distance of replicates to their group centroid: Mock 72 hpi 20.4, SARS-CoV-2 72 hpi 45.8; group centroids are 87.5 apart. The most spread group is 2.2× the least spread. Condition explains 89% of PC1 and 8% of PC2 variance. Samples closer to another group's centroid: none. Condition explains 89% of PC1 variance and 8% of PC2 variance, and no sample was misclustered (none). The SARS-CoV-2 group shows notably more within-group spread than the Mock group (45.8 vs 20.4), consistent with variable infection/response magnitude across replicates.

![PC association with design variables](figures/pc_association.png)

*File: `figures/pc_association.png` — Association (R²) between the first 6 principal components and design variables: condition.*

![Sample-sample correlation heatmap](figures/sample_correlation.png)

*File: `figures/sample_correlation.png` — Pearson correlation of log2(counts + 1) between all 6 samples (range 0.95–0.99 between different samples).*

No samples were excluded as outliers; all six were retained for differential expression.

## 4. Differential Expression

Low-count genes were filtered prior to DE (48.6% of genes removed, leaving 15,308 of 29,755). DESeq2 was run with design ~condition, contrast SARS-CoV-2 72 hpi vs Mock 72 hpi.

| Measure | Value |
|---|---|
| Contrast | SARS-CoV-2 72 hpi vs Mock 72 hpi |
| Design | `~condition` |
| Genes tested | 15,308 |
| Significant (padj < 0.05) | 3,212 (21.0%) |
| Up-regulated | 1,688 (11.0%) |
| Down-regulated | 1,524 (10.0%) |

![Volcano plot of differential expression](figures/volcano.png)

*File: `figures/volcano.png` — Volcano plot: 15308 genes tested; 3212 with padj < 0.05 (1688 up, 1524 down; design ~condition).*

![MA plot of differential expression](figures/ma_plot.png)

*File: `figures/ma_plot.png` — MA plot: 15308 genes tested; 3212 with padj < 0.05 (1688 up, 1524 down; design ~condition).*

The top upregulated genes include several SARS-CoV-2 viral transcripts, e.g. N (log2FC 11.64, padj < 1e-300), ORF6 (log2FC 11.97, padj < 1e-300) and S (log2FC 11.67, padj < 1e-300), reflecting direct detection of viral RNA in infected samples rather than host transcriptional change. Among host genes, CTGF (log2FC 3.43, padj 2.78e-72) and TM4SF1 (log2FC 4.01, padj 6.17e-112) are strongly upregulated. The top downregulated gene is UGT2A1 (log2FC -2.88, padj 1.05e-36).

| Gene | log2FC | padj | baseMean |
|---|---|---|---|
| ORF6 | 11.97 | < 1e-300 | 17457 |
| ORF3a | 11.92 | < 1e-300 | 41857 |
| ORF10 | 11.79 | < 1e-300 | 25478 |
| ORF8 | 11.68 | < 1e-300 | 46487 |
| S | 11.67 | < 1e-300 | 87812 |
| N | 11.64 | < 1e-300 | 214791 |
| ORF1ab/ORF1a | 11.50 | < 1e-300 | 33421 |
| ORF7a | 11.48 | < 1e-300 | 45383 |
| M | 11.43 | < 1e-300 | 62692 |
| E | 11.98 | 4.01e-300 | 16198 |

| Gene | log2FC | padj | baseMean |
|---|---|---|---|
| UGT2A1 | -2.88 | 1.05e-36 | 222 |
| NWD1 | -1.88 | 1.22e-31 | 656 |
| KDR | -1.56 | 1.17e-25 | 597 |
| SNTG1 | -1.64 | 2.79e-24 | 794 |
| CACNA2D3 | -1.81 | 1.21e-23 | 251 |
| WFIKKN2 | -2.15 | 3.40e-23 | 410 |
| SOSTDC1 | -1.54 | 5.30e-23 | 5260 |
| SLC17A8 | -1.61 | 3.98e-22 | 12023 |
| OPCML | -1.45 | 1.86e-21 | 636 |
| GPC6 | -1.41 | 8.71e-21 | 850 |

![Heatmap of top differentially expressed genes](figures/de_heatmap.png)

*File: `figures/de_heatmap.png` — Top 25 up-regulated (above the line) and 25 down-regulated genes with padj < 0.05, ranked by padj then |log2FC|; z-scored log2(counts + 1) across 6 samples.*

## 5. Gene Set Enrichment

| Library | Term | Overlap | padj |
|---|---|---|---|
| GO Biological Process 2023 | Regulation Of Cell Migration (GO:0030334) | 39/434 | 1.23e-08 |
| GO Biological Process 2023 | Positive Regulation Of Cell Migration (GO:0030335) | 29/272 | 6.19e-08 |
| GO Biological Process 2023 | Positive Regulation Of Cell Motility (GO:2000147) | 25/221 | 2.56e-07 |
| GO Biological Process 2023 | Regulation Of Wound Healing (GO:0061041) | 12/47 | 7.04e-07 |
| GO Biological Process 2023 | Regulation Of Cell Population Proliferation (GO:0042127) | 47/766 | 6.64e-06 |
| GO Biological Process 2023 | Extracellular Matrix Organization (GO:0030198) | 20/176 | 6.64e-06 |
| GO Biological Process 2023 | Regulation Of Blood Coagulation (GO:0030193) | 9/29 | 6.64e-06 |
| GO Biological Process 2023 | Integrin-Mediated Signaling Pathway (GO:0007229) | 14/85 | 6.64e-06 |
| GO Biological Process 2023 | Regulation Of Apoptotic Process (GO:0042981) | 42/705 | 5.05e-05 |
| GO Biological Process 2023 | Positive Regulation Of Peptidyl-Tyrosine Phosphorylation (GO:0050731) | 15/130 | 2.09e-04 |
| KEGG 2021 Human | Focal adhesion | 28/201 | 3.77e-11 |
| KEGG 2021 Human | ECM-receptor interaction | 15/88 | 5.01e-07 |
| KEGG 2021 Human | Complement and coagulation cascades | 14/85 | 1.73e-06 |
| KEGG 2021 Human | Regulation of actin cytoskeleton | 21/218 | 8.14e-06 |
| KEGG 2021 Human | AGE-RAGE signaling pathway in diabetic complications | 14/100 | 8.57e-06 |
| KEGG 2021 Human | Small cell lung cancer | 12/92 | 1.10e-04 |
| KEGG 2021 Human | PI3K-Akt signaling pathway | 25/354 | 1.10e-04 |
| KEGG 2021 Human | Proteoglycans in cancer | 18/205 | 1.19e-04 |
| KEGG 2021 Human | Tight junction | 15/169 | 5.82e-04 |
| KEGG 2021 Human | Amoebiasis | 11/102 | 0.001 |

![GO Biological Process enrichment, upregulated genes](figures/enrichment_upregulated_go_biological_process_2023.png)

*File: `figures/enrichment_upregulated_go_biological_process_2023.png` — Enrichr GO_Biological_Process_2023, upregulated genes (500 input genes): top 10 terms by adjusted p-value.*

![KEGG enrichment, upregulated genes](figures/enrichment_upregulated_kegg_2021_human.png)

*File: `figures/enrichment_upregulated_kegg_2021_human.png` — Enrichr KEGG_2021_Human, upregulated genes (500 input genes): top 10 terms by adjusted p-value.*

Upregulated genes are enriched for Regulation Of Cell Migration (GO:0030334) (39/434 genes, padj 1.23e-08), Extracellular Matrix Organization (GO:0030198) (20/176 genes, padj 6.64e-06), Focal adhesion (28/201 genes, padj 3.77e-11) and Complement and coagulation cascades (14/85 genes, padj 1.73e-06), together with Inflammatory Response (GO:0006954) (15/236 genes, padj 0.023) and Apoptosis (12/142 genes, padj 0.003).

| Library | Term | Overlap | padj |
|---|---|---|---|
| GO Biological Process 2023 | Metal Ion Transport (GO:0030001) | 15/171 | 0.002 |
| GO Biological Process 2023 | Inorganic Cation Import Across Plasma Membrane (GO:0098659) | 9/103 | 0.088 |
| GO Biological Process 2023 | Calcium Ion Transmembrane Import Into Cytosol (GO:0097553) | 8/83 | 0.088 |
| GO Biological Process 2023 | Monoatomic Ion Transport (GO:0006811) | 9/107 | 0.088 |
| GO Biological Process 2023 | Sodium Ion Transport (GO:0006814) | 8/89 | 0.088 |
| GO Biological Process 2023 | Monoatomic Cation Transmembrane Transport (GO:0098655) | 15/281 | 0.088 |
| GO Biological Process 2023 | Calcium-Independent Cell-Cell Adhesion Via Plasma Membrane Cell-Adhesion Molecules (GO:0016338) | 4/19 | 0.088 |
| GO Biological Process 2023 | Inorganic Cation Transmembrane Transport (GO:0098662) | 15/284 | 0.088 |
| GO Biological Process 2023 | Potassium Ion Transport (GO:0006813) | 9/122 | 0.094 |
| GO Biological Process 2023 | Calcium Ion Import Across Plasma Membrane (GO:0098703) | 5/37 | 0.094 |
| KEGG 2021 Human | Bile secretion | 9/90 | 0.013 |
| KEGG 2021 Human | Pancreatic secretion | 9/102 | 0.017 |
| KEGG 2021 Human | Proximal tubule bicarbonate reclamation | 4/23 | 0.056 |
| KEGG 2021 Human | Mineral absorption | 6/60 | 0.056 |
| KEGG 2021 Human | Aldosterone synthesis and secretion | 7/98 | 0.125 |
| KEGG 2021 Human | Adrenergic signaling in cardiomyocytes | 8/150 | 0.289 |
| KEGG 2021 Human | Oxytocin signaling pathway | 8/154 | 0.289 |
| KEGG 2021 Human | Adipocytokine signaling pathway | 5/69 | 0.289 |
| KEGG 2021 Human | PPAR signaling pathway | 5/74 | 0.313 |
| KEGG 2021 Human | Gastric acid secretion | 5/76 | 0.313 |

![GO Biological Process enrichment, downregulated genes](figures/enrichment_downregulated_go_biological_process_2023.png)

*File: `figures/enrichment_downregulated_go_biological_process_2023.png` — Enrichr GO_Biological_Process_2023, downregulated genes (384 input genes): top 10 terms by adjusted p-value.*

![KEGG enrichment, downregulated genes](figures/enrichment_downregulated_kegg_2021_human.png)

*File: `figures/enrichment_downregulated_kegg_2021_human.png` — Enrichr KEGG_2021_Human, downregulated genes (384 input genes): top 10 terms by adjusted p-value.*

Downregulated genes are enriched for Metal Ion Transport (GO:0030001) (15/171 genes, padj 0.002) and Bile secretion (9/90 genes, padj 0.013)/Pancreatic secretion (9/102 genes, padj 0.017) KEGG terms, which largely reflect ion-transport and secretory-epithelium functions rather than classic immune pathways.

## 6. Comparison with the Published Study

(Jacob et al., 2020; PMID: 33010822) reports that SARS-CoV-2 infection of hiPSC-derived choroid plexus organoids caused increased cell death and "transcriptional dysregulation indicative of an inflammatory response and cellular function deficits." This analysis is consistent with that description at a pathway level: the upregulated set shows significant enrichment of Inflammatory Response (GO:0006954) (15/236 genes, padj 0.023) and Apoptosis (12/142 genes, padj 0.003), and the downregulated set shows loss of ion-transport/secretory terms such as Metal Ion Transport (GO:0030001) (15/171 genes, padj 0.002), consistent with a "cellular function deficit" in a secretory epithelium. The abstract does not name specific genes, so no gene-level comparison table is included; the agreement here is at the level of the general pathway themes it describes.

## 7. Biological Interpretation

- Viral transcripts (N (log2FC 11.64, padj < 1e-300), ORF6 (log2FC 11.97, padj < 1e-300), S (log2FC 11.67, padj < 1e-300), among others) dominate the largest log2 fold-changes, confirming that the SARS-CoV-2 samples were productively infected, as expected from the study design.
- The upregulated host response is dominated by Focal adhesion (28/201 genes, padj 3.77e-11), Extracellular Matrix Organization (GO:0030198) (20/176 genes, padj 6.64e-06) and Regulation Of Cell Migration (GO:0030334) (39/434 genes, padj 1.23e-08) terms, alongside a smaller but significant Inflammatory Response (GO:0006954) (15/236 genes, padj 0.023) signature, indicating that ECM remodeling/adhesion changes are the largest-magnitude transcriptional theme, with inflammation present but less dominant by gene count.
- Apoptosis (12/142 genes, padj 0.003) enrichment among upregulated genes is consistent with the increased cell death reported in the source study, though this analysis only shows transcriptional correlates, not direct cell-death measurements.
- Downregulated genes are enriched for Metal Ion Transport (GO:0030001) (15/171 genes, padj 0.002), pointing to reduced expression of ion-channel/transporter genes that are not part of the inflammatory narrative and represent a result outside the paper's stated focus.
- ITG* genes: 9 of 21 significant (8 up, 1 down), several of which appear among the focal-adhesion/ECM upregulated genes, illustrate a broad integrin-pathway shift accompanying infection.

## 8. Limitations

- Each condition has only 3 (Mock) and 3 (SARS-CoV-2) replicates, limiting statistical power and the ability to detect subtler effects.
- The SARS-CoV-2 group shows higher PCA spread than the Mock group (45.8 vs 20.4), indicating more variable infection response across replicates that should be considered when interpreting individual gene calls.
- Enrichment analysis used the Enrichr default background rather than a study-specific background, and the lowest-magnitude viral-gene fold changes are not directly comparable to host gene changes because viral transcripts are absent in Mock samples by design.
- This is a re-analysis of the authors' processed count matrix (GSE157852_CPO_RawCounts.txt.gz) rather than a re-alignment from FASTQ, so any upstream processing choices made by the original authors are inherited as-is.
- No paired/batch covariate was available beyond condition; the design used is ~condition.

## Methods

### Data

Counts: authors' supplementary file `GSE157852_CPO_RawCounts.txt.gz` from GEO GSE157852 (MD5 `988011dda6a9b9c5c2830a862aaeba94`); values: raw integer counts; gene IDs: symbol; 0 duplicate gene IDs summed.

### Quality control

Library size is the total count per sample and genes detected the number of genes with at least one count, both computed before low-count filtering. PCA: singular value decomposition of gene-centred log2(counts + 1) of the filtered matrix (15,308 genes), without scaling. Sample correlation: Pearson correlation of log2(counts + 1) of the filtered matrix (15,308 genes).

### Low-count filtering

Genes were kept when at least 2 samples had ≥ 10 counts: 15,308 of 29,755 genes kept (48.6% removed).

### Differential expression

Differential expression was tested with PyDESeq2 0.5.4 using the design `~condition`, comparing SARS-CoV-2 72 hpi with Mock 72 hpi (reference level). PyDESeq2 fits a negative binomial generalised linear model per gene, with median-of-ratios size factors and dispersions shrunk towards a fitted trend, and tests the condition coefficient with a Wald test. P-values were adjusted with the Benjamini–Hochberg method after PyDESeq2's default independent filtering and Cook's-distance outlier handling; 15,308 genes received an adjusted p-value. Genes with padj < 0.05 are called significant. Log2 fold changes are unshrunken maximum-likelihood estimates.

### Enrichment (up-regulated)

Over-representation analysis of up-regulated genes used Enrichr through GSEApy 1.3.1 against GO_Biological_Process_2023, KEGG_2021_Human (human). Genes with padj < 0.05 and |log2FC| ≥ 1.0 (576 genes) were ranked by padj, then |log2FC|, and the top 500 submitted as 500 gene symbols. Enrichr's default background (all genes in each library) was used, not the genes tested here. Terms are reported by Enrichr's adjusted p-value (Benjamini–Hochberg). Exact input: `analysis/enrichment_input_upregulated.txt`.

### Enrichment (down-regulated)

Over-representation analysis of down-regulated genes used Enrichr through GSEApy 1.3.1 against GO_Biological_Process_2023, KEGG_2021_Human (human). Genes with padj < 0.05 and |log2FC| ≥ 1.0 (384 genes) were ranked by padj, then |log2FC|, and the top 384 submitted as 384 gene symbols. Enrichr's default background (all genes in each library) was used, not the genes tested here. Terms are reported by Enrichr's adjusted p-value (Benjamini–Hochberg). Exact input: `analysis/enrichment_input_downregulated.txt`.

### Software

| Software | Version |
|---|---|
| Python | 3.11.14 |
| NumPy | 2.4.6 |
| pandas | 2.3.3 |
| PyDESeq2 | 0.5.4 |
| Matplotlib | 3.11.2 |
| SciPy | 1.17.1 |
| GSEApy | 1.3.1 |

### Reproducibility

Every analysis step is logged in `analysis/tool_calls.jsonl`, and `analysis/replay.py` re-runs them without the LLM. This Methods section is generated from those steps, not written by the LLM.

## References

(Jacob et al., 2020; PMID: 33010822)

---

**Disclaimer:** This report was generated by an AI system. Large language models can produce inaccurate statements (hallucinations). All biological claims, gene annotations, pathway interpretations, and cited references should be independently verified before use in publications or clinical decisions. Biological roles listed alongside gene names are AI-generated summaries and must be verified against primary databases (UniProt, NCBI Gene).

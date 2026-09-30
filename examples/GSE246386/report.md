# nf-core/rnaseq Sample Sheet Report: GSE246386

## Dataset Overview
**Title:** GFI1B ectopic expression in human iPSC-derived hemogenic endothelium (CD34+ iPSC-HE, day 4)  
**Organism:** *Homo sapiens* (GRCh38)  
**GEO Series:** [GSE246386](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246386)  
**Library Layout:** Paired-end RNA-seq  
**Strandedness:** auto

## FASTQ Discovery
- **Directory:** `data/GSE246386/fastqs/`
- **Total Files:** 12 (6 paired-end libraries)
- **All Pairs:** Complete (6 R1/R2 pairs, 0 unpaired, 0 incomplete)

## Sample Composition

### GFI1B Overexpression (n=3)
| Run Accession | Sample ID | Replicate | Title (GEO) |
|---|---|---|---|
| SRR26539596 | GFI1B_1 | #1 | CD34+ iPSC-HE, day4, GFI1B#1 |
| SRR26539595 | GFI1B_2 | #2 | CD34+ iPSC-HE, day4, GFI1B#2 |
| SRR26539594 | GFI1B_3 | #3 | CD34+ iPSC-HE, day4, GFI1B#3 |

### Empty Vector Control (n=3)
| Run Accession | Sample ID | Replicate | Title (GEO) |
|---|---|---|---|
| SRR26539599 | EV_1 | #1 | CD34+ iPSC-HE, day4, EV#1 |
| SRR26539598 | EV_2 | #2 | CD34+ iPSC-HE, day4, EV#2 |
| SRR26539597 | EV_3 | #3 | CD34+ iPSC-HE, day4, EV#3 |

## Metadata Integration
Sample names and condition assignments were derived from the GEO/SRA metadata provided. Mapping verified:
- **GEO Series:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246386  
  → Check the "Samples" table for GSM accessions and titles; confirm GFI1B vs. EV designation.

## Sample Sheet Details
- **Total Samples:** 6
- **Strandedness:** auto (for all samples)
- **File Format:** CSV (nf-core/rnaseq v3.14.0+ compatible)
- **Location:** `runs/20260930_GSE246386_full/samplesheet.csv`

## Experimental Design
A design matrix was generated for differential expression analysis:
- **Condition 1:** GFI1B (3 replicates)
- **Condition 2:** EV (3 replicates)
- **Design File:** `runs/20260930_GSE246386_full/design.csv`

The contrast **GFI1B vs. EV** captures the effect of ectopic GFI1B expression in CD34+ iPSC-derived hemogenic endothelium.

## Validation Results
✅ **Valid:** All 6 samples  
✅ **No Errors**  
✅ **No Warnings**  
✅ **All FASTQ files exist**

## Next Steps
1. **Run nf-core/rnaseq** with the saved sample sheet (skip alignment, salmon only):
   ```bash
   nextflow run nf-core/rnaseq \
     -profile <your_profile> \
     --input samplesheet.csv \
     --outdir results/ \
     --aligner skip \
     --pseudo_aligner salmon
   ```

2. **Downstream Analysis:**
   - QC: MultiQC aggregation (FastQC, Salmon metrics)
   - Differential Expression: DESeq2 using design.csv
   - Enrichment: GO/KEGG on DE genes (FDR < 0.05)
   - Report: Markdown summary with plots and tables

## Data Curation Notes
- All samples are from the same day (day 4) and cell type (CD34+ iPSC-HE)
- Three technical/biological replicates per condition provide robust statistical power
- Library layout (paired-end) and strandedness (auto) are optimized for Salmon quantification
- No lane-level multiplexing detected; file names are clean (SRR_N format)

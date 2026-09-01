#!/usr/bin/env bash
# Download subsampled FASTQs for GSE225648 from SRA.
# Dataset: Drosophila melanogaster RNA-seq, 6 paired-end samples
# Study: "Reduced levels of ALS gene DCTN1 induce motor defects in Drosophila"
# Platform: Illumina HiSeq 2500
# Design: 2 conditions (ME18 knockdown vs CTRL), 3 reps each
#
# Requires: sra-tools (fastq-dump). Install via:
#   brew install sratoolkit
#
# Downloads first 50k read pairs per sample (~5MB per file).

set -euo pipefail

OUTDIR="$(dirname "$0")/fastqs"
MAX_SPOTS=50000

mkdir -p "$OUTDIR"

# GSM -> SRR mapping:
#   GSM7053673 (Act_ME18_Rep1) -> SRR23558740
#   GSM7053674 (Act_ME18_Rep2) -> SRR23558739
#   GSM7053675 (Act_ME18_Rep3) -> SRR23558738
#   GSM7053676 (Act_CTRL_Rep2) -> SRR23558737
#   GSM7053677 (Act_CTRL_Rep3) -> SRR23558736
#   GSM7053678 (Act_CTRL_Rep4) -> SRR23558735

ACCESSIONS=(
    SRR23558740
    SRR23558739
    SRR23558738
    SRR23558737
    SRR23558736
    SRR23558735
)

for acc in "${ACCESSIONS[@]}"; do
    if [ -f "${OUTDIR}/${acc}_1.fastq.gz" ] && [ -f "${OUTDIR}/${acc}_2.fastq.gz" ]; then
        echo "Skipping ${acc} (already downloaded)"
        continue
    fi
    echo "Downloading ${acc} (first ${MAX_SPOTS} spots)..."
    fastq-dump "$acc" \
        --outdir "$OUTDIR" \
        --split-files \
        --gzip \
        -X "$MAX_SPOTS"
    echo "  -> ${acc}_1.fastq.gz, ${acc}_2.fastq.gz"
done

echo ""
echo "Done. Files in ${OUTDIR}:"
ls -lh "$OUTDIR"

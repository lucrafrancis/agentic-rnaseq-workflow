#!/usr/bin/env bash
# Download nf-core/rnaseq test data (yeast, subsampled to 50k reads).
# Source: https://github.com/nf-core/test-datasets/tree/rnaseq
# Dataset: GSE110004 — 7 runs, 3 conditions (WT, RAP1_UNINDUCED, RAP1_IAA_30M)
# Mix of paired-end and single-end samples.

set -euo pipefail

BASE_URL="https://raw.githubusercontent.com/nf-core/test-datasets/rnaseq/testdata/GSE110004"
OUTDIR="$(dirname "$0")/fastqs"

mkdir -p "$OUTDIR"

for acc in SRR6357070 SRR6357071 SRR6357072 SRR6357073 SRR6357074 SRR6357075 SRR6357076; do
    echo "Downloading ${acc}..."
    curl -sL -o "${OUTDIR}/${acc}_1.fastq.gz" "${BASE_URL}/${acc}_1.fastq.gz"
    # R2 only exists for paired-end samples; skip 404s silently
    curl -sL -f -o "${OUTDIR}/${acc}_2.fastq.gz" "${BASE_URL}/${acc}_2.fastq.gz" 2>/dev/null || true
done

echo "Done. Files in ${OUTDIR}:"
ls -lh "$OUTDIR"

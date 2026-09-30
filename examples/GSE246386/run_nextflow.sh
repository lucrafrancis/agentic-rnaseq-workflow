#!/usr/bin/env bash
set -euo pipefail

nextflow run nf-core/rnaseq \
  -r 3.26.0 \
  --input runs/20260930_GSE246386_full/samplesheet.csv \
  --outdir runs/20260930_GSE246386_full/results \
  -profile docker \
  -c runs/20260930_GSE246386_full/custom.config \
  -params-file runs/20260930_GSE246386_full/nf_params.yml

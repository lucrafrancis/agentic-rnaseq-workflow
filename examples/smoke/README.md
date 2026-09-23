# Smoke tests

Cheap end-to-end runs that exercise the real LLM against the new code paths before
spending compute on nextflow. Each test lists the command, what you should see at the
approval step, and what to check afterwards. Anything that doesn't match is a finding.

Prompts deliberately state only the dataset and the comparison — file choice, sample
mapping and design are left to the agents, so these test their judgement, not just the
plumbing.

After every run, note `usage.jsonl` (tokens per agent) in the run directory.

---

## 1. GSE157852 — counts mode, single factor

```bash
uv run python run.py --analyze examples/smoke/GSE157852_counts/prompt.txt
```

Run directory: `runs/<date>_GSE157852_counts/`. No AWS or nextflow needed.

**Tests:** author file chosen over NCBI; column names that don't match GEO titles are
mapped correctly; a specific contrast chosen from a three-group design.

**At the approval step, expect:**

- Source: `GSE157852_CPO_RawCounts.txt.gz` (author), values `raw_integer_counts`,
  gene IDs `symbol`, ~29.7k genes.
- Mapping (the column suffix `S1`–`S9` follows GSM order):

  | file column | GSM | condition |
  |---|---|---|
  | CPO_Mock_72hpi_S1–S3 | GSM4776784–86 | Mock 72 hpi |
  | CPO_SARS-CoV-2_24hpi_S4–S6 | GSM4776787–89 | SARS-CoV-2 24 hpi |
  | CPO_SARS-CoV-2_72hpi_S7–S9 | GSM4776790–92 | SARS-CoV-2 72 hpi |

- Condition from the `treatment` characteristic. No covariates.
- Either keeping all 9 samples or dropping the 24 hpi group is acceptable (keeping them
  gives DESeq2 more samples for dispersion estimates) — but the agent should say which.

**Afterwards, check:**

- Contrast `["condition", "SARS-CoV-2 72 hpi", "Mock 72 hpi"]`, design `~condition`.
- Report states 3 vs 3 replicates and cites GSE157852.
- `analysis/replay.py` exists; optionally run it and compare `de_results.csv`.

---

## 2. GSE164073 — counts mode, hidden covariate

```bash
uv run python run.py --analyze examples/smoke/GSE164073_counts/prompt.txt
```

Run directory: `runs/<date>_GSE164073_counts/`.

**Tests:** mapping of prefixed column names; noticing — unprompted — that tissue must be
adjusted for. The prompt intentionally does not mention tissue.

**At the approval step, expect:**

- Source: `GSE164073_Eye_count_matrix.csv.gz` (author), values `raw_integer_counts`,
  gene IDs `symbol`.
- Mapping — `MW1`–`MW18` follow GSM order (GSM4996084–GSM4996101):

  | file columns | GSMs | tissue | infection |
  |---|---|---|---|
  | MW1–3_cornea_mock | GSM4996084–86 | cornea | mock |
  | MW4–6_cornea_CoV2 | GSM4996087–89 | cornea | SARS-CoV-2, MOI = 1.0 |
  | MW7–9_limbus_mock | GSM4996090–92 | limbus | mock |
  | MW10–12_limbus_CoV2 | GSM4996093–95 | limbus | SARS-CoV-2, MOI = 1.0 |
  | MW13–15_sclera_mock | GSM4996096–98 | sclera | mock |
  | MW16–18_sclera_CoV2 | GSM4996099–101 | sclera | SARS-CoV-2, MOI = 1.0 |

- Condition from `infection`; `tissue` listed as a covariate (`time_point` is constant,
  so it is dropped).

**Afterwards, check:**

- `run_deseq2` called with `covariates=["tissue"]` — design `~tissue + condition`.
  **Main failure to watch for:** design `~condition`, which pools three tissues without
  adjusting for them.
- Report mentions the tissue adjustment in Methods and 9 vs 9 samples.
- PCA likely separates samples by tissue more than by infection — the report should
  notice this.

---

## 3. GSE246386 — full pipeline (download → samplesheet → nextflow → analysis)

```bash
uv run python run.py examples/smoke/GSE246386_full/prompt.txt
```

Run directory: `runs/<date>_GSE246386_full/`. Needs Nextflow + Docker; runs salmon-only
locally. FASTQs go to `./data/GSE246386/fastqs/` (gitignored).

**Before running:** free disk space. Expect ~6 GB of FASTQs plus reference download,
Docker images and a `work/` directory that can grow to several times the input size
(earlier runs left 95 GB). Delete `work/` once you've checked the results.

**Tests:** resolving a recent series via its BioProject (no SRA link in GEO's search
record); runs whose sample labels ENA leaves blank (3 of 6 here, filled from GEO);
paired-end handling end to end; the full stage chain and `run_state.json`.

**At the download approval, expect:** 6 runs, 12 files, ~5.9 GB, every run labelled:

| run | GSM | title |
|---|---|---|
| SRR26539599 | GSM7868165 | CD34+ iPSC-HE, day4, EV#1 |
| SRR26539598 | GSM7868166 | CD34+ iPSC-HE, day4, EV#2 |
| SRR26539597 | GSM7868167 | CD34+ iPSC-HE, day4, EV#3 |
| SRR26539596 | GSM7868168 | CD34+ iPSC-HE, day4, GFI1B#1 |
| SRR26539595 | GSM7868169 | CD34+ iPSC-HE, day4, GFI1B#2 |
| SRR26539594 | GSM7868170 | CD34+ iPSC-HE, day4, GFI1B#3 |

**Worth testing resume here:** interrupt (Ctrl-C) during the download, then
`uv run python run.py --resume runs/<date>_GSE246386_full` — it should regenerate the
script without the LLM and skip files that already pass MD5.

**At the samplesheet approval, expect:** 6 rows, `fastq_1`/`fastq_2` paired correctly
(`_1`/`_2`), strandedness `auto`, sample names reflecting EV/GFI1B and replicate.
`design.csv`: EV × 3, GFI1B × 3.

**At the nextflow approval, expect:** `genome: GRCh38`, `skip_alignment: true`,
`pseudo_aligner: salmon`, resource limits matching this machine in `custom.config`.

**Afterwards, check:**

- Contrast GFI1B vs EV (EV as reference), design `~condition`, 3 vs 3.
- Pipeline warnings summary makes sense; `usage.jsonl` shows tokens per agent.
- Optional cross-check: the authors' raw counts (`GSE246386_raw_count.txt.gz`, Ensembl
  IDs) should correlate strongly with the Salmon gene counts — a quick sanity check
  that the preprocessing is sound.

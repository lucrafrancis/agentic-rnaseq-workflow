# agentic-rnaseq-workflow

Raw Anthropic API agent loop (no frameworks) that builds nf-core/rnaseq sample sheets from FASTQ directories, then submits pipeline runs. Tools return structured summary dicts the LLM reasons over.

## Architecture

Four stages, with human approval **between** stages (not inside the agent loop):

0. **Download agent** (`agents/download/`) — resolves GEO/SRA accessions via NCBI/ENA APIs, generates a download script (aspera > aria2c > curl), validates MD5 checksums. Runs whenever the prompt contains an accession (it checks existing files by MD5 and only downloads what's missing); `--skip-download` bypasses it. Rejecting the download script exits. GSE resolution falls back to the series SOFT record's BioProject/SRA link when GEO's search record has none (most series since ~2023), and labels every run with its GSM and title from GEO via the SRX (ENA often leaves these blank). SuperSeries are detected and their SubSeries listed; NCBI/ENA outages are reported as `lookup_failed`, not "no runs". The agent can download a subset of runs (`runs=`). **Hard stop:** runs without a sample ID and title are never downloaded or passed on (`unlabelled_runs`, enforced in the tool and in `run.py`). The selected runs' labels are written to `sample_metadata.csv` and handed to the samplesheet agent.
1. **Samplesheet agent** (`agents/samplesheet/`) — scans FASTQs, reads metadata, matches pairs, drafts/validates/saves a sample sheet, writes a report. Fully implemented and tested with real GEO data.
2. **Submission agent** (`agents/submission/`) — reads the prompt and configures nextflow params (genome, skip_alignment, etc.) via an open-ended tool. Has a `run_command` tool to check system resources (RAM, CPUs) before configuring — always sets resource limits based on the actual machine; with the docker profile the memory ceiling is Docker's memory, not the host's, and `configure_submission` (and the troubleshooter's `propose_fix`) reject a `max_memory` above the host RAM or Docker's memory. Writes `run_nextflow.sh` + `nf_params.yml` + `custom.config`. Pipeline params go in the YAML; resource limits (`max_memory`, `max_cpus`, `max_time`) go in `custom.config` as nextflow `resourceLimits` to avoid nf-schema validation warnings. Human approves the script before execution. Salmon-only runs reuse a shared reference (`agents/submission/reference.py`): `configure_submission` swaps `--genome` for the cached `gtf`/`transcript_fasta`/`salmon_index` in `data/reference/<genome>/nf-core-rnaseq-<revision>/` when a complete one exists, otherwise sets `save_reference: true`, and `run.py` copies what nf-core published to `results/genome/` into the cache after a successful run (temporary folder + rename; refuses incomplete or ambiguous files; `reference.json` records source run, Salmon version, MD5s and the original GTF/FASTA (`gtf_source`, `fasta_source`) so reusing runs still report the annotation). Delete a cache folder to force a rebuild.
3. **Analysis agent** (`agents/analysis/`) — downstream DE, enrichment, QC, and reporting on nf-core/rnaseq outputs. Python-only (PyDESeq2, gseapy, numpy PCA). Optionally accepts a paper PDF as context via Claude's native base64 content blocks. Fetches GEO metadata and PubMed abstracts for context/citations. Fully implemented and tested.

The full flow in `run.py`: download agent (if needed) → human approves download script → samplesheet agent → human approves sample sheet → submission agent configures params → human approves `run_nextflow.sh` → nextflow runs (with troubleshooting on failure) → post-run warning review → analysis agent runs downstream analysis.

## Running

```bash
uv run python run.py <prompt.txt>                # full pipeline (stages 0-3)
uv run python run.py <prompt.txt> --skip-download  # FASTQs already local
uv run python run.py --resume <run_dir>          # resume from existing run directory
uv run python run.py --analyze <run_dir>         # re-run analysis (stage 3) on an existing run
uv run python run.py --analyze <prompt.txt>      # analysis only: GEO counts if the prompt has an accession
uv run python run.py --analyze <prompt.txt> --skip-download  # analysis only, local count matrix
```

Completed stages are recorded in `run_state.json` (only after approval/validation), along with the source prompt path. `--resume` trusts that file, not artifact existence: download done → skip; interrupted download (metadata has `output_dir`) → script regenerated without the LLM, skipping files with valid MD5; samplesheet approved → skip, but an unapproved `samplesheet.csv` is re-presented for approval without re-running the agent; `params.json` → skip submission agent (the script is still re-approved); `salmon.merged.gene_counts.tsv` in `results/` → skip nextflow.

`--analyze <run_dir>` attaches to an existing run directory (via `resume_run`), reads the original `prompt.txt` as context, and runs the analysis agent (with or without `results/`). Output goes into the same run dir (`analysis/`), not a new one. `--analyze <prompt.txt>` creates a new analysis-only run (`mode: "analysis"` in `run_state.json`, so `--resume` never enters stages 1-2). If the prompt has a GEO accession, the download agent runs in **counts mode** (`agents/download/counts.py`): `list_geo_count_sources` → `preview_geo_file` → `fetch_geo_counts` → `save_geo_design`. It prefers author supplementary files with raw counts, falling back to NCBI-generated counts; the LLM picks the file, gene ID column and column→GSM mapping, the tools parse and validate (one column per GSM, numeric, non-negative, duplicate IDs summed, samples renamed to GEO titles, symbols from NCBI `gene_info` cached in `data/reference/`). The counts agent stops after the summary — analysis is a separate stage. Writes `counts.tsv`, `counts_metadata.json` (provenance incl. MD5) and `design.csv` (condition, varying GEO characteristics as covariates, and constant ones such as time point kept as facts); the human approves the mapping table before analysis. Raw GEO files are kept in `geo/`. Per-sample files (`_RAW.tar`) and xlsx are listed but not yet parsed. With `--skip-download` (or no accession) the prompt's local count matrix is used. Tool logs (`tool_calls.jsonl`, written by `core/loop.py`) record every path inside the repo relative to it, so committed logs don't expose the local home directory. Each analysis writes `analysis/replay.py` (repo paths relative to the repo root, which the script finds by walking up to `pyproject.toml`, so a copy under `examples/` runs from any clone), which re-runs the successful tool calls without the LLM (`python -m agents.analysis.replay <run_dir>` for older runs). An optional interactive prompt lets the user add extra instructions (saved as `analysis_prompt.txt` if provided); Enter continues without them, "skip" aborts.

The prompt file describes the data (FASTQ location, organism, metadata path, strandedness). See `examples/smoke/` for working prompts. Run directories are named `runs/<datestamp>_<project>/` (derived from prompt file's parent directory).

## Key design rules

- No frameworks (no LangChain, no CrewAI)
- Tools return `dict[str, Any]` summaries — the LLM never sees raw data
- Human approval belongs between stages, never inside the agent loop
- Keep code minimal — no bloat, no premature abstractions
- Tests are offline (no API calls)
- Tools must enforce correctness, not the LLM — paths are resolved to absolute in the tools themselves (`scan_fastqs`, `draft_samplesheet`), the samplesheet agent never relays paths or CSV (`scan_fastqs`/`match_pairs` keep files and pairs in the session, `draft_samplesheet` takes only `pair_id` + sample name + strandedness, `validate_samplesheet`/`save_samplesheet` take no arguments; validation checks every FASTQ exists, save refuses an unvalidated draft and always writes to the canonical run directory location, and `write_report` prepends a ⚠ banner if no sheet was saved), the sample sheet approval screen shows a code-built sample → condition → run → GSM → GEO title table, `SubmissionParams` routes resource limits to `custom.config` vs pipeline params to `nf_params.yml` automatically, `configure_submission` forces `pseudo_aligner: salmon` when `skip_alignment` is set (so counts are always produced), and `generate_figures`/`write_report` supply real figure paths and captions and append the disclaimer rather than trusting the LLM to get filenames or caveats right. Tools may make lossless mechanical corrections to LLM output, and reject anything else with a reason — never silently delete what the LLM wrote. The LLM is a lossy intermediary; don't trust it to preserve values faithfully between tool calls.

## Prerequisites

- Python 3.11+, uv
- `ANTHROPIC_API_KEY` environment variable
- Nextflow + Docker (for pipeline submission — not needed for samplesheet agent)

## Config

`core/config.py` — model name, max tokens/turns, pipeline version, paths. Single source of truth. `MODEL` (Haiku 4.5) for the download/counts, samplesheet and submission agents, `ANALYSIS_MODEL` (Sonnet 5, as in agentic-scrna-workflow; `ANALYSIS_MODEL=haiku` to override) for the analysis agent, `MODEL_SONNET` (Sonnet 4.5) for troubleshooting and warning review.

Agent loops and the troubleshooter use automatic prompt caching (top-level `cache_control`; Haiku 4.5 only caches prefixes ≥ 4096 tokens). Token usage and estimated cost per agent run are printed and appended to the run's `usage.jsonl` (prices in `config.PRICE_PER_MTOK`; a model without a price is logged with `estimated_cost_usd: null`).

## Analysis agent tools

`fetch_geo_metadata` → `fetch_abstract` → `scan_results` → `load_counts` → `inspect_counts` → `set_design` (if needed) → `compute_qc` → `read_multiqc` (nf-core runs) → `filter_low_counts` → `run_deseq2` → `get_top_genes` → `run_enrichment` → `query_genes` / `search_enrichment` → `summarize_findings` → `generate_figures` → `write_report`

- `fetch_geo_metadata` / `fetch_abstract` — NCBI E-utilities for GEO series metadata and PubMed abstracts. Provides cell type, organism, experimental context, and citable references. Called first when a GEO accession is in the prompt.
- `read_multiqc` — `scan_results` searches `results/` for every `multiqc_general_stats.txt` (the layout differs between nf-core versions and aligners, so nothing is hardcoded); the tool returns whatever columns that version reports for the analysis samples (the LLM picks the file if there are several) and makes them citable as `multiqc.<column>.min|median|max`.
- `inspect_counts` — detects whether data is raw counts vs normalised (TPM/FPKM/log). Prevents running DESeq2 on pre-normalised data. Also flags Salmon fractional counts (`salmon_fractional: true`) so downstream tools know rounding is safe.
- `run_deseq2` — takes optional `covariates` (design columns such as donor/batch/tissue; categorical): design = `~ covariates + factor`. Refuses confounded designs, unknown levels, constant covariates and non-identifier column names; returns the formula (also in `summarize_findings`). Rounds Salmon fractional counts to integers automatically (only when data looks like raw counts with minor fractional noise, not normalised). Reports `counts_rounded` and `pct_non_integer_before_rounding` in its return dict. Also returns `low_replication_warning` when any group has < 3 replicates — the agent must include this caveat in the report.
- `save_design` (samplesheet agent) writes `design.csv` for traceability; `set_design` (analysis agent) is the fallback when none exists
- `run_enrichment(direction, padj_max=0.05, lfc_min=1.0, max_genes=500)` selects genes itself from the DESeq2 results — the LLM never passes gene lists: padj and |log2FC| filters, split by direction, ranked by padj → |log2FC| → gene ID (deterministic; the same `_ranked` helper feeds the top-genes tables, heatmap and DESeq2 preview), capped, converted to symbols. Saves `analysis/enrichment_input_<label>.txt` (exact genes) and `analysis/enrichment_<label>.csv` (raw Enrichr results). Replay scripts therefore contain selection parameters, not hardcoded genes.
- `query_genes(symbols, prefix)` and `search_enrichment(query, direction)` — the LLM isn't limited to the top genes/terms: it can look up any gene or gene family (e.g. prefix `ITG`) and search every Enrichr term (rank, padj, genes behind it) before making a claim. Citable as `{{genes:PREFIX*}}` (family counts) and `{{term:up|down:<term>}}`. Both are read-only (skipped in replay).
- **Paper comparison:** gene symbols named in fetched abstracts (exact matches to genes in the data; all-caps or containing a digit) get this analysis's result in `summarize_findings["paper_genes"]` and the required `{{table:paper_genes}}`, so agreement with the paper is shown, not asserted.
- **Biological Interpretation is short and anchored:** `write_report` rejects an interpretation/discussion section over 250 words or with any paragraph/bullet containing no placeholder.
- `summarize_findings` returns `software_versions` (real installed versions), `pipeline_versions` (parsed from nf-core's `nf_core_rnaseq_software_mqc_versions.yml` — exact Nextflow, pipeline, and Salmon versions), and `data_source` ("nf-core" or "user-provided") to prevent hallucination in the Methods section
- `generate_figures` draws every figure before the report is written and returns each path with a factual caption (e.g. heatmap = top 25 up + 25 down by padj), so the LLM describes real figures. `pca_color_by` lets the agent choose extra PCA colourings from informative design columns only (identifier/constant columns are excluded; `condition` is accepted and ignored since it's always `pca.png`). `write_report` requires every generated figure to be placed, fixes image references losslessly (missing `figures/` prefix, or `![figures/x.png]` without the `(path)` part when it names a real figure), rejects references to figures that don't exist (nothing written), inserts a code-generated caption with the file path under each figure, and appends the standard disclaimer (an LLM-written disclaimer is left in place). `generate_report` remains only as a back-compat wrapper for old replay scripts.
- `summarize_findings` returns `data_source` ("nf-core/rnaseq", "GEO count matrix" or "user-provided") and `data_provenance` from files on disk (quantification method, or GEO file/source/value type/MD5), plus the filtering settings and enrichment inputs actually used, and `facts` — every citable value, rendered by `_facts()`, keyed by placeholder name: DE counts/percentages, QC (library sizes and genes detected from the pre-filter `compute_qc` snapshot; PCA variance and minimum sample correlation from the same computation as the figures), filtering, each enrichment term with its statistics (`enrichment.up.go_bp.1`), provenance, the reference genome/annotation nf-core actually used (`reference.*`, from `results/pipeline_info/params_*.json`, following the cache's `reference.json` `gtf_source`), PCA replicate spread and group separation on PC1–PC2 (`qc.pca_spread.<condition>`, `qc.pca_condition_r2_pc1`, `qc.pca_misclustered`, `qc.pca_summary`, also in the PCA caption; compared within each level of the DESeq2 covariates, with `qc.pca_r2_pc1.<covariate>`, so a dominant tissue effect doesn't flag half the samples as misclustered), versions rendered with the tool name (`PyDESeq2 0.5.4`; a name the LLM already wrote before the placeholder isn't repeated), and design columns (`design.<col>` for values shared by all samples such as time point, `design.<col>_values` for levels of varying ones).
- **The LLM never types values into the report.** `write_report` renders `{{fact.name}}` (from the same `_facts()`, so what the LLM reads and what gets printed can't differ), `{{gene:SYMBOL}}` (log2FC and padj from results), `{{cite:PMID}}` (only abstracts fetched this session) and `{{table:name}}` (code-generated tables: `qc`, `de_summary`, `top_up`, `top_down`, `enrichment_up`, `enrichment_down`, `versions`, `methods`). **The LLM never writes Methods:** `{{table:methods}}` is required and generated from what the tools did (data source, reference annotation, QC definitions, filtering, the DESeq2 model/test/rounding, enrichment selection and background, versions); a Methods heading containing anything else is rejected. Tables required by the analysis state must be placed. Any number typed in the prose or an LLM-written table is rejected (`_stray_numbers`, ported from agentic-scrna-workflow; names like CXCL8/GSE164073/PC1/GO IDs and digit-containing design labels such as "MOI = 1.0" are allowed), as is a gene named in the same sentence as an enrichment term placeholder that isn't one of that term's genes (`_term_gene_problems`; `{{gene:X}}` is exempt). Lossless fixes: a `%` or `.` typed after a fact that already ends in one is dropped, and `{{fact.x}}` is read as `{{x}}`. `{{table:x}}` must be on a line of its own (a table inside a sentence or bullet breaks it). Any problem rejects the whole report with every issue listed; after two rejections it is written with visible ⚠ markers and a warning. The LLM still writes every section and all prose. `fetch_abstract` is replayed (not skipped) so citations render in replays.
- Enrichment uses gseapy/Enrichr (network call) — mocked in tests
- PCA via numpy SVD on log2(counts+1), no scanpy dependency

## Troubleshooting and warnings

- On pipeline failure: conversational troubleshooting — LLM (Sonnet) reads the full log, explains the issue, proposes a fix via `propose_fix` tool. User can chat, ask questions, or redirect before approving. Max 3 attempts. All logged to `troubleshooting.jsonl`. Watchdog kills hung nextflow processes (SIGTERM → 10s grace → SIGKILL) so troubleshooting isn't blocked.
- On success: LLM scans log for WARN lines and produces a concise summary of anything affecting downstream analysis.
- Samplesheet agent cites GEO/SRA URLs when assigning conditions, with instructions on where to verify.

## Test datasets

- `examples/GSE164073/`, `examples/GSE157852/` — curated copies of counts-mode smoke runs (2026-09-30, analysis on Sonnet 5): prompt, design, provenance, both agents' tool logs, usage, report + figures, DE/enrichment results and `replay.py` (data paths repointed to `examples/<GSE>/` — the only edit; logs unchanged). Replaying either reproduces its report and `de_results.csv` exactly. Refresh by re-running the smoke test and copying the run folder without `geo/` and `run_state.json`.
- `examples/smoke/` — the three smoke-test prompts (below). Older ad-hoc examples (GSE245856, demo/yeast/Drosophila FASTQs) were moved out of git to `data/legacy_examples/` (gitignored, local data kept).

## Current status (2026-09-23)

**Goal right now:** validate the workflow end to end on real data, then build a benchmark of varied GEO datasets.

**Smoke tests** (`examples/smoke/`, README has commands, expected approval screens and post-run checks):

1. `GSE157852_counts` (counts mode, single factor) — **passing** on Haiku and Sonnet.
2. `GSE164073_counts` (counts mode, hidden `tissue` covariate) — **passing**; both models add `covariates=["tissue"]` unprompted.
3. `GSE246386_full` (full pipeline: download → samplesheet → nextflow salmon-only → analysis, 6 paired-end runs, 5.49 GB) — **in progress**. Download (12 files, MD5-verified) and samplesheet (approved; names match GEO titles) are done. The first samplesheet attempt exposed the `save_samplesheet` schema/signature mismatch and the LLM path relay (both fixed), and the first submission set `max_memory` 32 GB against Docker's 23.4 GB (prompt and tool check fixed). Next: delete the run's `params.json`, `--resume`, confirm `save_reference: true` in `nf_params.yml`, approve. This first nextflow run builds the GRCh38 Salmon index and should populate `data/reference/GRCh38/nf-core-rnaseq-3.26.0/`; the next Salmon-only run is the first real test of reuse.

**Model comparison** (same decisions on both datasets): Haiku ~$0.09–0.12 per counts-mode dataset, Sonnet ~$0.28–0.31. Sonnet writes longer reports and types more numbers by hand; Haiku made one rounding slip before the relevant facts existed. Decision: Haiku for the counts/download/samplesheet/submission agents, Sonnet 5 for the analysis agent (report prose quality; ~$0.35–0.45 per counts-mode run, not yet measured); samplesheet-agent model still open.

**Next steps**

- Finish and review the GSE246386 full run: nextflow, reference caching, post-run warnings, analysis.
- Download approval screen prints the whole ~400-line script and only shows the run → GSM → title table when a subset is selected; it should always show the labelled table, sizes and tool, and just give the script path.
- Build the benchmark — see **Benchmark plan** below.
- **To do (parked): review agent** — see **Review agent plan** below.

**Known gaps (not blocking)**

- Counts mode step 2/3: per-sample GEO files (`_RAW.tar`, e.g. GSE245856) and xlsx are listed but not parsed.
- Only one DE contrast is tracked per run (`SESSION.deseq_results` is overwritten).
- Analysis `set_design` has no approval step.
- Enrichr and NCBI rate limits (429s, captcha pages) are handled but can still slow or interrupt runs; ENA download speed varies a lot (aria2c `-x 16` may help).
- Nextflow `work/` directories are never cleaned up automatically.

## Review agent plan (parked, not built)

**Why:** code now guarantees every number, table, figure caption and the Methods section, but LLM prose still contains wrong claims (invented experimental details, citations for things the abstract doesn't say, gene roles from memory, going further than the paper). Code can't check prose, and adding code rules just moves the text elsewhere (e.g. rejecting "Filtering" headings), so we don't.

**Design:** a separate Sonnet agent runs after `write_report` succeeds. It never sees the analysis agent's reasoning. It checks every claim against `summarize_findings` facts, the generated `{{table:methods}}` block (the answer key for any methods talk elsewhere), `paper_genes`, GEO metadata and fetched abstracts, and it has its own `query_genes` / `search_enrichment` lookups. It reports each problem with a tool call: exact quote, type (contradicts data / unsupported / invented detail / wrong citation / claims more than the paper / methods disagree with Methods), and evidence. Code rejects quotes that aren't in the report, puts a ⚠ note next to each flagged sentence (deleting nothing), and adds a short "Review" summary at the end.

**Phase 1 — flag only (decided):** flags and ⚠ notes, no changes to the text. Measure the reviewer first: real errors caught and false flags raised on the test case below, then as a benchmark metric.

**Phase 2 — targeted fixes (later, once flags are trustworthy):** for each flag the analysis agent replaces the quoted sentence, removes it, or disputes the flag with a reason; code applies the change to that sentence only; the reviewer re-checks only the changed sentences; one round, and anything left gets a ⚠. Not whole-report rewrites (they change unflagged text), and not reviewer-written fixes (nobody checks the checker).

**Logging (both phases):** nothing is overwritten — `analysis/report_v1.md`, `review_v1.json`, `revisions.jsonl` (phase 2: flag → replaced/removed/disputed, old/new text, reason), `review_v2.json`, final `report.md`. Revisions are saved as text-replacement steps so `replay.py` reproduces the final report without either LLM. Reviewer cost goes in `usage.jsonl` (estimate $0.15–0.30 per report on Sonnet).

**First test case:** `runs/20260923_GSE246386_full` second report (15:32). Known errors: "aligned to the reference genome" (Salmon-only); integrins "down" (9 of 14 significant are up); BCR signalling "a hallmark of hematopoietic commitment" (the term is driven by SYK/BTK/PRKCB, which also drive "Platelet activation"); "sufficient to reprogram"/"recapitulates EHT" (abstract: "partial"); invented "24-hour window" and "absence of LSD1 activity"; ALOX15/IGFBP5/CCDC80/TGFBI "implicated in hematopoietic specification" and GIMAP4 as an endothelial marker (from memory); interferon genes as "markers of hematopoietic priming" (uncited); "minimal batch effects" (no batch information); counts called "reads"; a rewritten filtering rule in its own section.

## Benchmark plan (not built yet)

**Purpose:** show that the workflow is reliable on varied real GEO data (the evidence academic labs and small biotechs need before trusting it), find what to fix next, compare models on accuracy and cost, and serve as a regression suite after changes.

**Datasets:** ~8–10 human/mouse GEO series, chosen to cover the known hard cases, each tagged with what it tests:

- FASTQ path: paired vs single-end, multi-lane samples, stranded libraries, a recent series resolvable only via BioProject, runs with blank ENA labels, a subset request ("2 per condition"), a SuperSeries.
- Counts mode: clean author matrix, featureCounts with annotation columns, R-style headers, NCBI-generated only, normalised-only (must stop and explain), per-sample `_RAW.tar` (must stop cleanly until supported).
- Design: >2 groups needing a specific contrast, a hidden covariate (tissue/batch), a donor-paired design, a confounded design (must drop the covariate and say why), mouse.
- Must-stop cases: ambiguous conditions (agent should ask), unlabelled runs (hard stop).

GSE157852, GSE164073 and GSE246386 from `examples/smoke/` are the first entries.

**Expected answers** per dataset (`expected.json`): mode, runs or SubSeries chosen, count file, column → GSM mapping, condition field/labels, covariates, contrast (test/reference), samplesheet rows (sample → FASTQ pair, strandedness), and whether the correct outcome is to stop. Keep expected answers where agents cannot read them — restrict the samplesheet agent's `read_file`/`list_directory` to the run and data directories before building this.

**How it runs (never auto-approves):** each stage is scored on the artifact it proposes at its approval point, then the harness stops — stopping is not approving. Later stages are scored separately from gold inputs (expected `samplesheet.csv`, `counts.tsv` + `design.csv`), so an error in one stage doesn't cascade and each agent is measured in isolation. Stages 0–1 and counts mode need no AWS; run nextflow on only a few datasets.

**Scoring** comes from artifacts and tool logs, never from report prose: `download_metadata.json`, `counts_metadata.json`, `design.csv`, `samplesheet.csv`, and `run_deseq2`/`run_enrichment` arguments in `analysis/tool_calls.jsonl`. Pass/fail per check, plus report health (`write_report` rejections, ⚠ markers, typed-number rejections) and cost/turns from `usage.jsonl`. Spot-check report quality by hand.

**Output:** a results table (dataset × check × model, plus cost) committed under `benchmark/`; raw runs stay gitignored.

**Before building:** restrict file-reading tool paths (above).

## Future direction

- **Streamlit UI** — wrap the CLI in a web app for non-coders. Local mode (user has nextflow/Docker), with cloud submission as a later addition.
- **Multi-provider LLM support** — abstract the agent loop to support OpenAI alongside Anthropic. Thin adapter layer over the current `loop.py` pattern.
- **scRNA-seq workflow** — separate repo, using nf-core for preprocessing with custom downstream analysis. Planned port of the agentic pattern.

## Don't

- Don't add approval logic inside agent tool calls
- Don't auto-approve anything — all approval requires interactive human input
- Don't trust the LLM to relay values (paths, filenames) faithfully between tools — enforce in the tool code
- Don't commit large data files (count matrices, FASTQ, results) — they belong in `.gitignore`. Exception: the curated example runs under `examples/GSE*/` (a few MB each: counts, DE results, figures), so readers can see full outputs and `replay.py` runs from a clone
- Don't make changes without asking first

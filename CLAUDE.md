# agentic-rnaseq-workflow

Raw Anthropic API agent loop (no frameworks) that builds nf-core/rnaseq sample sheets from FASTQ directories, then submits pipeline runs. Tools return structured summary dicts the LLM reasons over.

## Architecture

Four stages, with human approval **between** stages (not inside the agent loop):

0. **Download agent** (`agents/download/`) — resolves GEO/SRA accessions via NCBI/ENA APIs, generates a download script (aspera > aria2c > curl), validates MD5 checksums. Runs whenever the prompt contains an accession (it checks existing files by MD5 and only downloads what's missing); `--skip-download` bypasses it. Rejecting the download script exits. GSE resolution falls back to the series SOFT record's BioProject/SRA link when GEO's search record has none (most series since ~2023), and labels every run with its GSM and title from GEO via the SRX (ENA often leaves these blank). SuperSeries are detected and their SubSeries listed; NCBI/ENA outages are reported as `lookup_failed`, not "no runs". The agent can download a subset of runs (`runs=`). **Hard stop:** runs without a sample ID and title are never downloaded or passed on (`unlabelled_runs`, enforced in the tool and in `run.py`). The selected runs' labels are written to `sample_metadata.csv` and handed to the samplesheet agent.
1. **Samplesheet agent** (`agents/samplesheet/`) — scans FASTQs, reads metadata, matches pairs, drafts/validates/saves a sample sheet, writes a report. Fully implemented and tested with real GEO data.
2. **Submission agent** (`agents/submission/`) — reads the prompt and configures nextflow params (genome, skip_alignment, etc.) via an open-ended tool. Has a `run_command` tool to check system resources (RAM, CPUs) before configuring — always sets resource limits based on the actual machine. Writes `run_nextflow.sh` + `nf_params.yml` + `custom.config`. Pipeline params go in the YAML; resource limits (`max_memory`, `max_cpus`, `max_time`) go in `custom.config` as nextflow `resourceLimits` to avoid nf-schema validation warnings. Human approves the script before execution.
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

`--analyze <run_dir>` attaches to an existing run directory (via `resume_run`), reads the original `prompt.txt` as context, and runs the analysis agent (with or without `results/`). Output goes into the same run dir (`analysis/`), not a new one. `--analyze <prompt.txt>` creates a new analysis-only run (`mode: "analysis"` in `run_state.json`, so `--resume` never enters stages 1-2). If the prompt has a GEO accession, the download agent runs in **counts mode** (`agents/download/counts.py`): `list_geo_count_sources` → `preview_geo_file` → `fetch_geo_counts` → `save_geo_design`. It prefers author supplementary files with raw counts, falling back to NCBI-generated counts; the LLM picks the file, gene ID column and column→GSM mapping, the tools parse and validate (one column per GSM, numeric, non-negative, duplicate IDs summed, samples renamed to GEO titles, symbols from NCBI `gene_info` cached in `data/reference/`). The counts agent stops after the summary — analysis is a separate stage. Writes `counts.tsv`, `counts_metadata.json` (provenance incl. MD5) and `design.csv` (condition, varying GEO characteristics as covariates, and constant ones such as time point kept as facts); the human approves the mapping table before analysis. Raw GEO files are kept in `geo/`. Per-sample files (`_RAW.tar`) and xlsx are listed but not yet parsed. With `--skip-download` (or no accession) the prompt's local count matrix is used. Each analysis writes `analysis/replay.py`, which re-runs the successful tool calls without the LLM (`python -m agents.analysis.replay <run_dir>` for older runs). An optional interactive prompt lets the user add extra instructions (saved as `analysis_prompt.txt` if provided); Enter continues without them, "skip" aborts.

The prompt file describes the data (FASTQ location, organism, metadata path, strandedness). See `examples/` for working prompts. Run directories are named `runs/<datestamp>_<project>/` (derived from prompt file's parent directory).

## Key design rules

- No frameworks (no LangChain, no CrewAI)
- Tools return `dict[str, Any]` summaries — the LLM never sees raw data
- Human approval belongs between stages, never inside the agent loop
- Keep code minimal — no bloat, no premature abstractions
- Tests are offline (no API calls)
- Tools must enforce correctness, not the LLM — paths are resolved to absolute in the tools themselves (`scan_fastqs`, `draft_samplesheet`), `save_samplesheet` always writes to the canonical run directory location, `SubmissionParams` routes resource limits to `custom.config` vs pipeline params to `nf_params.yml` automatically, `configure_submission` forces `pseudo_aligner: salmon` when `skip_alignment` is set (so counts are always produced), and `generate_figures`/`write_report` supply real figure paths and captions and append the disclaimer rather than trusting the LLM to get filenames or caveats right. Tools may make lossless mechanical corrections to LLM output, and reject anything else with a reason — never silently delete what the LLM wrote. The LLM is a lossy intermediary; don't trust it to preserve values faithfully between tool calls.

## Prerequisites

- Python 3.11+, uv
- `ANTHROPIC_API_KEY` environment variable
- Nextflow + Docker (for pipeline submission — not needed for samplesheet agent)

## Config

`core/config.py` — model name, max tokens/turns, pipeline version, paths. Single source of truth. `MODEL` (Haiku 4.5) for agent loops, `MODEL_SONNET` for troubleshooting and warning review.

Agent loops and the troubleshooter use automatic prompt caching (top-level `cache_control`; Haiku 4.5 only caches prefixes ≥ 4096 tokens). Token usage per agent run is printed and appended to the run's `usage.jsonl`.

## Analysis agent tools

`fetch_geo_metadata` → `fetch_abstract` → `scan_results` → `load_counts` → `inspect_counts` → `set_design` (if needed) → `compute_qc` → `filter_low_counts` → `run_deseq2` → `get_top_genes` → `run_enrichment` → `summarize_findings` → `generate_figures` → `write_report`

- `fetch_geo_metadata` / `fetch_abstract` — NCBI E-utilities for GEO series metadata and PubMed abstracts. Provides cell type, organism, experimental context, and citable references. Called first when a GEO accession is in the prompt.
- `inspect_counts` — detects whether data is raw counts vs normalised (TPM/FPKM/log). Prevents running DESeq2 on pre-normalised data. Also flags Salmon fractional counts (`salmon_fractional: true`) so downstream tools know rounding is safe.
- `run_deseq2` — takes optional `covariates` (design columns such as donor/batch/tissue; categorical): design = `~ covariates + factor`. Refuses confounded designs, unknown levels, constant covariates and non-identifier column names; returns the formula (also in `summarize_findings`). Rounds Salmon fractional counts to integers automatically (only when data looks like raw counts with minor fractional noise, not normalised). Reports `counts_rounded` and `pct_non_integer_before_rounding` in its return dict. Also returns `low_replication_warning` when any group has < 3 replicates — the agent must include this caveat in the report.
- `save_design` (samplesheet agent) writes `design.csv` for traceability; `set_design` (analysis agent) is the fallback when none exists
- `run_enrichment(direction, padj_max=0.05, lfc_min=1.0, max_genes=500)` selects genes itself from the DESeq2 results — the LLM never passes gene lists: padj and |log2FC| filters, split by direction, ranked by padj → |log2FC| → gene ID (deterministic; the same `_ranked` helper feeds the top-genes tables, heatmap and DESeq2 preview), capped, converted to symbols. Saves `analysis/enrichment_input_<label>.txt` (exact genes) and `analysis/enrichment_<label>.csv` (raw Enrichr results). Replay scripts therefore contain selection parameters, not hardcoded genes.
- `summarize_findings` returns `software_versions` (real installed versions), `pipeline_versions` (parsed from nf-core's `nf_core_rnaseq_software_mqc_versions.yml` — exact Nextflow, pipeline, and Salmon versions), and `data_source` ("nf-core" or "user-provided") to prevent hallucination in the Methods section
- `generate_figures` draws every figure before the report is written and returns each path with a factual caption (e.g. heatmap = top 25 up + 25 down by padj), so the LLM describes real figures. `pca_color_by` lets the agent choose extra PCA colourings from informative design columns only (identifier/constant columns are excluded; `condition` is accepted and ignored since it's always `pca.png`). `write_report` requires every generated figure to be placed, fixes image references losslessly (missing `figures/` prefix, or `![figures/x.png]` without the `(path)` part when it names a real figure), rejects references to figures that don't exist (nothing written), inserts a code-generated caption with the file path under each figure, and appends the standard disclaimer (an LLM-written disclaimer is left in place). `generate_report` remains only as a back-compat wrapper for old replay scripts.
- `summarize_findings` returns `data_source` ("nf-core/rnaseq", "GEO count matrix" or "user-provided") and `data_provenance` from files on disk (quantification method, or GEO file/source/value type/MD5), plus the filtering settings and enrichment inputs actually used, and `facts` — every citable value, rendered by `_facts()`, keyed by placeholder name: DE counts/percentages, QC (library sizes and genes detected from the pre-filter `compute_qc` snapshot; PCA variance and minimum sample correlation from the same computation as the figures), filtering, each enrichment term with its statistics (`enrichment.up.go_bp.1`), provenance, versions, and design columns (`design.<col>` for values shared by all samples such as time point, `design.<col>_values` for levels of varying ones).
- **The LLM never types values into the report.** `write_report` renders `{{fact.name}}` (from the same `_facts()`, so what the LLM reads and what gets printed can't differ), `{{gene:SYMBOL}}` (log2FC and padj from results), `{{cite:PMID}}` (only abstracts fetched this session) and `{{table:name}}` (code-generated tables: `qc`, `de_summary`, `top_up`, `top_down`, `enrichment_up`, `enrichment_down`, `versions`). Tables required by the analysis state must be placed. Any problem rejects the whole report with every issue listed; after two rejections it is written with visible ⚠ markers and a warning. The LLM still writes every section and all prose. `fetch_abstract` is replayed (not skipped) so citations render in replays.
- Enrichment uses gseapy/Enrichr (network call) — mocked in tests
- PCA via numpy SVD on log2(counts+1), no scanpy dependency

## Troubleshooting and warnings

- On pipeline failure: conversational troubleshooting — LLM (Sonnet) reads the full log, explains the issue, proposes a fix via `propose_fix` tool. User can chat, ask questions, or redirect before approving. Max 3 attempts. All logged to `troubleshooting.jsonl`. Watchdog kills hung nextflow processes (SIGTERM → 10s grace → SIGKILL) so troubleshooting isn't blocked.
- On success: LLM scans log for WARN lines and produces a concise summary of anything affecting downstream analysis.
- Samplesheet agent cites GEO/SRA URLs when assigning conditions, with instructions on where to verify.

## Test datasets

- `examples/GSE245856/` — VPA treatment in HEK293T cells, 6 samples (3 CTRL, 3 VPA). `design.csv` and `prompt.txt` are in git; `counts.csv` (47k genes) is gitignored and must be generated locally. Used for analysis agent testing.

## Current status (2026-09-23)

**Goal right now:** validate the workflow end to end on real data, then build a benchmark of varied GEO datasets before scaling up runs (AWS / API credits being sourced).

**Smoke tests** (`examples/smoke/`, README has commands, expected approval screens and post-run checks):

1. `GSE157852_counts` (counts mode, single factor) — **passing** on Haiku and Sonnet (runs `runs/20260923_GSE157852_counts_2` / `_3`).
2. `GSE164073_counts` (counts mode, hidden `tissue` covariate) — **passing**; both models add `covariates=["tissue"]` unprompted (`_2` / `_3`).
3. `GSE246386_full` (full pipeline: download → samplesheet → nextflow salmon-only → analysis, 6 paired-end runs, 5.49 GB) — **in progress**: `runs/20260923_GSE246386_full`, FASTQs downloading via aria2c into `data/GSE246386/fastqs/` (ENA ~0.5–0.8 MB/s, ~2–3 h). The run waits at the samplesheet approval. If interrupted: `uv run python run.py --resume runs/20260923_GSE246386_full`. This is the first real test of the new FASTQ path (BioProject resolution, GEO run labels, `sample_metadata.csv` handover, real `--resume`). The 19 MB `SRR26539596_1.fastq.gz` in the data folder is a stale partial that aria2c overwrites when it reaches that file.

**Model comparison** (same decisions on both datasets): Haiku ~$0.09–0.12 per counts-mode dataset, Sonnet ~$0.28–0.31. Sonnet writes longer reports and types more numbers by hand; Haiku made one rounding slip before the relevant facts existed. Decision so far: keep Haiku as default (`AGENT_MODEL=sonnet` to compare); samplesheet-agent model still open.

**Next steps**

- Finish and review the GSE246386 full run (check samplesheet against README §3 before approving).
- Download approval screen prints the whole ~400-line script and only shows the run → GSM → title table when a subset is selected; it should always show the labelled table, sizes and tool, and just give the script path.
- Build the benchmark: ~8–10 human/mouse GEO datasets with expected mapping/design/samplesheet; score agents (and Haiku vs Sonnet) with `usage.jsonl` for cost.

**Known gaps (not blocking)**

- Counts mode step 2/3: per-sample GEO files (`_RAW.tar`, e.g. GSE245856) and xlsx are listed but not parsed.
- Untraced-numbers audit in `write_report` (flag numbers typed in prose that match no fact) — not built.
- Only one DE contrast is tracked per run (`SESSION.deseq_results` is overwritten).
- Samplesheet agent relays `match_pairs` output into `draft_samplesheet` (gated by approval, but still an LLM relay); analysis `set_design` has no approval step.
- Enrichr and NCBI rate limits (429s, captcha pages) are handled but can still slow or interrupt runs; ENA download speed varies a lot (aria2c `-x 16` may help).
- Nextflow `work/` directories are never cleaned up automatically.

## Future direction

- **Streamlit UI** — wrap the CLI in a web app for non-coders. Local mode (user has nextflow/Docker), with cloud submission as a later addition.
- **Multi-provider LLM support** — abstract the agent loop to support OpenAI alongside Anthropic. Thin adapter layer over the current `loop.py` pattern.
- **scRNA-seq workflow** — separate repo, using nf-core for preprocessing with custom downstream analysis. Planned port of the agentic pattern.

## Don't

- Don't add approval logic inside agent tool calls
- Don't auto-approve anything — all approval requires interactive human input
- Don't trust the LLM to relay values (paths, filenames) faithfully between tools — enforce in the tool code
- Don't commit large data files (count matrices, FASTQ, results) — they belong in `.gitignore`
- Don't make changes without asking first

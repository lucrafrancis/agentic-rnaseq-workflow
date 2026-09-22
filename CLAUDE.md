# agentic-rnaseq-workflow

Raw Anthropic API agent loop (no frameworks) that builds nf-core/rnaseq sample sheets from FASTQ directories, then submits pipeline runs. Tools return structured summary dicts the LLM reasons over.

## Architecture

Four stages, with human approval **between** stages (not inside the agent loop):

0. **Download agent** (`agents/download/`) — resolves GEO/SRA accessions via NCBI/ENA APIs, generates a download script (aspera > aria2c > curl), validates MD5 checksums. Only runs when the prompt contains an accession but no FASTQ path.
1. **Samplesheet agent** (`agents/samplesheet/`) — scans FASTQs, reads metadata, matches pairs, drafts/validates/saves a sample sheet, writes a report. Fully implemented and tested with real GEO data.
2. **Submission agent** (`agents/submission/`) — reads the prompt and configures nextflow params (genome, skip_alignment, etc.) via an open-ended tool. Has a `run_command` tool to check system resources (RAM, CPUs) before configuring — always sets resource limits based on the actual machine. Writes `run_nextflow.sh` + `nf_params.yml` + `custom.config`. Pipeline params go in the YAML; resource limits (`max_memory`, `max_cpus`, `max_time`) go in `custom.config` as nextflow `resourceLimits` to avoid nf-schema validation warnings. Human approves the script before execution.
3. **Analysis agent** (`agents/analysis/`) — downstream DE, enrichment, QC, and reporting on nf-core/rnaseq outputs. Python-only (PyDESeq2, gseapy, numpy PCA). Optionally accepts a paper PDF as context via Claude's native base64 content blocks. Fetches GEO metadata and PubMed abstracts for context/citations. Fully implemented and tested.

The full flow in `run.py`: download agent (if needed) → human approves download script → samplesheet agent → human approves sample sheet → submission agent configures params → human approves `run_nextflow.sh` → nextflow runs (with troubleshooting on failure) → post-run warning review → analysis agent runs downstream analysis.

## Running

```bash
uv run python run.py <prompt.txt>                # full pipeline (stages 0-3)
uv run python run.py --resume <run_dir>          # resume from existing run directory
uv run python run.py --analyze <run_dir>         # analysis only (stage 3)
```

`--resume` skips stages whose artifacts already exist: `download_metadata.json` → skip download, `samplesheet.csv` → skip samplesheet agent, `params.json` → skip submission agent, `results/` non-empty → skip nextflow, go straight to analysis.

`--analyze` attaches to an existing run directory (via `resume_run`), reads the original `prompt.txt` as context, and runs the analysis agent. Output goes into the same run dir (`analysis/`), not a new one. An optional interactive prompt lets the user add extra instructions (saved as `analysis_prompt.txt` if provided); Enter continues without them, "skip" aborts.

The prompt file describes the data (FASTQ location, organism, metadata path, strandedness). See `examples/` for working prompts. Run directories are named `runs/<datestamp>_<project>/` (derived from prompt file's parent directory).

## Key design rules

- No frameworks (no LangChain, no CrewAI)
- Tools return `dict[str, Any]` summaries — the LLM never sees raw data
- Human approval belongs between stages, never inside the agent loop
- Keep code minimal — no bloat, no premature abstractions
- Tests are offline (no API calls)
- Tools must enforce correctness, not the LLM — paths are resolved to absolute in the tools themselves (`scan_fastqs`, `draft_samplesheet`), `save_samplesheet` always writes to the canonical run directory location, `SubmissionParams` routes resource limits to `custom.config` vs pipeline params to `nf_params.yml` automatically, `configure_submission` forces `pseudo_aligner: salmon` when `skip_alignment` is set (so counts are always produced), and `generate_report` returns actual figure paths and appends disclaimers rather than trusting the LLM to get filenames or caveats right. The LLM is a lossy intermediary; don't trust it to preserve values faithfully between tool calls.

## Prerequisites

- Python 3.11+, uv
- `ANTHROPIC_API_KEY` environment variable
- Nextflow + Docker (for pipeline submission — not needed for samplesheet agent)

## Config

`core/config.py` — model name, max tokens/turns, pipeline version, paths. Single source of truth. `MODEL` (Haiku 4.5) for agent loops, `MODEL_SONNET` for troubleshooting and warning review.

## Analysis agent tools

`fetch_geo_metadata` → `fetch_abstract` → `scan_results` → `load_counts` → `inspect_counts` → `set_design` (if needed) → `compute_qc` → `filter_low_counts` → `run_deseq2` → `get_top_genes` → `run_enrichment` → `summarize_findings` → `generate_report`

- `fetch_geo_metadata` / `fetch_abstract` — NCBI E-utilities for GEO series metadata and PubMed abstracts. Provides cell type, organism, experimental context, and citable references. Called first when a GEO accession is in the prompt.
- `inspect_counts` — detects whether data is raw counts vs normalised (TPM/FPKM/log). Prevents running DESeq2 on pre-normalised data. Also flags Salmon fractional counts (`salmon_fractional: true`) so downstream tools know rounding is safe.
- `run_deseq2` — rounds Salmon fractional counts to integers automatically (only when data looks like raw counts with minor fractional noise, not normalised). Reports `counts_rounded` and `pct_non_integer_before_rounding` in its return dict. Also returns `low_replication_warning` when any group has < 3 replicates — the agent must include this caveat in the report.
- `save_design` (samplesheet agent) writes `design.csv` for traceability; `set_design` (analysis agent) is the fallback when none exists
- `run_enrichment` takes a `label` param (e.g. "upregulated", "downregulated") — results accumulate across calls, each label gets its own figures
- `summarize_findings` returns `software_versions` (real installed versions), `pipeline_versions` (parsed from nf-core's `nf_core_rnaseq_software_mqc_versions.yml` — exact Nextflow, pipeline, and Salmon versions), and `data_source` ("nf-core" or "user-provided") to prevent hallucination in the Methods section
- `generate_report` returns a `figures` dict with actual file paths — the LLM must use only those, not invent filenames. Appends a hallucination disclaimer automatically (including a gene annotation verification caveat). Strips trailing `---` before appending to prevent doubling.
- Enrichment uses gseapy/Enrichr (network call) — mocked in tests
- PCA via numpy SVD on log2(counts+1), no scanpy dependency

## Troubleshooting and warnings

- On pipeline failure: conversational troubleshooting — LLM (Sonnet) reads the full log, explains the issue, proposes a fix via `propose_fix` tool. User can chat, ask questions, or redirect before approving. Max 3 attempts. All logged to `troubleshooting.jsonl`. Watchdog kills hung nextflow processes (SIGTERM → 10s grace → SIGKILL) so troubleshooting isn't blocked.
- On success: LLM scans log for WARN lines and produces a concise summary of anything affecting downstream analysis.
- Samplesheet agent cites GEO/SRA URLs when assigning conditions, with instructions on where to verify.

## Test datasets

- `examples/GSE245856/` — VPA treatment in HEK293T cells, 6 samples (3 CTRL, 3 VPA). `design.csv` and `prompt.txt` are in git; `counts.csv` (47k genes) is gitignored and must be generated locally. Used for analysis agent testing.

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

# agentic-rnaseq-workflow

Raw Anthropic API agent loop (no frameworks) that builds nf-core/rnaseq sample sheets from FASTQ directories, then submits pipeline runs. Tools return structured summary dicts the LLM reasons over.

## Architecture

Four stages, with human approval **between** stages (not inside the agent loop):

0. **Download agent** (`agents/download/`) — resolves GEO/SRA accessions via NCBI/ENA APIs, generates a download script (aspera > aria2c > curl), validates MD5 checksums. Only runs when the prompt contains an accession but no FASTQ path.
1. **Samplesheet agent** (`agents/samplesheet/`) — scans FASTQs, reads metadata, matches pairs, drafts/validates/saves a sample sheet, writes a report. Fully implemented and tested with real GEO data.
2. **Submission agent** (`agents/submission/`) — reads the prompt and configures nextflow params (genome, skip_alignment, etc.) via an open-ended tool. Has a `run_command` tool to check system resources (RAM, CPUs) before configuring — always sets resource limits based on the actual machine. Writes `run_nextflow.sh` + `nf_params.yml` + `custom.config`. Pipeline params go in the YAML; resource limits (`max_memory`, `max_cpus`, `max_time`) go in `custom.config` as nextflow `resourceLimits` to avoid nf-schema validation warnings. Human approves the script before execution.
3. **Analysis agent** (`agents/analysis/`) — downstream DE, enrichment, QC, and reporting on nf-core/rnaseq outputs. Python-only (PyDESeq2, gseapy, numpy PCA). Optionally accepts a paper PDF as context via Claude's native base64 content blocks. Fully implemented and tested.

The full flow in `run.py`: download agent (if needed) → human approves download script → samplesheet agent → human approves sample sheet → submission agent configures params → human approves `run_nextflow.sh` → nextflow runs (with troubleshooting on failure) → post-run warning review → analysis agent runs downstream analysis.

## Running

```bash
uv run python run.py <prompt.txt>         # full pipeline (stages 1-3)
uv run python run.py --analyze <results>  # analysis only (stage 3)
```

The prompt file describes the data (FASTQ location, organism, metadata path, strandedness). See `examples/` for working prompts. Run directories are named `runs/<datestamp>_<project>/` (derived from prompt file's parent directory).

## Key design rules

- No frameworks (no LangChain, no CrewAI)
- Tools return `dict[str, Any]` summaries — the LLM never sees raw data
- Human approval belongs between stages, never inside the agent loop
- Keep code minimal — no bloat, no premature abstractions
- Tests are offline (no API calls)
- Tools must enforce correctness, not the LLM — paths are resolved to absolute in the tools themselves (`scan_fastqs`, `draft_samplesheet`), `save_samplesheet` always writes to the canonical run directory location, and `SubmissionParams` routes resource limits to `custom.config` vs pipeline params to `nf_params.yml` automatically. The LLM is a lossy intermediary; don't trust it to preserve values faithfully between tool calls.

## Prerequisites

- Python 3.11+, uv
- `ANTHROPIC_API_KEY` environment variable
- Nextflow + Docker (for pipeline submission — not needed for samplesheet agent)

## Config

`core/config.py` — model name, max tokens/turns, pipeline version, paths. Single source of truth. `MODEL` (Haiku 4.5) for agent loops, `MODEL_SONNET` for troubleshooting and warning review.

## Analysis agent tools

`scan_results` → `load_counts` → `set_design` (if needed) → `compute_qc` → `filter_low_counts` → `run_deseq2` → `get_top_genes` → `run_enrichment` → `summarize_findings` → `generate_report`

- `save_design` (samplesheet agent) writes `design.csv` for traceability; `set_design` (analysis agent) is the fallback when none exists
- Enrichment uses gseapy/Enrichr (network call) — mocked in tests
- PCA via numpy SVD on log2(counts+1), no scanpy dependency

## Troubleshooting and warnings

- On pipeline failure: conversational troubleshooting — LLM (Sonnet) reads the full log, explains the issue, proposes a fix via `propose_fix` tool. User can chat, ask questions, or redirect before approving. Max 3 attempts. All logged to `troubleshooting.jsonl`. Watchdog kills hung nextflow processes (SIGTERM → 10s grace → SIGKILL) so troubleshooting isn't blocked.
- On success: LLM scans log for WARN lines and produces a concise summary of anything affecting downstream analysis.
- Samplesheet agent cites GEO/SRA URLs when assigning conditions, with instructions on where to verify.

## Future direction

- **Streamlit UI** — wrap the CLI in a web app for non-coders. Local mode (user has nextflow/Docker), with cloud submission as a later addition.
- **Multi-provider LLM support** — abstract the agent loop to support OpenAI alongside Anthropic. Thin adapter layer over the current `loop.py` pattern.

## Don't

- Don't add approval logic inside agent tool calls
- Don't auto-approve anything — all approval requires interactive human input
- Don't trust the LLM to relay values (paths, filenames) faithfully between tools — enforce in the tool code
- Don't make changes without asking first

# agentic-rnaseq-workflow

Raw Anthropic API agent loop (no frameworks) that builds nf-core/rnaseq sample sheets from FASTQ directories, then submits pipeline runs. Tools return structured summary dicts the LLM reasons over.

## Architecture

Three stages, with human approval **between** stages (not inside the agent loop):

1. **Samplesheet agent** (`agents/samplesheet/`) — scans FASTQs, reads metadata, matches pairs, drafts/validates/saves a sample sheet, writes a report. Fully implemented and tested with real GEO data.
2. **Submission handler** (`agents/submission/`) — builds nextflow command, presents for human approval, submits and streams live output. Wired into `run.py`, tested end-to-end.
3. **Analysis agent** (`agents/analysis/`) — stub only. Next to implement.

The full flow in `run.py`: samplesheet agent → human approves sample sheet → submission builds nextflow command → human approves submission → nextflow runs with live streaming output.

## Running

```bash
uv run python run.py <prompt.txt>
```

The prompt file describes the data (FASTQ location, organism, metadata path, strandedness). See `examples/` for working prompts.

## Key design rules

- No frameworks (no LangChain, no CrewAI)
- Tools return `dict[str, Any]` summaries — the LLM never sees raw data
- Human approval belongs between stages, never inside the agent loop
- Keep code minimal — no bloat, no premature abstractions
- Tests are offline (no API calls)
- Tools must enforce correctness, not the LLM — paths are resolved to absolute in the tools themselves (`scan_fastqs`, `draft_samplesheet`), and `save_samplesheet` always writes to the canonical run directory location. The LLM is a lossy intermediary; don't trust it to preserve values faithfully between tool calls.

## Prerequisites

- Python 3.11+, uv
- `ANTHROPIC_API_KEY` environment variable
- Nextflow + Docker (for pipeline submission — not needed for samplesheet agent)

## Config

`core/config.py` — model name, max tokens/turns, pipeline version, paths. Single source of truth. Currently using Haiku 4.5 for cheaper test runs.

## Future direction

- **Analysis agent** — next to build. Downstream analysis on nf-core/rnaseq outputs (DE, pathway enrichment, QC summary).
- **Streamlit UI** — wrap the CLI in a web app for non-coders. Local mode (user has nextflow/Docker), with cloud submission as a later addition.
- **Multi-provider LLM support** — abstract the agent loop to support OpenAI alongside Anthropic. Thin adapter layer over the current `loop.py` pattern.

## Don't

- Don't add approval logic inside agent tool calls
- Don't auto-approve anything — all approval requires interactive human input
- Don't trust the LLM to relay values (paths, filenames) faithfully between tools — enforce in the tool code
- Don't make changes without asking first

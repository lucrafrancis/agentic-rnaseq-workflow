# agentic-rnaseq-workflow

Raw Anthropic API agent loop (no frameworks) that builds nf-core/rnaseq sample sheets from FASTQ directories, then submits pipeline runs. Tools return structured summary dicts the LLM reasons over.

## Architecture

Three stages, with human approval **between** stages (not inside the agent loop):

1. **Samplesheet agent** (`agents/samplesheet/`) — scans FASTQs, reads metadata, matches pairs, drafts/validates/saves a sample sheet, writes a report. Fully implemented.
2. **Submission handler** (`agents/submission/`) — builds nextflow command, presents for human approval, submits and monitors. Skeleton implemented, not yet wired into `run.py`.
3. **Analysis agent** (`agents/analysis/`) — stub only.

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

## Prerequisites

- Python 3.11+, uv
- `ANTHROPIC_API_KEY` environment variable
- Nextflow + Docker (for pipeline submission — not needed for samplesheet agent)

## Config

`core/config.py` — model name, max tokens/turns, pipeline version, paths. Single source of truth.

## Don't

- Don't add approval logic inside agent tool calls
- Don't auto-approve anything — all approval requires interactive human input
- Don't make changes without asking first

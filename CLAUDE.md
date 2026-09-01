# Instructions: Set up agentic-rnaseq-workflow

Set up a new repo called `agentic-rnaseq-workflow`. Follow the same philosophy as my existing `agentic-scrna-workflow` repo (https://github.com/lucrafrancis/agentic-scrna-workflow) — raw Anthropic API loop, no framework, tools that return summaries the model reasons over. Use uv + pyproject.toml for dependency management.

## Repo structure

```
agentic-rnaseq-workflow/
├── agents/
│   ├── samplesheet/       # Stage 1: sample sheet generation agent
│   │   ├── loop.py        # agent loop (talks to Claude, runs tools)
│   │   ├── tools.py       # tools: scan files, read metadata, draft sheet, etc.
│   │   ├── schemas.py     # tool descriptions + name-to-function map
│   │   └── prompts.py     # system prompt for the sample sheet agent
│   ├── submission/        # Stage 2: nf-core/rnaseq submission handler
│   │   ├── submit.py      # builds nextflow command, submits, monitors
│   │   ├── params.py      # default params + prompt-based overrides
│   │   └── errors.py      # parse nextflow logs, diagnose common failures
│   └── analysis/          # Stage 3: post-pipeline analysis agent (stub for now)
│       └── __init__.py
├── core/
│   ├── config.py          # model name, paths, shared settings
│   ├── approval.py        # human-in-the-loop approval flow (present draft, wait for confirm/edit)
│   └── session.py         # run directory management, state between stages
├── tests/
│   ├── fixtures/          # example FASTQs (tiny/empty), metadata files, expected sample sheets
│   ├── test_samplesheet.py
│   └── test_submission.py
├── run.py                 # entry point — orchestrates the three stages
├── pyproject.toml
├── .gitignore
├── .python-version
├── LICENSE (MIT)
└── README.md
```

## What to implement now vs. stub

**Implement:** repo scaffolding, all files with module docstrings and class/function signatures. `core/` modules with basic working logic (config loading, run directory creation, the approval flow that presents a draft and waits for user input). The samplesheet `tools.py` should have function signatures with docstrings describing what each tool does — `scan_fastqs(directory)`, `read_metadata(filepath)`, `match_pairs(fastq_list)`, `draft_samplesheet(matches, metadata)`, `validate_samplesheet(sheet)` — but implementation can be minimal stubs.

**Stub only:** `agents/analysis/`, the actual agent loop logic in samplesheet (we'll build this next).

## Key design decisions

- Each tool returns a structured summary string (not raw data) that the agent reads to decide next steps
- The sample sheet agent presents: draft sheet (first 10 rows), reasoning for its decisions, path to full CSV, and summary stats (total samples, pairs found, any warnings)
- Submission handler has two modes: default (sensible nf-core/rnaseq defaults) and prompt-driven (user provides instructions, agent parses them and may make additional decisions like matching genome references to existing data)
- Human approval is required before sample sheet is finalised and before nf-core submission in prompt-driven mode
- Tests should be offline by default (no API calls) — mock the LLM responses, test the tools and parsing logic directly

## Don't

- Don't implement the actual Claude API loop yet
- Don't add a framework (no LangChain, no CrewAI)
- Don't implement the analysis agent beyond a stub
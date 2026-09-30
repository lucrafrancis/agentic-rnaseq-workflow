# agentic-rnaseq-workflow

An LLM agent workflow that turns a short plain-English prompt ("GSE164073, compare
SARS-CoV-2 infected against mock") into a complete bulk RNA-seq analysis. It finds the
data on GEO/SRA, builds the nf-core/rnaseq sample sheet, configures and runs the
pipeline, and writes a differential expression report with figures, enrichment and
citations.

It is written against the raw Anthropic API with no agent framework. The design goal is
simple: **the LLM makes the judgement calls, and code makes sure every value is correct.**
A person approves each stage's output before the next stage starts.

<img src="examples/GSE164073/analysis/figures/pca_tissue.png" width="480" alt="PCA of GSE164073 coloured by tissue">

*From the [GSE164073 example](examples/GSE164073/analysis/report.md). The prompt only
asked to compare infected against mock. The agent read the GEO metadata, saw that the 18
samples come from three eye tissues (which dominate the PCA), and fitted
`~tissue + condition` without being asked.*

## How it works

1. **Write a prompt**: a few lines naming the GEO accession, organism and the comparison
   you want, plus anything extra the agents should know (see
   [`examples/smoke/`](examples/smoke/README.md) and [Usage](#usage)).
2. **Run it.** The agents take it from there, stopping for your approval between stages.
3. **Before the analysis starts**, you can type extra instructions and give the path to a
   PDF of the paper (both optional). The analysis agent reads the PDF for context. It also
   fetches the GEO series metadata and the linked paper's PubMed abstract itself; only
   fetched abstracts can be cited in the report.

```
prompt.txt
   │
   ▼
0. Download agent ──▶ GEO/SRA/ENA lookup, download script, MD5 checks
   │                    ── human approves download script ──
   ▼
1. Samplesheet agent ─▶ scan FASTQs, pair reads, label samples from GEO
   │                    ── human approves sample → condition → run → GSM table ──
   ▼
2. Submission agent ──▶ check machine resources, write nf-core/rnaseq params + config
   │                    ── human approves run_nextflow.sh ──
   ▼
   nextflow run nf-core/rnaseq   (on failure: an LLM troubleshooter proposes fixes)
   │
   ▼
   │                    ── optional: extra instructions, paper PDF ──
   ▼
3. Analysis agent ────▶ GEO metadata + PubMed abstract, QC, PyDESeq2, Enrichr
                         enrichment, figures, report.md
```

There is also a **counts mode** (`--analyze`). It skips FASTQ processing: the agent picks
the author's count matrix on GEO, maps each column to a GSM and builds a design table.
You approve the mapping, and then the analysis runs.

## Design principles

Most of the code exists because an LLM is a lossy intermediary. It is good at deciding
*what* to do and unreliable at copying values from one step to the next.

- **Tools return summaries, not data.** Each tool returns a small dict for the LLM to
  reason over. The LLM never sees a count matrix or a FASTQ list.
- **The LLM never relays values.** It refers to FASTQ pairs by ID, not by path. It picks
  DE thresholds rather than gene lists. The code keeps track of files and paths.
- **The LLM never types a number into the report.** It writes prose with placeholders
  (`{{de.n_up}}`, `{{gene:SOD2}}`, `{{cite:34022129}}`, `{{table:top_up}}`), and
  code fills them in from the actual results. Code also generates the Methods section
  from what the tools actually did. A report is rejected and sent back to the agent, with
  every problem listed, if it types a number itself, cites a figure, gene, paper or fact
  that doesn't exist, or names a gene as part of an enrichment term it isn't in.
- **Tools enforce correctness.** Examples: the sample sheet cannot be saved until it
  passes validation, confounded DESeq2 designs are refused, runs without a GEO sample
  label are never downloaded, and memory limits above what Docker has are rejected.
  Tools may make lossless corrections to LLM output and reject anything else with a
  reason; they never silently drop what the LLM wrote.
- **Humans approve between stages, never inside the agent loop.** Nothing is
  auto-approved.
- **Reproducible.** Every tool call is logged, and every analysis writes a `replay.py`
  that reruns it without the LLM. `--resume` continues an interrupted run.

## Examples

Both examples were run in counts mode. Each folder holds the prompt, the design table,
the provenance of the GEO file, both agents' tool logs, token usage, the report with its
figures, and a `replay.py`.

| Example | What it shows |
|---|---|
| [GSE164073: SARS-CoV-2 in eye tissue](examples/GSE164073/analysis/report.md) | A hidden covariate: the agent adds `tissue` to the model unprompted, and QC compares replicates within each tissue |
| [GSE157852: SARS-CoV-2 in choroid plexus organoids](examples/GSE157852/analysis/report.md) | Choosing one contrast from three groups; viral transcripts among the top "genes", recognised as such |

The reports are unedited. Numbers, tables, figure captions and Methods come from code; the
prose is the LLM's and can still overstate things (see [Limitations](#status-and-limitations)).

## Setup

Tested on macOS (Apple Silicon) only.

What you need depends on what you want to run:

| To… | You need |
|---|---|
| Replay the examples | Python 3.11+, [uv](https://docs.astral.sh/uv/), network access (Enrichr and PubMed) |
| Analyse a GEO count matrix (counts mode) | + an `ANTHROPIC_API_KEY` |
| Run the full pipeline from FASTQs | + Nextflow, Docker Desktop, and plenty of disk |

```bash
git clone https://github.com/lucrafrancis/agentic-rnaseq-workflow
cd agentic-rnaseq-workflow
uv sync

# Replay an example without the LLM (no API key needed)
uv run python examples/GSE164073/analysis/replay.py

# Counts mode: GEO count matrix → DE report
export ANTHROPIC_API_KEY=...
uv run python run.py --analyze examples/smoke/GSE164073_counts/prompt.txt
```

Replays write to `runs/`. They query Enrichr and PubMed live, so enrichment results can
change if Enrichr's gene-set libraries are updated.

### Full pipeline

Tested with Nextflow 26.04, OpenJDK 26, Docker 29.7 and nf-core/rnaseq 3.26.0 (pinned in
`core/config.py`).

- **Nextflow** (`brew install nextflow`) needs Java 17 or later.
- **Docker Desktop.** The pipeline can only use the memory Docker is given (Settings →
  Resources), not all of your Mac's RAM. The submission agent sets resource limits from
  what Docker actually has.
- **Disk.** Budget for the FASTQs (about 6 GB for the smoke test), the reference and
  Docker images, and Nextflow's `work/` folder, which can grow to several times the input
  size. Delete `work/` once you've checked the results.
- **Faster downloads (optional).** The download script uses Aspera if installed, then
  aria2c (`brew install aria2`), then curl.

```bash
uv run python run.py examples/smoke/GSE246386_full/prompt.txt
```

The first Salmon-only run builds a GRCh38 index and caches it in `data/reference/`; later
runs reuse it. The full path has been run stage by stage on real data; a complete
end-to-end run is in progress.

## Usage

| Command | What it does |
|---|---|
| `run.py prompt.txt` | full pipeline (stages 0–3) |
| `run.py prompt.txt --skip-download` | FASTQs already local |
| `run.py --resume runs/<dir>` | continue an interrupted run |
| `run.py --analyze runs/<dir>` | re-run the analysis on an existing run |
| `run.py --analyze prompt.txt` | analysis only, from a GEO or local count matrix |

A prompt names the data and the comparison. The agents work out the rest (which file,
which samples, the design). For example,
[`examples/smoke/GSE164073_counts/prompt.txt`](examples/smoke/GSE164073_counts/prompt.txt):

```
GSE164073 — SARS-CoV-2 infection of human ocular surface tissue.
Homo sapiens. Start from the processed count matrix on GEO (no FASTQ processing).
Compare SARS-CoV-2 infected against mock.
Run the full downstream analysis: QC, differential expression, enrichment, and report.
```

Add anything else the agents should know in plain words: a subset of samples, which
covariates matter, strandedness, a lower fold-change threshold. Just before the analysis
starts, `run.py` also asks for optional extra instructions and the path to a PDF of the
paper (Enter skips both).

[`examples/smoke/`](examples/smoke/README.md) holds the test prompts, with what each
approval screen should show and what to check afterwards.

## Models and cost

- **Claude Haiku 4.5** runs the download/counts, samplesheet and submission agents.
- **Claude Sonnet 5** runs the analysis agent, which writes the report.
- **Claude Sonnet 4.5** runs the pipeline troubleshooter and the warning review.

Both examples cost about **$0.25** each in API usage (about $0.02 for the counts agent,
the rest for the analysis). Every run records its token use and estimated cost per agent
in `usage.jsonl`, using the prices in `core/config.py`. Adding a paper PDF raises the
analysis cost, since the whole PDF is sent with every request (cached after the first).
Set `ANALYSIS_MODEL=haiku` to run the analysis agent on Haiku instead (about $0.13 per
run), or `AGENT_MODEL` to change the other agents.

## Project layout

```
run.py                  stage orchestration, approvals, resume
core/                   agent loop, session state, approval prompts, config
agents/download/        GEO/SRA/ENA resolution, download script, counts mode
agents/samplesheet/     FASTQ scanning, pairing, sample sheet validation
agents/submission/      nf-core params, resource checks, reference cache, troubleshooter
agents/analysis/        QC, DESeq2, enrichment, figures, report rendering, replay
examples/               finished example runs, and smoke-test prompts
tests/                  offline test suite (no API calls; network mocked)
```

```bash
uv run pytest
```

## Status and limitations

This is a working prototype, not a validated clinical or production tool. Check the
results before relying on them.

- Counts mode has been tested end to end on real GEO series. The full FASTQ → nf-core
  path has been exercised stage by stage on real data; a complete end-to-end run is in
  progress.
- Code guarantees the numbers, tables, figure captions and Methods, but not the prose.
  The LLM can still overstate or invent claims in the interpretation, e.g. calling counts
  "reads".
- Only one DE contrast per run. Per-sample GEO files (`_RAW.tar`) and xlsx count
  tables are not parsed yet.
- Enrichment uses the Enrichr web API, so it depends on the network and on Enrichr's
  rate limits.

**Planned:**
- A separate review agent that flags unsupported claims in the report prose.
- A benchmark of about 10 varied GEO datasets with expected answers, scored per stage
  and per model.
- A Streamlit UI for people who don't use the command line.

## License

MIT

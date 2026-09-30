"""Central configuration. Everything tunable lives here so nothing else has to change.

The model name in particular is intentionally a single constant: swapping models is a
one-line edit and nothing else in the codebase depends on it.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

# --- LLM ---------------------------------------------------------------------
MODEL_SONNET = "claude-sonnet-4-5-20250929"
# Agent-loop model. AGENT_MODEL overrides it for one run, e.g. AGENT_MODEL=sonnet to
# compare against the default. Recorded per agent in the run's usage.jsonl.
_MODEL_ALIASES = {"haiku": "claude-haiku-4-5-20251001", "sonnet": MODEL_SONNET, "sonnet5": "claude-sonnet-5"}
_requested = os.environ.get("AGENT_MODEL") or "haiku"
MODEL = _MODEL_ALIASES.get(_requested, _requested)
# The analysis agent writes the report, where prose quality matters most, so it runs on
# Sonnet 5 (as in agentic-scrna-workflow). ANALYSIS_MODEL=haiku overrides it for one run.
_requested_analysis = os.environ.get("ANALYSIS_MODEL") or "sonnet5"
ANALYSIS_MODEL = _MODEL_ALIASES.get(_requested_analysis, _requested_analysis)
# USD per million tokens, for the cost estimate in each run's usage.jsonl (Anthropic
# first-party API rates). Cache writes (5-minute TTL) cost 1.25x the input price and cache
# reads 0.1x. A model missing here is logged without a cost.
PRICE_PER_MTOK = {
    "claude-haiku-4-5-20251001": {"input": 1.00, "output": 5.00},
    "claude-sonnet-5": {"input": 2.00, "output": 10.00},
    MODEL_SONNET: {"input": 3.00, "output": 15.00},  # Sonnet 4.5
}
CACHE_WRITE_MULTIPLIER = 1.25
CACHE_READ_MULTIPLIER = 0.1
MAX_TOKENS = 16384
MAX_TURNS = 40

# --- Paths -------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
RUNS_DIR = ROOT / "runs"
REFERENCE_DIR = ROOT / "data" / "reference"  # shared cache (e.g. NCBI gene_info), gitignored

# --- nf-core/rnaseq defaults ------------------------------------------------
NFCORE_PIPELINE = "nf-core/rnaseq"
NFCORE_REVISION = "3.26.0"


@dataclass(frozen=True)
class RunPaths:
    """Per-run output locations, namespaced by a timestamped run name.

    A run on project X at 14:30 writes everything under runs/<timestamp>_<name>/,
    so runs never collide. The Session derives the name and holds an instance of this.
    """

    name: str

    @property
    def dir(self) -> Path:
        return RUNS_DIR / self.name

    @property
    def state_file(self) -> Path:
        return self.dir / "run_state.json"

    @property
    def samplesheet(self) -> Path:
        return self.dir / "samplesheet.csv"

    @property
    def nextflow_log(self) -> Path:
        return self.dir / "nextflow.log"

    @property
    def usage_log(self) -> Path:
        """Token usage per agent run — for tracking API spend."""
        return self.dir / "usage.jsonl"

    @property
    def tool_log(self) -> Path:
        return self.dir / "tool_calls.jsonl"

    @property
    def params_file(self) -> Path:
        return self.dir / "params.json"

    @property
    def nextflow_script(self) -> Path:
        return self.dir / "run_nextflow.sh"

    @property
    def nf_params(self) -> Path:
        return self.dir / "nf_params.yml"

    @property
    def nf_config(self) -> Path:
        return self.dir / "custom.config"

    @property
    def troubleshooting_log(self) -> Path:
        return self.dir / "troubleshooting.jsonl"

    @property
    def download_script(self) -> Path:
        return self.dir / "download.sh"

    @property
    def download_metadata(self) -> Path:
        return self.dir / "download_metadata.json"

    @property
    def geo_dir(self) -> Path:
        """Raw GEO supplementary files downloaded in counts mode (kept for provenance)."""
        return self.dir / "geo"

    @property
    def geo_sources(self) -> Path:
        return self.dir / "geo_sources.json"

    @property
    def counts_matrix(self) -> Path:
        return self.dir / "counts.tsv"

    @property
    def counts_metadata(self) -> Path:
        return self.dir / "counts_metadata.json"

    @property
    def sample_metadata(self) -> Path:
        """Run -> GSM -> title for downloaded runs, handed to the samplesheet agent."""
        return self.dir / "sample_metadata.csv"

    @property
    def design(self) -> Path:
        return self.dir / "design.csv"

    @property
    def analysis_dir(self) -> Path:
        return self.dir / "analysis"

    @property
    def analysis_figures(self) -> Path:
        return self.analysis_dir / "figures"

    @property
    def analysis_report(self) -> Path:
        return self.analysis_dir / "report.md"

    @property
    def de_results(self) -> Path:
        return self.analysis_dir / "de_results.csv"

    @property
    def analysis_tool_log(self) -> Path:
        return self.analysis_dir / "tool_calls.jsonl"

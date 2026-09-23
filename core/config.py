"""Central configuration. Everything tunable lives here so nothing else has to change.

The model name in particular is intentionally a single constant: swapping models is a
one-line edit and nothing else in the codebase depends on it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

# --- LLM ---------------------------------------------------------------------
MODEL = "claude-haiku-4-5-20251001"
MODEL_SONNET = "claude-sonnet-4-5-20250929"
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

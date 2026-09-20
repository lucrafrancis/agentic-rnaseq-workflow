"""Run-scoped state shared across stages.

Each run gets a timestamped directory under runs/. The session tracks which stages have
completed and holds the paths that all stages write to. Stages reach shared state through
the module-level SESSION singleton rather than passing it as an argument.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any

from core.config import RunPaths

if TYPE_CHECKING:
    import pandas as pd


class Session:
    def __init__(self) -> None:
        self.name: str | None = None
        self.paths: RunPaths | None = None
        self.fastq_dir: Path | None = None
        self.stages_completed: list[str] = []

        # Analysis state — populated by Stage 3 tools
        self.counts_df: pd.DataFrame | None = None
        self.gene_names: pd.Series | None = None
        self.tpm_df: pd.DataFrame | None = None
        self.design_df: pd.DataFrame | None = None
        self.deseq_results: pd.DataFrame | None = None
        self.enrichment_results: dict[str, Any] | None = None
        self.results_dir: Path | None = None
        self.multiqc_stats: pd.DataFrame | None = None

    def begin_run(self, project_name: str) -> None:
        """Create a new timestamped run directory and initialise paths."""
        datestamp = datetime.now().strftime("%Y%m%d")
        base = f"{datestamp}_{project_name}"
        self.name = base
        self.paths = RunPaths(self.name)
        n = 1
        while self.paths.dir.exists():
            n += 1
            self.name = f"{base}_{n}"
            self.paths = RunPaths(self.name)
        self.paths.dir.mkdir(parents=True, exist_ok=True)
        self.stages_completed = []

    def mark_stage_complete(self, stage: str) -> None:
        if stage not in self.stages_completed:
            self.stages_completed.append(stage)

    def require_paths(self) -> RunPaths:
        if self.paths is None:
            raise RuntimeError("No run started. Call begin_run first.")
        return self.paths


SESSION = Session()

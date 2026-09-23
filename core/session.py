"""Run-scoped state shared across stages.

Each run gets a timestamped directory under runs/. The session tracks which stages have
completed and holds the paths that all stages write to. Stages reach shared state through
the module-level SESSION singleton rather than passing it as an argument.
"""

from __future__ import annotations

import json
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
        self.source_prompt: str | None = None  # where the prompt was read from, for relative paths
        self.mode: str = "pipeline"  # "pipeline" (FASTQ -> nextflow -> analysis) or "analysis" (counts only)

        # Analysis state — populated by Stage 3 tools
        self.counts_df: pd.DataFrame | None = None
        self.gene_names: pd.Series | None = None
        self.tpm_df: pd.DataFrame | None = None
        self.design_df: pd.DataFrame | None = None
        self.deseq_results: pd.DataFrame | None = None
        self.deseq_design: str | None = None
        self.enrichment_results: dict[str, Any] | None = None
        self.results_dir: Path | None = None
        self.multiqc_stats: pd.DataFrame | None = None

    def begin_run(self, project_name: str, source_prompt: Path | None = None, mode: str = "pipeline") -> None:
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
        self.source_prompt = str(source_prompt.resolve()) if source_prompt else None
        self.mode = mode
        self._save_state()

    def resume_run(self, run_dir: Path) -> None:
        """Attach to an existing run directory for resumption."""
        run_dir = run_dir.resolve()
        if not run_dir.is_dir():
            raise FileNotFoundError(f"Run directory not found: {run_dir}")
        from core.config import RUNS_DIR
        if run_dir.parent != RUNS_DIR.resolve():
            raise ValueError(f"Run directory must be under {RUNS_DIR}, got: {run_dir}")
        self.name = run_dir.name
        self.paths = RunPaths(self.name)
        state_file = self.paths.state_file
        state = json.loads(state_file.read_text()) if state_file.is_file() else {}
        self.stages_completed = state.get("stages_completed", [])
        self.source_prompt = state.get("source_prompt")
        self.mode = state.get("mode", "pipeline")

    def mark_stage_complete(self, stage: str) -> None:
        """Record a stage as done (approved/validated). Persisted so --resume can trust it."""
        if stage not in self.stages_completed:
            self.stages_completed.append(stage)
            self._save_state()

    def is_complete(self, stage: str) -> bool:
        return stage in self.stages_completed

    def _save_state(self) -> None:
        state = {"stages_completed": self.stages_completed, "source_prompt": self.source_prompt, "mode": self.mode}
        self.require_paths().state_file.write_text(json.dumps(state, indent=2) + "\n")

    def require_paths(self) -> RunPaths:
        if self.paths is None:
            raise RuntimeError("No run started. Call begin_run first.")
        return self.paths


SESSION = Session()

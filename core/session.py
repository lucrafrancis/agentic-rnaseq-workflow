"""Run-scoped state shared across stages.

Each run gets a timestamped directory under runs/. The session tracks which stages have
completed and holds the paths that all stages write to. Stages reach shared state through
the module-level SESSION singleton rather than passing it as an argument.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

from core.config import RunPaths


class Session:
    def __init__(self) -> None:
        self.name: str | None = None
        self.paths: RunPaths | None = None
        self.fastq_dir: Path | None = None
        self.stages_completed: list[str] = []

    def begin_run(self, project_name: str) -> None:
        """Create a new timestamped run directory and initialise paths."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.name = f"{timestamp}_{project_name}"
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

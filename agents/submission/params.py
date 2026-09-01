"""Default parameters and prompt-based overrides for nf-core/rnaseq.

Sensible defaults for a standard paired-end RNA-seq run. In prompt-driven mode the agent
can override any of these based on user instructions (e.g. a specific genome, aligner,
or output directory).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from core.config import NFCORE_PIPELINE, NFCORE_REVISION


@dataclass
class SubmissionParams:
    """Parameters for an nf-core/rnaseq submission."""

    pipeline: str = NFCORE_PIPELINE
    revision: str = NFCORE_REVISION
    input_samplesheet: str = ""
    outdir: str = ""
    genome: str = "GRCh38"
    profile: str = "docker"
    extra_args: dict[str, Any] = field(default_factory=dict)

    def to_nextflow_args(self) -> list[str]:
        """Build the nextflow run command arguments."""
        args = [
            "nextflow", "run", self.pipeline,
            "-r", self.revision,
            "--input", self.input_samplesheet,
            "--outdir", self.outdir,
            "--genome", self.genome,
            "-profile", self.profile,
        ]
        for key, value in self.extra_args.items():
            if isinstance(value, bool):
                if value:
                    args.append(f"--{key}")
            else:
                args.extend([f"--{key}", str(value)])
        return args

    def to_command_string(self) -> str:
        return " ".join(self.to_nextflow_args())

    def to_dict(self) -> dict[str, Any]:
        return {
            "pipeline": self.pipeline,
            "revision": self.revision,
            "input": self.input_samplesheet,
            "outdir": self.outdir,
            "genome": self.genome,
            "profile": self.profile,
            **self.extra_args,
        }


def default_params(samplesheet_path: str, outdir: str) -> SubmissionParams:
    """Construct params with sensible defaults for a standard run."""
    return SubmissionParams(
        input_samplesheet=samplesheet_path,
        outdir=outdir,
    )

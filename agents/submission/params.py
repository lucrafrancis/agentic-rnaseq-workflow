"""Default parameters and prompt-based overrides for nf-core/rnaseq.

Sensible defaults for a standard paired-end RNA-seq run. In prompt-driven mode the agent
can override any of these based on user instructions (e.g. a specific genome, aligner,
or output directory).
"""

from __future__ import annotations

import json
import shlex
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

    def to_script(self, params_file: str) -> str:
        """Generate a bash script that references a params YAML file."""
        return "\n".join([
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            "",
            f"nextflow run {self.pipeline} \\",
            f"  -r {self.revision} \\",
            f"  --input {shlex.quote(self.input_samplesheet)} \\",
            f"  --outdir {shlex.quote(self.outdir)} \\",
            f"  -profile {self.profile} \\",
            f"  -params-file {shlex.quote(params_file)}",
            "",
        ])

    def write_nf_params(self, path: Path) -> None:
        """Write pipeline params (genome + extras) to a YAML file.

        YAML preserves types — booleans stay booleans, so nf-schema
        won't reject them as strings.
        """
        params: dict[str, Any] = {"genome": self.genome}
        params.update(self.extra_args)
        lines = []
        for key, value in params.items():
            if isinstance(value, bool):
                lines.append(f"{key}: {str(value).lower()}")
            elif isinstance(value, (int, float)):
                lines.append(f"{key}: {value}")
            else:
                lines.append(f"{key}: \"{value}\"")
        path.write_text("\n".join(lines) + "\n")

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


    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "SubmissionParams":
        """Reconstruct from a to_dict() round-trip (e.g. params.json)."""
        _KNOWN = {"pipeline", "revision", "input", "outdir", "genome", "profile"}
        return cls(
            pipeline=data.get("pipeline", NFCORE_PIPELINE),
            revision=data.get("revision", NFCORE_REVISION),
            input_samplesheet=data.get("input", ""),
            outdir=data.get("outdir", ""),
            genome=data.get("genome", "GRCh38"),
            profile=data.get("profile", "docker"),
            extra_args={k: v for k, v in data.items() if k not in _KNOWN},
        )

    @classmethod
    def load(cls, params_file: Path) -> "SubmissionParams":
        """Load from a params.json file."""
        return cls.from_dict(json.loads(params_file.read_text()))


def default_params(samplesheet_path: str, outdir: str) -> SubmissionParams:
    """Construct params with sensible defaults for a standard run."""
    return SubmissionParams(
        input_samplesheet=samplesheet_path,
        outdir=outdir,
    )

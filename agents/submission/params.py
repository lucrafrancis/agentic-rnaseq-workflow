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


_NEXTFLOW_CONFIG_KEYS = {"max_memory", "max_cpus", "max_time"}


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

    @property
    def pipeline_params(self) -> dict[str, Any]:
        """Params that belong in the YAML params file (nf-schema validated)."""
        return {k: v for k, v in self.extra_args.items() if k not in _NEXTFLOW_CONFIG_KEYS}

    @property
    def config_params(self) -> dict[str, Any]:
        """Params that belong in nextflow config (max_memory, max_cpus, max_time)."""
        return {k: v for k, v in self.extra_args.items() if k in _NEXTFLOW_CONFIG_KEYS}

    def to_script(
        self,
        params_file: str,
        config_file: str | None = None,
        resume: bool = False,
    ) -> str:
        """Generate a bash script that references a params YAML file."""
        lines = [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            "",
            f"nextflow run {self.pipeline} \\",
            f"  -r {self.revision} \\",
        ]
        if resume:
            lines.append("  -resume \\")
        lines += [
            f"  --input {shlex.quote(self.input_samplesheet)} \\",
            f"  --outdir {shlex.quote(self.outdir)} \\",
            f"  -profile {self.profile} \\",
        ]
        if config_file:
            lines.append(f"  -c {shlex.quote(config_file)} \\")
        lines.append(f"  -params-file {shlex.quote(params_file)}")
        lines.append("")
        return "\n".join(lines)

    def write_nf_config(self, path: Path) -> None:
        """Write nextflow resource limits config (max_memory, max_cpus, max_time)."""
        cfg = self.config_params
        if not cfg:
            return
        key_map = {"max_memory": "memory", "max_cpus": "cpus", "max_time": "time"}
        entries = []
        for key, value in cfg.items():
            nf_key = key_map.get(key, key)
            if isinstance(value, bool):
                entries.append(f"        {nf_key}: {str(value).lower()}")
            elif isinstance(value, (int, float)):
                entries.append(f"        {nf_key}: {value}")
            else:
                entries.append(f"        {nf_key}: '{value}'")
        lines = [
            "process {",
            "    resourceLimits = [",
            ",\n".join(entries),
            "    ]",
            "}",
            "",
        ]
        path.write_text("\n".join(lines))

    def write_nf_params(self, path: Path) -> None:
        """Write pipeline params (genome + extras) to a YAML file.

        YAML preserves types — booleans stay booleans, so nf-schema
        won't reject them as strings.
        """
        params: dict[str, Any] = {"genome": self.genome}
        params.update(self.pipeline_params)
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

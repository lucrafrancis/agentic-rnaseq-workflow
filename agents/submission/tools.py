"""Submission agent tools — configure nextflow parameters for nf-core/rnaseq.

The LLM reads the user's prompt and freely decides what params to set. The tool
accepts an open-ended dict so nothing is restricted by the Python side.
"""

from __future__ import annotations

import json
from typing import Any

from agents.submission.params import SubmissionParams
from core.config import NFCORE_PIPELINE, NFCORE_REVISION
from core.session import SESSION


def configure_submission(
    genome: str,
    profile: str = "docker",
    extra_params: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Save nextflow submission parameters to params.json.

    The LLM decides genome, profile, and any extra nf-core/rnaseq params
    (skip_alignment, max_memory, max_cpus, etc.) based on the user's prompt.
    Samplesheet path and outdir are set from the session — the LLM doesn't
    need to provide them.
    """
    paths = SESSION.require_paths()

    clean = {}
    for k, v in (extra_params or {}).items():
        if isinstance(v, str) and v.lower() in ("true", "false"):
            v = v.lower() == "true"
        clean[k] = v

    params = SubmissionParams(
        input_samplesheet=str(paths.samplesheet),
        outdir=str(paths.dir / "results"),
        genome=genome,
        profile=profile,
        extra_args=clean,
    )

    paths.params_file.write_text(json.dumps(params.to_dict(), indent=2) + "\n")

    return {
        "command": params.to_command_string(),
        "params_file": str(paths.params_file),
        "genome": genome,
        "profile": profile,
        "extra_params": extra_params or {},
    }

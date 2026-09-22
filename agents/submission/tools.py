"""Submission agent tools — configure nextflow parameters for nf-core/rnaseq.

The LLM reads the user's prompt and freely decides what params to set. The tool
accepts an open-ended dict so nothing is restricted by the Python side.
"""

from __future__ import annotations

import json
import platform
import shutil
import subprocess
from typing import Any

from agents.submission.params import SubmissionParams
from core.session import SESSION


def check_resources() -> dict[str, Any]:
    """Check system resources (RAM, CPUs, disk, Docker) for pipeline configuration."""
    info: dict[str, Any] = {"platform": platform.system()}

    info["cpus"] = _cpu_count()
    info["memory_gb"] = _memory_gb()
    info["disk"] = _disk_info()
    info["docker"] = _docker_info()

    return info


def _cpu_count() -> int | None:
    import os
    return os.cpu_count()


def _memory_gb() -> float | None:
    system = platform.system()
    try:
        if system == "Darwin":
            out = subprocess.run(
                ["sysctl", "-n", "hw.memsize"],
                capture_output=True, text=True, timeout=10,
            )
            if out.returncode == 0:
                return round(int(out.stdout.strip()) / (1024 ** 3), 1)
        elif system == "Linux":
            with open("/proc/meminfo") as f:
                for line in f:
                    if line.startswith("MemTotal:"):
                        kb = int(line.split()[1])
                        return round(kb / (1024 ** 2), 1)
    except (OSError, ValueError):
        pass
    return None


def _disk_info() -> dict[str, Any]:
    try:
        usage = shutil.disk_usage(".")
        return {
            "total_gb": round(usage.total / (1024 ** 3), 1),
            "free_gb": round(usage.free / (1024 ** 3), 1),
        }
    except OSError:
        return {"error": "unable to check disk"}


def _docker_info() -> dict[str, Any]:
    try:
        result = subprocess.run(
            ["docker", "info", "--format",
             "{{.MemTotal}}\n{{.NCPU}}\n{{.ServerVersion}}"],
            capture_output=True, text=True, timeout=15,
        )
        if result.returncode != 0:
            return {"available": False, "error": result.stderr.strip()[:200]}
        lines = result.stdout.strip().splitlines()
        return {
            "available": True,
            "memory_bytes": int(lines[0]) if lines[0].isdigit() else lines[0],
            "cpus": int(lines[1]) if len(lines) > 1 and lines[1].isdigit() else None,
            "version": lines[2] if len(lines) > 2 else None,
        }
    except FileNotFoundError:
        return {"available": False, "error": "docker not found on PATH"}
    except (OSError, subprocess.TimeoutExpired):
        return {"available": False, "error": "docker info timed out"}


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
        "params_file": str(paths.params_file),
        "genome": genome,
        "profile": profile,
        "extra_params": extra_params or {},
    }

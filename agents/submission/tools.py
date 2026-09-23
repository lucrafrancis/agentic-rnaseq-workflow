"""Submission agent tools — configure nextflow parameters for nf-core/rnaseq.

The LLM reads the user's prompt and freely decides what params to set. The tool
accepts an open-ended dict so nothing is restricted by the Python side.
"""

from __future__ import annotations

import json
import platform
import re
import shutil
import subprocess
from typing import Any

from agents.submission.params import SubmissionParams
from agents.submission.reference import find_reference, salmon_only
from core.config import NFCORE_REVISION
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
            "memory_gb": round(int(lines[0]) / (1024 ** 3), 1) if lines[0].isdigit() else None,
            "cpus": int(lines[1]) if len(lines) > 1 and lines[1].isdigit() else None,
            "version": lines[2] if len(lines) > 2 else None,
        }
    except FileNotFoundError:
        return {"available": False, "error": "docker not found on PATH"}
    except (OSError, subprocess.TimeoutExpired):
        return {"available": False, "error": "docker info timed out"}


_MEMORY_UNITS_GB = {"B": 1 / 1024 ** 3, "KB": 1 / 1024 ** 2, "MB": 1 / 1024, "GB": 1, "TB": 1024}


def _limit_gb(value: Any) -> float | None:
    """'32.GB', '32 GB', '32GB', '500.MB' -> GB (Nextflow units are binary). None if unparseable."""
    m = re.fullmatch(r"\s*(\d+(?:\.\d+)?)\s*\.?\s*([KMGT]?B)\s*", str(value), re.IGNORECASE)
    return float(m.group(1)) * _MEMORY_UNITS_GB[m.group(2).upper()] if m else None


def resource_limit_error(profile: str, params: dict[str, Any]) -> str | None:
    """Why max_memory can't work on this machine, or None if it can (or isn't set).

    A limit above the host's RAM, or above Docker's memory with the docker profile,
    gets tasks killed mid-run, so it is rejected rather than trusted.
    """
    if "max_memory" not in params:
        return None
    requested = _limit_gb(params["max_memory"])
    if requested is None:
        return f"max_memory '{params['max_memory']}' is not a memory size like '22.GB'."
    host = _memory_gb()
    if host and requested > host:
        return f"max_memory {requested:g} GB exceeds the machine's {host} GB of RAM."
    if "docker" in profile.split(","):
        docker = _docker_info()
        if docker.get("memory_gb") and requested > docker["memory_gb"]:
            return (f"max_memory {requested:g} GB exceeds the {docker['memory_gb']} GB available to Docker; "
                    "containers are killed above that. Set it below Docker's memory.")
    return None


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

    # skip_alignment without a pseudo-aligner produces no counts — useless.
    if clean.get("skip_alignment") and "pseudo_aligner" not in clean:
        clean["pseudo_aligner"] = "salmon"

    error = resource_limit_error(profile, clean)
    if error:
        return {"error": "resource_limit", "message": error}

    # Salmon-only runs reuse a cached index for this genome/revision, or build and save one.
    # Skipped when the LLM set reference files itself.
    ref_dir = None
    if salmon_only(clean) and not {"fasta", "gtf", "transcript_fasta", "salmon_index"} & clean.keys():
        ref_dir = find_reference(genome, NFCORE_REVISION)
        if ref_dir is None:
            clean.setdefault("save_reference", True)

    params = SubmissionParams(
        input_samplesheet=str(paths.samplesheet),
        outdir=str(paths.dir / "results"),
        genome=genome,
        profile=profile,
        extra_args=clean,
        reference=str(ref_dir) if ref_dir else None,
    )

    paths.params_file.write_text(json.dumps(params.to_dict(), indent=2) + "\n")

    return {
        "params_file": str(paths.params_file),
        "genome": genome,
        "profile": profile,
        "extra_params": clean,
        "reference": (f"reusing cached {genome} reference and Salmon index" if ref_dir
                      else "building the reference; saved for reuse after a successful run"
                      if clean.get("save_reference") else "nf-core default (--genome)"),
    }

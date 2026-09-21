"""Build and submit nf-core/rnaseq nextflow commands.

Two modes:
  - Default: sensible nf-core/rnaseq defaults, human approves the command before launch.
  - Prompt-driven: user provides instructions, the agent parses them and may make
    additional decisions (e.g. matching genome references). Human approves before launch
    in this mode too.
"""

from __future__ import annotations

import json
import subprocess
import threading
import time
from pathlib import Path

from agents.submission.errors import diagnose_log, format_diagnoses
from agents.submission.params import SubmissionParams, default_params
from core.approval import present_for_approval
from core.session import SESSION


def build_submission(
    samplesheet_path: str,
    outdir: str | None = None,
    params_overrides: dict | None = None,
) -> SubmissionParams:
    """Construct submission parameters from a sample sheet and optional overrides.

    Parameters
    ----------
    samplesheet_path:
        Path to the validated sample sheet CSV.
    outdir:
        Output directory for nf-core results. Defaults to a subdirectory of the run dir.
    params_overrides:
        Optional dict of parameter overrides (genome, profile, aligner, etc.).

    Returns
    -------
    SubmissionParams ready for approval and launch.
    """
    paths = SESSION.require_paths()
    if outdir is None:
        outdir = str(paths.dir / "results")

    params = default_params(samplesheet_path, outdir)
    if params_overrides:
        for key, value in params_overrides.items():
            if hasattr(params, key):
                setattr(params, key, value)
            else:
                params.extra_args[key] = value

    params_file = paths.params_file
    params_file.write_text(json.dumps(params.to_dict(), indent=2) + "\n")

    return params


def submit_and_monitor(params: SubmissionParams) -> dict:
    """Present the nextflow command for approval, then submit and monitor.

    Returns a summary dict with the outcome: success/failure, log path, and any
    diagnosed errors.
    """
    paths = SESSION.require_paths()
    command = params.to_command_string()

    result = present_for_approval(
        title="nf-core/rnaseq submission",
        preview=command,
        summary_stats=params.to_dict(),
    )
    if not result.approved:
        return {"submitted": False, "reason": result.reason or "User rejected submission."}

    log_path = paths.nextflow_log
    try:
        proc = subprocess.Popen(
            params.to_nextflow_args(),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            cwd=str(paths.dir),
        )

        _HANG_TIMEOUT = 60
        pipeline_failed = False

        def _watchdog():
            """Kill nextflow if it hangs after a detected failure."""
            time.sleep(_HANG_TIMEOUT)
            if proc.poll() is None:
                print(f"\nNextflow hung for {_HANG_TIMEOUT}s after failure. Terminating...")
                proc.terminate()

        with log_path.open("w") as log_file:
            for line in proc.stdout:
                print(line, end="")
                log_file.write(line)
                log_file.flush()
                if not pipeline_failed and (
                    "Pipeline completed with errors" in line
                    or "Pipeline failed" in line
                ):
                    pipeline_failed = True
                    threading.Thread(target=_watchdog, daemon=True).start()
        proc.wait()
        print(f"\nNextflow exited (code {proc.returncode}).")

        if proc.returncode == 0:
            SESSION.mark_stage_complete("submission")
            return {
                "submitted": True,
                "success": True,
                "log_path": str(log_path),
                "outdir": params.outdir,
            }

        diagnoses = diagnose_log(str(log_path))
        return {
            "submitted": True,
            "success": False,
            "returncode": proc.returncode,
            "log_path": str(log_path),
            "diagnosis": format_diagnoses(diagnoses),
        }

    except FileNotFoundError:
        return {
            "submitted": False,
            "error": "nextflow_not_found",
            "message": "Nextflow is not installed or not on PATH.",
        }

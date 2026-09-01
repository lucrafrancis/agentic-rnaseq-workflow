"""Entrypoint: uv run python run.py prompt.txt

Reads a human-written prompt from a text file and hands it to the samplesheet agent.
The prompt file is copied into the run directory as part of the audit trail.

After the samplesheet agent finishes, the human reviews the sample sheet and approves,
edits, or rejects it. On approval, the pipeline submission is built and presented for
a second approval before launching nextflow.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

from core.approval import present_for_approval
from core.session import SESSION


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit("usage: python run.py <prompt.txt>")

    prompt_file = Path(sys.argv[1])
    if not prompt_file.is_file():
        sys.exit(f"File not found: {prompt_file}")

    prompt = prompt_file.read_text().strip()
    if not prompt:
        sys.exit("Prompt file is empty.")

    SESSION.begin_run("rnaseq")
    shutil.copy2(prompt_file, SESSION.paths.dir / "prompt.txt")
    print(f"Run directory: {SESSION.paths.dir}")

    # --- Stage 1: samplesheet generation ---
    from agents.samplesheet.loop import run_samplesheet_agent
    run_samplesheet_agent(prompt)

    paths = SESSION.require_paths()
    samplesheet = paths.samplesheet
    report = paths.dir / "report.md"

    if not samplesheet.is_file():
        sys.exit("Samplesheet agent did not produce a sample sheet. Check the logs.")

    preview = samplesheet.read_text()
    if report.is_file():
        print(f"\nReport: {report}")

    result = present_for_approval(
        title="Sample sheet",
        preview=preview,
        file_path=samplesheet,
        summary_stats={"samples": preview.count("\n") - 1},
    )

    if not result.approved:
        print(f"Sample sheet rejected. Reason: {result.reason or 'none given'}")
        sys.exit(1)

    if result.edited:
        preview = samplesheet.read_text()
        print(f"Using edited sample sheet ({preview.count(chr(10)) - 1} samples).")

    # --- Stage 2: pipeline submission ---
    from agents.submission.submit import build_submission, submit_and_monitor

    params = build_submission(samplesheet_path=str(samplesheet))
    outcome = submit_and_monitor(params)

    if not outcome.get("submitted"):
        print(f"\nSubmission skipped: {outcome.get('reason') or outcome.get('message')}")
        sys.exit(1)

    if outcome.get("success"):
        print(f"\nPipeline finished successfully.")
        print(f"  Results: {outcome['outdir']}")
    else:
        print(f"\nPipeline failed (exit code {outcome.get('returncode')}).")
        print(f"  Log: {outcome.get('log_path')}")
        if outcome.get("diagnosis"):
            print(f"  Diagnosis:\n{outcome['diagnosis']}")
        sys.exit(1)


if __name__ == "__main__":
    main()

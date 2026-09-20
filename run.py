"""Entrypoint:
  uv run python run.py <prompt.txt>         — full pipeline (stages 1-3)
  uv run python run.py --analyze <results>  — analysis only (stage 3)

Three stages with human approval between them:
  1. Samplesheet agent — scans FASTQs, builds sample sheet
  2. Submission handler — builds nextflow command, submits pipeline
  3. Analysis agent — downstream DE, enrichment, and reporting
"""

from __future__ import annotations

import base64
import shutil
import sys
from pathlib import Path

from core.approval import present_for_approval
from core.session import SESSION


def main() -> None:
    if len(sys.argv) == 3 and sys.argv[1] == "--analyze":
        results_dir = Path(sys.argv[2])
        if not results_dir.is_dir():
            sys.exit(f"Not a directory: {results_dir}")
        SESSION.begin_run("analysis")
        print(f"Run directory: {SESSION.paths.dir}")
        _run_analysis(str(results_dir))
        return

    if len(sys.argv) != 2:
        sys.exit("usage: python run.py <prompt.txt>\n       python run.py --analyze <results_dir>")

    prompt_file = Path(sys.argv[1])
    if not prompt_file.is_file():
        sys.exit(f"File not found: {prompt_file}")

    prompt = prompt_file.read_text().strip()
    if not prompt:
        sys.exit("Prompt file is empty.")

    project_name = prompt_file.resolve().parent.name.replace(" ", "_")
    SESSION.begin_run(project_name)
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

    # --- Stage 3: downstream analysis ---
    _run_analysis(outcome["outdir"])


def _run_analysis(results_dir: str) -> None:
    """Prompt for analysis instructions and run the analysis agent."""
    paths = SESSION.require_paths()

    print("\n" + "=" * 60)
    print("  Pipeline completed. Ready for downstream analysis.")
    print("=" * 60)

    try:
        analysis_input = input("\nAnalysis prompt (or 'skip'): ").strip()
    except EOFError:
        print("\nNo terminal input. Skipping analysis.")
        return

    if not analysis_input or analysis_input.lower() == "skip":
        print("Analysis skipped.")
        return

    # Optional paper PDF
    pdf_path = None
    try:
        pdf_input = input("Paper PDF path (optional, Enter to skip): ").strip()
        if pdf_input:
            pdf_path = Path(pdf_input)
            if not pdf_path.is_file():
                print(f"PDF not found: {pdf_path}. Continuing without it.")
                pdf_path = None
    except EOFError:
        pass

    # Build the user message
    context_parts = [analysis_input, f"\nThe nf-core/rnaseq results are at: {results_dir}"]
    if paths.design.is_file():
        context_parts.append(f"A design CSV is available at: {paths.design}")

    if pdf_path:
        pdf_b64 = base64.standard_b64encode(pdf_path.read_bytes()).decode("ascii")
        user_message = [
            {
                "type": "document",
                "source": {"type": "base64", "media_type": "application/pdf", "data": pdf_b64},
            },
            {"type": "text", "text": "\n".join(context_parts)},
        ]
    else:
        user_message = "\n".join(context_parts)

    from agents.analysis.loop import run_analysis_agent
    run_analysis_agent(user_message)
    SESSION.mark_stage_complete("analysis")

    if paths.analysis_report.is_file():
        print(f"\nAnalysis report: {paths.analysis_report}")
    print("Done.")


if __name__ == "__main__":
    main()

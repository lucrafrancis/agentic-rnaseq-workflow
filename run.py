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
import json
import re
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

    # --- Stage 0: data download (if GEO/SRA accession in prompt) ---
    accession = _extract_accession(prompt)
    if accession:
        _download_stage(prompt)

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

    # --- Stage 2: pipeline submission (with troubleshooting retry loop) ---
    outcome = _submit_with_troubleshooting(str(samplesheet))

    # --- Stage 3: downstream analysis ---
    _run_analysis(outcome["outdir"])


def _extract_accession(prompt: str) -> str | None:
    """Extract a GEO/SRA accession from the prompt text."""
    m = re.search(r"\b(GSE\d+|SRP\d+|ERP\d+|DRP\d+|PRJNA\d+)\b", prompt, re.IGNORECASE)
    return m.group(1) if m else None


def _download_stage(prompt: str) -> None:
    """Stage 0: resolve accession and download FASTQs if needed."""
    from agents.download.loop import run_download_agent
    from agents.download.tools import execute_download, validate_downloads

    print("\n--- Stage 0: Data Download ---")
    run_download_agent(prompt)

    paths = SESSION.require_paths()

    if not paths.download_script.is_file():
        print("No downloads needed — files already present.")
        return

    script_content = paths.download_script.read_text()
    metadata = json.loads(paths.download_metadata.read_text())

    n_downloads = metadata.get("n_downloads", "?")
    download_bytes = metadata.get("download_bytes", 0)
    size_str = f"{download_bytes / (1024**3):.2f} GB" if download_bytes else "unknown"

    result = present_for_approval(
        title="Download script",
        preview=script_content,
        file_path=paths.download_script,
        summary_stats={"files_to_download": n_downloads, "estimated_size": size_str},
    )

    if not result.approved:
        print(f"Download rejected: {result.reason or 'none given'}")
        sys.exit(1)

    print("\nDownloading...")
    exec_result = execute_download(str(paths.download_script))
    if not exec_result["success"]:
        sys.exit("Download failed. Check output above.")

    print("Validating checksums...")
    val_result = validate_downloads(metadata["runs"], metadata["output_dir"])
    if not val_result["all_valid"]:
        print("\n⚠️  Checksum validation failed:")
        for f in val_result["files"]:
            if f["status"] not in ("pass", "no_checksum"):
                print(f"  {f['filename']}: {f['status']}")
        sys.exit(1)

    print(f"✓ All {val_result['n_pass']} file(s) validated (MD5).")


def _submit_with_troubleshooting(samplesheet_path: str) -> dict:
    """Run Stage 2 with LLM-assisted troubleshooting on failure."""
    from agents.submission.submit import build_submission, submit_and_monitor
    from agents.submission.troubleshoot import (
        MAX_RETRIES,
        apply_parameter_fix,
        diagnose_and_propose,
        wait_for_user_action,
    )

    params = build_submission(samplesheet_path=samplesheet_path)
    history: list[dict] = []

    for attempt in range(1, MAX_RETRIES + 1):
        outcome = submit_and_monitor(params)

        if not outcome.get("submitted"):
            print(f"\nSubmission skipped: {outcome.get('reason') or outcome.get('message')}")
            sys.exit(1)

        if outcome.get("success"):
            print(f"\nPipeline finished successfully.")
            print(f"  Results: {outcome['outdir']}")

            from agents.submission.troubleshoot import review_warnings
            print("Checking for warnings...")
            warnings = review_warnings(outcome.get("log_path", ""))
            if warnings:
                print(f"\n⚠ Pipeline warnings:\n{warnings}")
            else:
                print("No warnings.")

            return outcome

        print(f"\nPipeline failed (exit code {outcome.get('returncode')}).")
        print(f"  Log: {outcome.get('log_path')}")

        if attempt == MAX_RETRIES:
            print(f"\nExhausted {MAX_RETRIES} attempts. Manual intervention needed.")
            if history:
                print("\nTroubleshooting history:")
                for h in history:
                    print(f"  Attempt {h['attempt']}: {h['root_cause']} → {h['fix_description']}")
            sys.exit(1)

        print("Analysing failure...")
        try:
            proposal = diagnose_and_propose(outcome, params, attempt, history)
        except Exception as exc:
            print(f"Troubleshooting failed: {exc}")
            sys.exit(1)

        if proposal is None:
            print("Pipeline aborted by user.")
            sys.exit(1)

        category = proposal.get("category", "unfixable")
        if category == "unfixable":
            print("Issue is beyond automated troubleshooting. Manual intervention needed.")
            sys.exit(1)

        if category == "parameter_change" and proposal.get("parameter_changes"):
            apply_parameter_fix(params, proposal["parameter_changes"])
            print("Parameters updated. Retrying...")
        elif category in ("user_action", "config_change"):
            if not wait_for_user_action():
                print("Pipeline aborted by user.")
                sys.exit(1)
            print("Retrying...")

        history.append({
            "attempt": attempt,
            "root_cause": proposal.get("root_cause", "unknown"),
            "fix_description": proposal.get("fix_description", "none"),
            "category": category,
            "outcome": "retrying",
        })

    sys.exit(1)


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

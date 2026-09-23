"""Entrypoint:
  uv run python run.py <prompt.txt>           — full pipeline (stages 0-3)
  uv run python run.py --resume <run_dir>     — resume from an existing run directory
  uv run python run.py --analyze <run_dir>    — re-run analysis (stage 3) on an existing run
  uv run python run.py --analyze <prompt.txt> — analysis only, from a count matrix in the prompt
  --skip-download                             — skip stage 0 (FASTQs already local)

Stages with human approval between them:
  0. Download agent — resolve accessions, download FASTQs (if needed)
  1. Samplesheet agent — scans FASTQs, builds sample sheet
  2. Submission handler — builds nextflow command, submits pipeline
  3. Analysis agent — downstream DE, enrichment, and reporting

--resume skips stages recorded as complete (approved/validated) in run_state.json.
"""

from __future__ import annotations

import argparse
import base64
import csv
import json
import re
import shutil
import sys
from pathlib import Path

from core.approval import present_for_approval
from core.session import SESSION


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Agentic nf-core/rnaseq workflow.",
        epilog=(
            "examples:\n"
            "  run.py prompt.txt                  full pipeline\n"
            "  run.py prompt.txt --skip-download  FASTQs already local; ignore accessions\n"
            "  run.py --resume runs/<run_dir>     continue an interrupted run\n"
            "  run.py --analyze runs/<run_dir>    re-run analysis on an existing run\n"
            "  run.py --analyze prompt.txt        analysis only, from a count matrix in the prompt"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("target", type=Path, help="prompt file, or a run directory with --resume/--analyze")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--resume", action="store_true", help="resume the run directory TARGET")
    mode.add_argument("--analyze", action="store_true", help="analysis only (TARGET: run directory or prompt file)")
    parser.add_argument("--skip-download", action="store_true", help="skip stage 0 even if the prompt has an accession")
    args = parser.parse_args()
    target: Path = args.target

    if args.analyze and target.is_dir():
        SESSION.resume_run(target)
        print(f"Resuming run for analysis: {target}")
        prompt_file = target / "prompt.txt"
        original_prompt = prompt_file.read_text().strip() if prompt_file.is_file() else ""
        if not original_prompt:
            print("Warning: no prompt.txt in run directory.")
        if SESSION.mode == "analysis":
            _run_analysis_only(original_prompt, skip_download=args.skip_download)
        else:
            results_dir = target / "results"
            _run_analysis(str(results_dir) if results_dir.is_dir() else None, original_prompt=original_prompt)
        return

    if args.resume:
        if not target.is_dir():
            sys.exit(f"Not a directory: {target}")
        SESSION.resume_run(target)
        prompt_file = target / "prompt.txt"
        if not prompt_file.is_file():
            sys.exit(f"No prompt.txt in {target}. Cannot resume.")
        print(f"Resuming run: {target}")
        prompt = prompt_file.read_text().strip()
        if SESSION.mode == "analysis":
            _run_analysis_only(prompt, skip_download=args.skip_download)
        else:
            _run_pipeline(prompt, skip_download=args.skip_download)
        return

    if not target.is_file():
        sys.exit(f"File not found: {target}")
    prompt = target.read_text().strip()
    if not prompt:
        sys.exit("Prompt file is empty.")

    project_name = target.resolve().parent.name.replace(" ", "_")
    SESSION.begin_run(project_name, source_prompt=target, mode="analysis" if args.analyze else "pipeline")
    shutil.copy2(target, SESSION.paths.dir / "prompt.txt")
    print(f"Run directory: {SESSION.paths.dir}")

    if args.analyze:
        _run_analysis_only(prompt, skip_download=args.skip_download)
    else:
        _run_pipeline(prompt, skip_download=args.skip_download)


def _run_analysis_only(prompt: str, *, skip_download: bool = False) -> None:
    """Analysis-only runs: fetch a GEO count matrix if the prompt has an accession, then
    analyse. Never touches the samplesheet/nextflow stages."""
    if SESSION.is_complete("download"):
        print("Stage 0: count matrix already approved, skipping.")
    elif skip_download or not _extract_accession(prompt):
        print("Stage 0: using the local count matrix described in the prompt.")
        SESSION.mark_stage_complete("download")
    else:
        _counts_stage(prompt)
    _run_analysis(None, original_prompt=prompt)


def _counts_stage(prompt: str) -> None:
    """Stage 0 (counts mode): fetch a processed count matrix + design from GEO, then approve."""
    paths = SESSION.require_paths()
    print("\n--- Stage 0: GEO Count Matrix ---")
    if paths.counts_metadata.is_file() and paths.design.is_file():
        print("Found an unapproved count matrix, presenting for approval.")
    else:
        from agents.download.loop import run_counts_agent
        run_counts_agent(prompt)
        if not (paths.counts_metadata.is_file() and paths.design.is_file()):
            sys.exit("Counts agent did not produce a count matrix and design. Check the logs above.")

    meta = json.loads(paths.counts_metadata.read_text())
    design = {row["sample"]: row for row in csv.DictReader(paths.design.open())}
    extra = [c for c in next(iter(design.values())) if c not in ("sample", "condition", "gsm", "title")]
    covariates = [c for c in extra if len({row[c] for row in design.values()}) > 1]
    constants = {c: next(iter(design.values()))[c] for c in extra if c not in covariates}
    lines = [
        f"Source: {meta['filename']} ({meta['source']})",
        f"Values: {meta['value_type']}  |  genes: {meta['n_genes']}  |  gene IDs: {meta['gene_id_type']}"
        + (f"  |  duplicate IDs summed: {meta['n_duplicates_summed']}" if meta["n_duplicates_summed"] else ""),
    ]
    if meta["value_type"] != "raw_integer_counts":
        lines.append("WARNING: values are not raw integer counts — DESeq2 may not be appropriate.")
    lines += ["", f"{'file column':<28} {'GSM':<12} {'sample':<28} condition" + "".join(f" | {c}" for c in covariates)]
    for s in meta["samples"]:
        row = design.get(s["sample"], {})
        lines.append(
            f"{s['file_column']:<28} {s['gsm']:<12} {s['sample']:<28} {row.get('condition', '?')}"
            + "".join(f" | {row.get(c, '')}" for c in covariates)
        )
    if constants:
        lines.append("\nSame for all samples: " + "; ".join(f"{k} = {v}" for k, v in constants.items()))
    if meta["dropped_columns"]:
        lines.append(f"\nDropped file columns: {', '.join(meta['dropped_columns'][:15])}")
    if meta["unmapped_gsms"]:
        lines.append(f"GSMs not included: {', '.join(meta['unmapped_gsms'][:15])}")

    result = present_for_approval(
        title="GEO count matrix and design",
        preview="\n".join(lines),
        file_path=paths.design,
        summary_stats={"samples": len(meta["samples"]), "counts": str(paths.counts_matrix)},
    )
    if not result.approved:
        print(f"Count matrix rejected. Reason: {result.reason or 'none given'}")
        sys.exit("Adjust the prompt (e.g. which file or samples to use) and start a new run.")
    SESSION.mark_stage_complete("download")


def _run_pipeline(prompt: str, *, skip_download: bool = False) -> None:
    """Run the full pipeline, skipping stages already recorded as complete."""
    paths = SESSION.require_paths()

    # --- Stage 0: data download (if accession present) ---
    if SESSION.is_complete("download"):
        print("Stage 0: download already complete, skipping.")
    elif skip_download:
        print("Stage 0: skipped (--skip-download).")
        SESSION.mark_stage_complete("download")
    elif _extract_accession(prompt):
        _download_stage(prompt)

    # --- Stage 1: samplesheet generation ---
    samplesheet = paths.samplesheet
    if SESSION.is_complete("samplesheet"):
        print(f"Stage 1: samplesheet already approved ({samplesheet}), skipping.")
    else:
        if samplesheet.is_file():
            print(f"Stage 1: found unapproved samplesheet ({samplesheet}), presenting for approval.")
        else:
            from agents.samplesheet.loop import run_samplesheet_agent
            run_samplesheet_agent(prompt + _sample_metadata_note())

            if not samplesheet.is_file():
                sys.exit("Samplesheet agent did not produce a sample sheet. Check the logs.")

            report = paths.dir / "report.md"
            if report.is_file():
                print(f"\nReport: {report}")

        from agents.samplesheet.tools import sample_label_table

        preview = samplesheet.read_text()
        labels = sample_label_table(paths)
        result = present_for_approval(
            title="Sample sheet",
            preview=f"{labels}\n\n{preview}" if labels else preview,
            file_path=samplesheet,
            summary_stats={"samples": len(preview.strip().splitlines()) - 1},
        )

        if not result.approved:
            print(f"Sample sheet rejected. Reason: {result.reason or 'none given'}")
            sys.exit(1)

        if result.edited:
            preview = samplesheet.read_text()
            print(f"Using edited sample sheet ({len(preview.strip().splitlines()) - 1} samples).")

        SESSION.mark_stage_complete("samplesheet")

    # --- Stage 2: submission configuration + pipeline run ---
    results_dir = paths.dir / "results"
    quant_complete = False
    if results_dir.is_dir():
        for subdir in ("star_salmon", "salmon"):
            quant_path = results_dir / subdir / "salmon.merged.gene_counts.tsv"
            if quant_path.is_file():
                quant_complete = True
                break
    if quant_complete:
        print(f"Stage 2: results exist ({results_dir}), skipping submission.")
        outcome = {"outdir": str(results_dir)}
    else:
        if not paths.params_file.is_file():
            from agents.submission.loop import run_submission_agent

            samplesheet_content = samplesheet.read_text()
            n_samples = samplesheet_content.count("\n") - 1
            submission_prompt = (
                f"{prompt}\n\n"
                f"Approved samplesheet ({n_samples} samples): {samplesheet}\n"
                f"Preview:\n{samplesheet_content}"
            )
            print("\n--- Stage 2: Submission Configuration ---")
            run_submission_agent(submission_prompt)
        else:
            print(f"Stage 2: params.json exists, skipping submission agent.")

        outcome = _submit_with_troubleshooting(str(samplesheet))

    # --- Stage 3: downstream analysis ---
    _run_analysis(outcome["outdir"], original_prompt=prompt)


def _sample_metadata_note() -> str:
    """Write the downloaded runs' sample labels to sample_metadata.csv and return a prompt
    note pointing the samplesheet agent at it. Labels come from download_metadata.json
    (GEO/ENA), never from the LLM. Empty string if nothing was downloaded."""
    from agents.download.tools import selected_runs

    paths = SESSION.require_paths()
    if not paths.download_metadata.is_file():
        return ""
    metadata = json.loads(paths.download_metadata.read_text())
    with paths.sample_metadata.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["run_accession", "sample", "title", "library_layout"])
        for r in selected_runs(metadata):
            writer.writerow([r["run_accession"], r.get("sample_alias", ""), r.get("sample_title", ""),
                             r.get("library_layout", "")])
    note = (
        f"\n\nSample metadata for the downloaded runs (run accession -> GEO sample -> title), "
        f"from GEO/ENA: {paths.sample_metadata}"
    )
    if metadata.get("output_dir"):
        note += f"\nThe FASTQs were downloaded to: {metadata['output_dir']}"
    return note


def _extract_accession(prompt: str) -> str | None:
    """Extract a GEO/SRA accession from the prompt text."""
    m = re.search(r"\b(GSE\d+|SRP\d+|ERP\d+|DRP\d+|PRJNA\d+)\b", prompt, re.IGNORECASE)
    return m.group(1) if m else None


def _download_stage(prompt: str) -> None:
    """Stage 0: resolve accession and download FASTQs if needed.

    If a previous attempt got as far as choosing an output directory, the script is
    regenerated directly from download_metadata.json (no LLM) — files already present
    with a valid MD5 are skipped, so an interrupted download picks up where it left off.
    """
    from agents.download.loop import run_download_agent
    from agents.download.tools import (
        execute_download,
        generate_download_script,
        selected_runs,
        unlabelled_runs,
        validate_downloads,
    )

    print("\n--- Stage 0: Data Download ---")
    paths = SESSION.require_paths()
    previous = json.loads(paths.download_metadata.read_text()) if paths.download_metadata.is_file() else {}
    paths.download_script.unlink(missing_ok=True)  # never approve a stale script

    if previous.get("output_dir"):
        print(f"Resuming download into {previous['output_dir']} (skipping files with valid MD5).")
        result = generate_download_script(previous["output_dir"], runs=previous.get("selected_runs"))
        if "error" in result:
            sys.exit(f"Could not regenerate download script: {result['message']}")
    else:
        run_download_agent(prompt)

    if not paths.download_metadata.is_file():
        sys.exit("Download agent failed to resolve the accession. Check the logs above.")

    # Hard stop: never proceed with runs that can't be tied to a sample.
    unlabelled = unlabelled_runs(selected_runs(json.loads(paths.download_metadata.read_text())))
    if unlabelled:
        sys.exit(
            f"Runs with no sample ID/title: {', '.join(unlabelled[:20])}. They can't be assigned "
            "to conditions, so the workflow stops here. Use the GEO series accession if one "
            "exists, or provide FASTQs with a metadata file and run with --skip-download."
        )

    if not paths.download_script.is_file():
        print("No downloads needed — files already present.")
        SESSION.mark_stage_complete("download")
        return

    script_content = paths.download_script.read_text()
    metadata = json.loads(paths.download_metadata.read_text())

    n_downloads = metadata.get("n_downloads", "?")
    download_bytes = metadata.get("download_bytes", 0)
    size_str = f"{download_bytes / (1024**3):.2f} GB" if download_bytes else "unknown"

    chosen = selected_runs(metadata)
    stats = {"files_to_download": n_downloads, "estimated_size": size_str}
    if metadata.get("selected_runs") is not None:
        stats["runs_selected"] = f"{len(chosen)} of {len(metadata['runs'])}"
        stats["selection"] = "".join(
            f"\n    {r['run_accession']}  {r.get('sample_alias', '')}  {r.get('sample_title', '')}" for r in chosen
        )

    result = present_for_approval(
        title="Download script",
        preview=script_content,
        file_path=paths.download_script,
        summary_stats=stats,
    )

    if not result.approved:
        print(f"Download rejected. Reason: {result.reason or 'none given'}")
        sys.exit("If the FASTQs are already local, re-run with --skip-download.")

    print("\nDownloading...")
    exec_result = execute_download(str(paths.download_script))
    if not exec_result["success"]:
        sys.exit("Download failed. Check output above.")

    output_dir = metadata.get("output_dir")
    if not output_dir:
        print("No output directory recorded in metadata. Skipping validation.")
        SESSION.mark_stage_complete("download")
        return

    print("Validating checksums...")
    val_result = validate_downloads(chosen, output_dir)
    if not val_result["all_valid"]:
        print("\n⚠️  Checksum validation failed:")
        for f in val_result["files"]:
            if f["status"] not in ("pass", "no_checksum"):
                print(f"  {f['filename']}: {f['status']}")
        sys.exit(1)

    print(f"✓ All {val_result['n_pass']} file(s) validated (MD5).")
    SESSION.mark_stage_complete("download")


def _submit_with_troubleshooting(samplesheet_path: str) -> dict:
    """Run Stage 2 with LLM-assisted troubleshooting on failure."""
    from agents.submission.params import SubmissionParams
    from agents.submission.submit import submit_and_monitor
    from agents.submission.troubleshoot import (
        MAX_RETRIES,
        apply_parameter_fix,
        diagnose_and_propose,
        wait_for_user_action,
    )

    paths = SESSION.require_paths()
    if paths.params_file.is_file():
        params = SubmissionParams.load(paths.params_file)
    else:
        from agents.submission.submit import build_submission
        params = build_submission(samplesheet_path=samplesheet_path)
    history: list[dict] = []

    for attempt in range(1, MAX_RETRIES + 1):
        outcome = submit_and_monitor(params, resume=(attempt > 1))

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


def _run_analysis(results_dir: str | None, *, original_prompt: str = "") -> None:
    """Run the analysis agent. Uses the original prompt as base context;
    optionally accepts extra instructions interactively (saved if provided).

    results_dir is None for analysis-only runs, where the prompt points at the count
    matrix (and optionally a design file) instead of nf-core outputs.
    """
    paths = SESSION.require_paths()

    # Fall back to prompt.txt in the run directory if not passed explicitly
    if not original_prompt:
        prompt_file = paths.dir / "prompt.txt"
        if prompt_file.is_file():
            original_prompt = prompt_file.read_text().strip()

    print("\n" + "=" * 60)
    print("  Ready for downstream analysis." if results_dir is None else "  Pipeline completed. Ready for downstream analysis.")
    print("=" * 60)

    # Extra instructions (optional)
    analysis_input = ""
    try:
        analysis_input = input("\nAdditional analysis instructions (Enter to continue, 'skip' to skip): ").strip()
    except EOFError:
        pass

    if analysis_input.lower() == "skip":
        print("Analysis skipped.")
        return

    if analysis_input:
        (paths.dir / "analysis_prompt.txt").write_text(analysis_input + "\n")
        print(f"Saved analysis prompt to {paths.dir / 'analysis_prompt.txt'}")

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
    context_parts = []
    if original_prompt:
        context_parts.append(original_prompt)
    if analysis_input:
        context_parts.append(f"\nAdditional instructions: {analysis_input}")
    if results_dir:
        context_parts.append(f"\nThe nf-core/rnaseq results are at: {results_dir}")
    else:
        context_parts.append("\nNo nf-core results for this run — use the count matrix described above.")
    if SESSION.source_prompt:
        context_parts.append(
            f"The prompt was read from {SESSION.source_prompt}; relative paths in it may be "
            f"relative to that file's directory or to the current directory ({Path.cwd()})."
        )
    if paths.counts_matrix.is_file():
        context_parts.append(f"The approved count matrix (fetched from GEO) is at: {paths.counts_matrix}")
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
    from agents.analysis.replay import write_replay_script

    log = paths.analysis_tool_log
    log_start = len(log.read_text().splitlines()) if log.is_file() else 0
    run_analysis_agent(user_message)
    SESSION.mark_stage_complete("analysis")

    if paths.analysis_report.is_file():
        print(f"\nAnalysis report: {paths.analysis_report}")
    replay = write_replay_script(log, paths.analysis_dir / "replay.py", start_line=log_start)
    if replay:
        print(f"Replay script: {replay}")
    print("Done.")


if __name__ == "__main__":
    main()

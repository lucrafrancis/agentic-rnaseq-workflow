"""Entrypoint: uv run python run.py <fastq_dir> [--metadata <path>]

Orchestrates the three stages of the RNA-seq workflow:
  1. Sample sheet generation (agent-driven)
  2. nf-core/rnaseq submission
  3. Post-pipeline analysis (stub)
"""

from __future__ import annotations

import argparse
import sys

from core.session import SESSION


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Agentic nf-core/rnaseq workflow",
    )
    parser.add_argument("fastq_dir", help="Directory containing FASTQ files")
    parser.add_argument("--metadata", help="Path to metadata CSV/TSV file", default=None)
    parser.add_argument("--name", help="Project name for the run directory", default="rnaseq")
    args = parser.parse_args()

    SESSION.begin_run(args.name)
    print(f"Run directory: {SESSION.paths.dir}")

    # Stage 1: Sample sheet generation
    print("\n--- Stage 1: Sample sheet generation ---")
    from agents.samplesheet.loop import run_samplesheet_agent

    try:
        result = run_samplesheet_agent(args.fastq_dir, args.metadata)
    except NotImplementedError:
        print("Sample sheet agent loop not yet implemented.")
        print("Run the tools manually or implement agents/samplesheet/loop.py.")
        return

    if "error" in result:
        sys.exit(f"Sample sheet generation failed: {result['error']}")

    samplesheet_path = result["samplesheet_path"]
    SESSION.mark_stage_complete("samplesheet")
    print(f"Sample sheet: {samplesheet_path}")

    # Stage 2: nf-core/rnaseq submission
    print("\n--- Stage 2: nf-core/rnaseq submission ---")
    from agents.submission.submit import build_submission, submit_and_monitor

    params = build_submission(samplesheet_path)
    outcome = submit_and_monitor(params)

    if not outcome.get("success"):
        print(f"Submission outcome: {outcome}")
        return

    # Stage 3: Post-pipeline analysis (stub)
    print("\n--- Stage 3: Post-pipeline analysis ---")
    print("Not yet implemented.")

    print(f"\nDone. Run directory: {SESSION.paths.dir}")


if __name__ == "__main__":
    main()

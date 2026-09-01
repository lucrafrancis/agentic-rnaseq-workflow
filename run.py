"""Entrypoint: uv run python run.py prompt.txt

Reads a human-written prompt from a text file and hands it to the samplesheet agent.
The prompt file is copied into the run directory as part of the audit trail.

After the agent finishes, the human reviews the report and sample sheet before
proceeding to nf-core submission (not yet wired up).
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

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

    from agents.samplesheet.loop import run_samplesheet_agent
    run_samplesheet_agent(prompt)

    print(f"\n{'='*60}")
    print("  Stage 1 complete")
    print(f"{'='*60}")
    print(f"  Run directory: {SESSION.paths.dir}")
    print(f"  Report:        {SESSION.paths.dir / 'report.md'}")
    print(f"\n  Review the report and sample sheet before proceeding.")


if __name__ == "__main__":
    main()

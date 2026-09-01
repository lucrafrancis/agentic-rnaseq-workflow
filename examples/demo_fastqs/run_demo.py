"""End-to-end demo of the samplesheet tools against realistic Illumina filenames.

Run: uv run python examples/demo_fastqs/run_demo.py
"""

from __future__ import annotations

import json
from pathlib import Path

from agents.samplesheet.tools import (
    draft_samplesheet,
    match_pairs,
    read_metadata,
    scan_fastqs,
    stage_fastqs,
    validate_samplesheet,
)

HERE = Path(__file__).parent
FASTQ_DIR = HERE / "fastqs"
METADATA = HERE / "metadata.csv"
STAGING_DIR = HERE / "staged"


def pp(label: str, result: dict) -> None:
    print(f"\n{'='*60}")
    print(f"  {label}")
    print(f"{'='*60}")
    print(json.dumps(result, indent=2))


def main() -> None:
    scan = scan_fastqs(str(FASTQ_DIR))
    pp("scan_fastqs", scan)

    meta = read_metadata(str(METADATA))
    pp("read_metadata", meta)

    pairs = match_pairs(scan["files"])
    pp("match_pairs", pairs)

    staged = stage_fastqs(pairs["matched"], str(FASTQ_DIR), str(STAGING_DIR))
    pp("stage_fastqs", staged)

    meta_lookup = {row["sample"]: row for row in meta["rows"]}
    draft = draft_samplesheet(staged["staged_pairs"], metadata=meta_lookup)
    pp("draft_samplesheet", draft)

    valid = validate_samplesheet(draft["csv_content"])
    pp("validate_samplesheet", valid)

    if valid["valid"]:
        print("\n\nSample sheet is valid and ready for approval.")
    else:
        print("\n\nSample sheet has errors — see above.")


if __name__ == "__main__":
    main()

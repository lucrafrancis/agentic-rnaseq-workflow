"""Human-in-the-loop approval flow.

Presents a draft artifact (sample sheet, nextflow command) to the user and waits for
them to approve, edit, or reject it. The caller supplies the draft and a description of
what it is; this module handles the presentation and input loop.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass
class ApprovalResult:
    """What the user decided."""

    approved: bool
    edited: bool = False
    reason: str | None = None


def present_for_approval(
    title: str,
    preview: str,
    file_path: Path | None = None,
    summary_stats: dict | None = None,
) -> ApprovalResult:
    """Show a draft to the user and ask for confirmation.

    Parameters
    ----------
    title:
        What the draft is (e.g. "Sample sheet" or "Nextflow submission command").
    preview:
        The first N rows or the command string — enough to judge correctness.
    file_path:
        Path to the full artifact on disk, so the user can inspect it.
    summary_stats:
        Key numbers (total samples, pairs found, warnings) shown alongside the preview.

    Returns
    -------
    ApprovalResult with the user's decision.
    """
    print(f"\n{'='*60}")
    print(f"  {title}")
    print(f"{'='*60}\n")
    print(preview)

    if summary_stats:
        print("\nSummary:")
        for key, value in summary_stats.items():
            print(f"  {key}: {value}")

    if file_path:
        print(f"\nFull file: {file_path}")

    if not sys.stdin.isatty():
        print("\n[non-interactive: auto-approved]")
        return ApprovalResult(approved=True)

    print()
    while True:
        choice = input("[a]pprove / [e]dit / [r]eject: ").strip().lower()
        if choice in ("a", "approve"):
            return ApprovalResult(approved=True)
        if choice in ("e", "edit"):
            if file_path:
                print(f"Edit the file at {file_path}, then press Enter to continue.")
                input()
                return ApprovalResult(approved=True, edited=True)
            else:
                print("No file to edit. Approve or reject.")
                continue
        if choice in ("r", "reject"):
            reason = input("Reason (optional): ").strip() or None
            return ApprovalResult(approved=False, reason=reason)
        print("Please enter 'a', 'e', or 'r'.")

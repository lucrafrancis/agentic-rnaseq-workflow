"""The agent loop for sample sheet generation.

Stub — the actual Claude API loop will be built next. This module defines the interface
that run.py calls into; the implementation will follow the same pattern as the scrna
workflow's loop.py: send conversation + tool schemas to Claude, dispatch tool calls,
append results, repeat until the agent finishes or MAX_TURNS is hit.
"""

from __future__ import annotations


def run_samplesheet_agent(fastq_dir: str, metadata_path: str | None = None) -> dict:
    """Drive sample sheet generation to completion.

    Parameters
    ----------
    fastq_dir:
        Path to the directory containing FASTQ files.
    metadata_path:
        Optional path to a metadata CSV/TSV file.

    Returns
    -------
    dict with 'samplesheet_path' on success, or 'error' on failure.
    """
    raise NotImplementedError("Agent loop not yet implemented.")

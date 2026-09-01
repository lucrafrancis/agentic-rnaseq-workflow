"""Test fixtures.

Everything here runs offline: tiny empty FASTQs and metadata files stand in for real
data so the samplesheet tools can be exercised in milliseconds without downloads or
the API. Outputs are redirected to a temp dir and the SESSION singleton is reset
between tests.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from core import config
from core.session import SESSION


@pytest.fixture
def fastq_dir(tmp_path: Path) -> Path:
    """A directory with three paired-end samples (empty .fastq.gz files)."""
    d = tmp_path / "fastqs"
    d.mkdir()
    for sample in ("sampleA", "sampleB", "sampleC"):
        (d / f"{sample}_R1.fastq.gz").touch()
        (d / f"{sample}_R2.fastq.gz").touch()
    return d


@pytest.fixture
def single_end_dir(tmp_path: Path) -> Path:
    """A directory with single-end FASTQs (R1 only)."""
    d = tmp_path / "fastqs_se"
    d.mkdir()
    for sample in ("se_sample1", "se_sample2"):
        (d / f"{sample}_R1.fastq.gz").touch()
    return d


@pytest.fixture
def metadata_csv(tmp_path: Path) -> Path:
    """A simple metadata CSV mapping samples to strandedness."""
    path = tmp_path / "metadata.csv"
    path.write_text(
        "sample,strandedness,condition\n"
        "sampleA,reverse,treated\n"
        "sampleB,reverse,control\n"
        "sampleC,forward,treated\n"
    )
    return path


@pytest.fixture
def metadata_tsv(tmp_path: Path) -> Path:
    """A TSV metadata file."""
    path = tmp_path / "metadata.tsv"
    path.write_text(
        "sample\tstrandedness\n"
        "sampleA\treverse\n"
        "sampleB\tforward\n"
    )
    return path


@pytest.fixture(autouse=True)
def isolate(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Redirect outputs to a temp dir and reset the SESSION singleton for each test."""
    monkeypatch.setattr(config, "RUNS_DIR", tmp_path / "runs")
    SESSION.__init__()

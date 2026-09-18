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


@pytest.fixture
def count_matrix_tsv(tmp_path: Path) -> Path:
    """A small synthetic count matrix TSV matching nf-core output format."""
    import random
    random.seed(42)

    path = tmp_path / "salmon.merged.gene_counts.tsv"
    header = "gene_id\tgene_name\tWT_REP1\tWT_REP2\tKO_REP1\tKO_REP2"
    rows = [header]
    for i in range(20):
        gene_id = f"GENE{i:03d}"
        gene_name = f"Gene{i}"
        wt1 = random.randint(50, 500)
        wt2 = random.randint(50, 500)
        # First 5 genes upregulated in KO, next 5 downregulated, rest unchanged
        if i < 5:
            ko1 = int(wt1 * 3.0 + random.randint(-10, 10))
            ko2 = int(wt2 * 3.0 + random.randint(-10, 10))
        elif i < 10:
            ko1 = max(0, int(wt1 * 0.3 + random.randint(-10, 10)))
            ko2 = max(0, int(wt2 * 0.3 + random.randint(-10, 10)))
        else:
            ko1 = wt1 + random.randint(-20, 20)
            ko2 = wt2 + random.randint(-20, 20)
        rows.append(f"{gene_id}\t{gene_name}\t{wt1}\t{wt2}\t{ko1}\t{ko2}")
    path.write_text("\n".join(rows) + "\n")
    return path


@pytest.fixture
def design_csv(tmp_path: Path) -> Path:
    """A design CSV mapping samples to conditions."""
    path = tmp_path / "design.csv"
    path.write_text(
        "sample,condition\n"
        "WT_REP1,WT\n"
        "WT_REP2,WT\n"
        "KO_REP1,KO\n"
        "KO_REP2,KO\n"
    )
    return path


@pytest.fixture
def nfcore_results_dir(tmp_path: Path, count_matrix_tsv: Path) -> Path:
    """A minimal nf-core/rnaseq output directory structure."""
    import shutil

    results = tmp_path / "results"
    star_salmon = results / "star_salmon"
    star_salmon.mkdir(parents=True)
    shutil.copy2(count_matrix_tsv, star_salmon / "salmon.merged.gene_counts.tsv")

    # TPM file with same structure
    tpm_path = star_salmon / "salmon.merged.gene_tpm.tsv"
    tpm_path.write_text(count_matrix_tsv.read_text())

    # MultiQC stub
    mqc_dir = results / "multiqc" / "star_salmon" / "multiqc_report_data"
    mqc_dir.mkdir(parents=True)
    (mqc_dir / "multiqc_general_stats.txt").write_text(
        "Sample\tstar-mapped_percent\tsamtools_stats-reads_mapped_percent\n"
        "WT_REP1\t90.0\t95.0\n"
        "WT_REP2\t91.0\t96.0\n"
        "KO_REP1\t89.0\t94.0\n"
        "KO_REP2\t92.0\t97.0\n"
    )
    return results


@pytest.fixture(autouse=True)
def isolate(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Redirect outputs to a temp dir and reset the SESSION singleton for each test."""
    monkeypatch.setattr(config, "RUNS_DIR", tmp_path / "runs")
    SESSION.__init__()

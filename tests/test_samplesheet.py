"""Tests for the sample sheet tools — offline, no API calls."""

from __future__ import annotations

from pathlib import Path

from agents.samplesheet.tools import (
    draft_samplesheet,
    match_pairs,
    read_metadata,
    scan_fastqs,
    validate_samplesheet,
)


class TestScanFastqs:
    def test_finds_paired_files(self, fastq_dir: Path):
        result = scan_fastqs(str(fastq_dir))
        assert result["n_files"] == 6
        assert all(f.endswith(".fastq.gz") for f in result["files"])

    def test_returns_error_for_bad_dir(self):
        result = scan_fastqs("/nonexistent/path")
        assert "error" in result

    def test_empty_directory(self, tmp_path: Path):
        d = tmp_path / "empty"
        d.mkdir()
        result = scan_fastqs(str(d))
        assert result["n_files"] == 0
        assert result["files"] == []


class TestReadMetadata:
    def test_reads_csv(self, metadata_csv: Path):
        result = read_metadata(str(metadata_csv))
        assert "sample" in result["columns"]
        assert "strandedness" in result["columns"]
        assert result["n_rows_shown"] == 3

    def test_reads_tsv(self, metadata_tsv: Path):
        result = read_metadata(str(metadata_tsv))
        assert "sample" in result["columns"]
        assert result["n_rows_shown"] == 2

    def test_returns_error_for_missing_file(self):
        result = read_metadata("/nonexistent/file.csv")
        assert "error" in result


class TestMatchPairs:
    def test_matches_three_pairs(self):
        files = [
            "sampleA_R1.fastq.gz", "sampleA_R2.fastq.gz",
            "sampleB_R1.fastq.gz", "sampleB_R2.fastq.gz",
            "sampleC_R1.fastq.gz", "sampleC_R2.fastq.gz",
        ]
        result = match_pairs(files)
        assert result["n_matched_pairs"] == 3
        assert result["n_unpaired"] == 0
        assert result["n_incomplete"] == 0

    def test_detects_unpaired(self):
        files = ["sampleA_R1.fastq.gz", "sampleB_R1.fastq.gz", "sampleB_R2.fastq.gz"]
        result = match_pairs(files)
        assert result["n_matched_pairs"] == 1
        assert result["n_incomplete"] == 1

    def test_handles_underscore_numbering(self):
        files = ["sample_1.fastq.gz", "sample_2.fastq.gz"]
        result = match_pairs(files)
        assert result["n_matched_pairs"] == 1


class TestDraftSamplesheet:
    def test_basic_draft(self):
        matches = [
            {"sample": "sampleA", "fastq_1": "/data/sampleA_R1.fq.gz", "fastq_2": "/data/sampleA_R2.fq.gz"},
            {"sample": "sampleB", "fastq_1": "/data/sampleB_R1.fq.gz", "fastq_2": "/data/sampleB_R2.fq.gz"},
        ]
        result = draft_samplesheet(matches)
        assert result["n_samples"] == 2
        assert "sample,fastq_1,fastq_2,strandedness" in result["preview"]
        assert "auto" in result["csv_content"]

    def test_metadata_overrides_strandedness(self):
        matches = [
            {"sample": "sampleA", "fastq_1": "a_R1.fq.gz", "fastq_2": "a_R2.fq.gz"},
        ]
        metadata = {"sampleA": {"strandedness": "reverse"}}
        result = draft_samplesheet(matches, metadata)
        assert "reverse" in result["csv_content"]


class TestValidateSamplesheet:
    def test_valid_sheet(self):
        sheet = (
            "sample,fastq_1,fastq_2,strandedness\n"
            "s1,/a_R1.fq.gz,/a_R2.fq.gz,auto\n"
            "s2,/b_R1.fq.gz,/b_R2.fq.gz,reverse\n"
        )
        result = validate_samplesheet(sheet)
        assert result["valid"] is True
        assert result["n_samples"] == 2
        assert result["errors"] == []

    def test_missing_columns(self):
        sheet = "sample,fastq_1\ns1,/a_R1.fq.gz\n"
        result = validate_samplesheet(sheet)
        assert result["valid"] is False
        assert any("Missing columns" in e for e in result["errors"])

    def test_duplicate_samples_warned(self):
        sheet = (
            "sample,fastq_1,fastq_2,strandedness\n"
            "s1,/a_R1.fq.gz,/a_R2.fq.gz,auto\n"
            "s1,/b_R1.fq.gz,/b_R2.fq.gz,auto\n"
        )
        result = validate_samplesheet(sheet)
        assert len(result["warnings"]) > 0

    def test_invalid_strandedness(self):
        sheet = (
            "sample,fastq_1,fastq_2,strandedness\n"
            "s1,/a_R1.fq.gz,/a_R2.fq.gz,banana\n"
        )
        result = validate_samplesheet(sheet)
        assert result["valid"] is False

    def test_empty_sheet(self):
        result = validate_samplesheet("")
        assert result["valid"] is False

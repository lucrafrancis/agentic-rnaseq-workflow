"""Tests for the download tools — offline, no API calls."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

import pytest

from agents.download.tools import (
    _compute_md5,
    _parse_ena_run,
    _resolve_gse,
    check_existing_files,
    generate_download_script,
    resolve_accession,
    validate_downloads,
)
from core.session import SESSION


@pytest.fixture
def run_dir(tmp_path: Path) -> Path:
    """Initialize a SESSION run for download tests."""
    SESSION.begin_run("test_download")
    return SESSION.require_paths().dir


def _write_metadata(runs: list[dict], accession: str = "GSE000") -> None:
    paths = SESSION.require_paths()
    paths.download_metadata.write_text(
        json.dumps({"accession": accession, "runs": runs}, indent=2)
    )


# --- _parse_ena_run ---


class TestParseEnaRun:
    def test_single_end(self):
        row = {
            "run_accession": "SRR12345",
            "fastq_ftp": "ftp.sra.ebi.ac.uk/vol1/fastq/SRR123/045/SRR12345/SRR12345.fastq.gz",
            "fastq_md5": "abc123",
            "fastq_bytes": "1000000",
            "library_layout": "SINGLE",
        }
        result = _parse_ena_run(row)
        assert result["run_accession"] == "SRR12345"
        assert len(result["files"]) == 1
        assert result["files"][0]["filename"] == "SRR12345.fastq.gz"
        assert result["files"][0]["md5"] == "abc123"
        assert result["files"][0]["bytes"] == 1000000

    def test_paired_end(self):
        row = {
            "run_accession": "SRR12345",
            "fastq_ftp": (
                "ftp.sra.ebi.ac.uk/vol1/SRR12345_1.fastq.gz;"
                "ftp.sra.ebi.ac.uk/vol1/SRR12345_2.fastq.gz"
            ),
            "fastq_md5": "abc123;def456",
            "fastq_bytes": "1000000;1100000",
            "library_layout": "PAIRED",
        }
        result = _parse_ena_run(row)
        assert len(result["files"]) == 2
        assert result["files"][0]["filename"] == "SRR12345_1.fastq.gz"
        assert result["files"][1]["filename"] == "SRR12345_2.fastq.gz"
        assert result["total_bytes"] == 2100000

    def test_empty_ftp(self):
        row = {
            "run_accession": "SRR12345",
            "fastq_ftp": "",
            "fastq_md5": "",
            "fastq_bytes": "",
        }
        result = _parse_ena_run(row)
        assert result["files"] == []


# --- resolve_accession ---


class TestResolveAccession:
    def test_invalid_accession_type(self, run_dir: Path):
        result = resolve_accession("INVALID123")
        assert result["error"] == "unknown_accession_type"

    @patch("agents.download.tools._resolve_gse")
    def test_gse_no_runs(self, mock_resolve, run_dir: Path):
        mock_resolve.return_value = []
        result = resolve_accession("GSE000000")
        assert result["error"] == "no_runs_found"

    @patch("agents.download.tools._resolve_gse")
    def test_gse_success(self, mock_resolve, run_dir: Path):
        mock_resolve.return_value = [
            {
                "run_accession": "SRR001",
                "library_layout": "SINGLE",
                "files": [
                    {"filename": "SRR001.fastq.gz", "url": "ftp://test", "md5": "abc", "bytes": 1000}
                ],
                "total_bytes": 1000,
            },
        ]
        result = resolve_accession("GSE123456")
        assert result["n_runs"] == 1
        assert result["n_files"] == 1
        assert result["accession"] == "GSE123456"

        meta = json.loads(SESSION.require_paths().download_metadata.read_text())
        assert meta["accession"] == "GSE123456"
        assert len(meta["runs"]) == 1

    @patch("agents.download.tools._resolve_study")
    def test_srp_accession(self, mock_resolve, run_dir: Path):
        mock_resolve.return_value = [
            {
                "run_accession": "SRR002",
                "library_layout": "PAIRED",
                "files": [],
                "total_bytes": 0,
            },
        ]
        result = resolve_accession("SRP282091")
        assert result["n_runs"] == 1


# --- _resolve_gse ---


class TestResolveGse:
    @patch("agents.download.tools._http_get_json")
    @patch("agents.download.tools._http_get")
    def test_success(self, mock_get, mock_get_json):
        mock_get_json.side_effect = [
            {"esearchresult": {"idlist": ["12345"]}},
            [
                {
                    "run_accession": "SRR001",
                    "fastq_ftp": "ftp.sra.ebi.ac.uk/SRR001.fastq.gz",
                    "fastq_md5": "abc123",
                    "fastq_bytes": "1000",
                    "library_layout": "SINGLE",
                }
            ],
        ]
        mock_get.return_value = "Run,SRAStudy\nSRR001,SRP001\n"

        runs = _resolve_gse("GSE123456")
        assert len(runs) == 1
        assert runs[0]["run_accession"] == "SRR001"

    @patch("agents.download.tools._http_get_json")
    def test_no_uids(self, mock_get_json):
        mock_get_json.return_value = {"esearchresult": {"idlist": []}}
        assert _resolve_gse("GSE000000") == []


# --- check_existing_files ---


class TestCheckExistingFiles:
    def test_no_metadata(self, run_dir: Path, tmp_path: Path):
        result = check_existing_files(str(tmp_path))
        assert result["error"] == "no_metadata"

    def test_all_missing(self, run_dir: Path, tmp_path: Path):
        _write_metadata([
            {
                "run_accession": "SRR001",
                "files": [{"filename": "SRR001.fastq.gz", "md5": "abc", "bytes": 100}],
            }
        ])
        search = tmp_path / "empty"
        search.mkdir()
        result = check_existing_files(str(search))
        assert result["n_missing"] == 1
        assert result["all_valid"] is False

    def test_all_valid(self, run_dir: Path, tmp_path: Path):
        fq = tmp_path / "SRR001.fastq.gz"
        fq.write_bytes(b"test content")
        md5 = _compute_md5(fq)

        _write_metadata([
            {
                "run_accession": "SRR001",
                "files": [{"filename": "SRR001.fastq.gz", "md5": md5, "bytes": 12}],
            }
        ])
        result = check_existing_files(str(tmp_path))
        assert result["n_valid"] == 1
        assert result["all_valid"] is True

    def test_corrupted(self, run_dir: Path, tmp_path: Path):
        fq = tmp_path / "SRR001.fastq.gz"
        fq.write_bytes(b"wrong content")

        _write_metadata([
            {
                "run_accession": "SRR001",
                "files": [{"filename": "SRR001.fastq.gz", "md5": "wrong_md5", "bytes": 13}],
            }
        ])
        result = check_existing_files(str(tmp_path))
        assert result["n_corrupted"] == 1
        assert result["all_valid"] is False


# --- generate_download_script ---


class TestGenerateDownloadScript:
    def test_no_metadata(self, run_dir: Path):
        result = generate_download_script("/tmp/out")
        assert result["error"] == "no_metadata"

    def test_generates_script(self, run_dir: Path, tmp_path: Path):
        _write_metadata([
            {
                "run_accession": "SRR001",
                "files": [
                    {
                        "filename": "SRR001.fastq.gz",
                        "url": "ftp://ftp.sra.ebi.ac.uk/test/SRR001.fastq.gz",
                        "md5": "abc123",
                        "bytes": 1000,
                    }
                ],
            }
        ])
        result = generate_download_script(str(tmp_path / "output"))
        assert result["n_files"] == 1
        assert result["script_path"] is not None
        assert Path(result["script_path"]).is_file()

        script = Path(result["script_path"]).read_text()
        assert "wget" in script
        assert "SRR001.fastq.gz" in script
        assert "abc123" in script

    def test_skips_valid_files(self, run_dir: Path, tmp_path: Path):
        out = tmp_path / "output"
        out.mkdir()
        fq = out / "SRR001.fastq.gz"
        fq.write_bytes(b"test")
        md5 = _compute_md5(fq)

        _write_metadata([
            {
                "run_accession": "SRR001",
                "files": [
                    {"filename": "SRR001.fastq.gz", "url": "ftp://test", "md5": md5, "bytes": 4}
                ],
            }
        ])
        result = generate_download_script(str(out))
        assert result["n_files"] == 0
        assert result["n_skipped"] == 1
        assert result["script_path"] is None

    def test_partial_download(self, run_dir: Path, tmp_path: Path):
        """Downloads only missing files when some are already valid."""
        out = tmp_path / "output"
        out.mkdir()
        fq1 = out / "SRR001.fastq.gz"
        fq1.write_bytes(b"valid")
        md5_valid = _compute_md5(fq1)

        _write_metadata([
            {
                "run_accession": "SRR001",
                "files": [
                    {"filename": "SRR001.fastq.gz", "url": "ftp://test/1", "md5": md5_valid, "bytes": 5}
                ],
            },
            {
                "run_accession": "SRR002",
                "files": [
                    {"filename": "SRR002.fastq.gz", "url": "ftp://test/2", "md5": "xyz", "bytes": 2000}
                ],
            },
        ])
        result = generate_download_script(str(out))
        assert result["n_files"] == 1
        assert result["n_skipped"] == 1


# --- validate_downloads ---


class TestValidateDownloads:
    def test_all_pass(self, tmp_path: Path):
        fq = tmp_path / "SRR001.fastq.gz"
        fq.write_bytes(b"good data")
        md5 = _compute_md5(fq)

        runs = [
            {
                "run_accession": "SRR001",
                "files": [{"filename": "SRR001.fastq.gz", "md5": md5, "bytes": 9}],
            }
        ]
        result = validate_downloads(runs, str(tmp_path))
        assert result["all_valid"] is True
        assert result["n_pass"] == 1

    def test_missing_file(self, tmp_path: Path):
        runs = [
            {
                "run_accession": "SRR001",
                "files": [{"filename": "SRR001.fastq.gz", "md5": "abc", "bytes": 100}],
            }
        ]
        result = validate_downloads(runs, str(tmp_path))
        assert result["all_valid"] is False
        assert result["n_missing"] == 1

    def test_checksum_mismatch(self, tmp_path: Path):
        fq = tmp_path / "SRR001.fastq.gz"
        fq.write_bytes(b"bad data")

        runs = [
            {
                "run_accession": "SRR001",
                "files": [{"filename": "SRR001.fastq.gz", "md5": "definitely_wrong", "bytes": 8}],
            }
        ]
        result = validate_downloads(runs, str(tmp_path))
        assert result["all_valid"] is False
        assert result["n_fail"] == 1
        assert result["files"][0]["status"] == "fail"

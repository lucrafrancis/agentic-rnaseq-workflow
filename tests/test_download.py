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
    """Write download metadata; runs get default sample labels unless they set their own."""
    for i, run in enumerate(runs):
        run.setdefault("sample_alias", f"GSM{i}")
        run.setdefault("sample_title", f"sample {i}")
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
            {"esearchresult": {"idlist": ["200123456"]}},
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
        mock_get.side_effect = lambda url: (
            "1. Some study\nSRA Run Selector: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP001\n"
            if "db=gds" in url else "^SAMPLE = GSM1\n!Sample_title = s1\n"
        )

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
        assert "curl" in script
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


# --- lookup failures vs. genuinely empty results ---


class TestLookupFailures:
    @patch("agents.download.tools._http_get_json")
    def test_lookup_failure_is_not_no_runs(self, mock_get_json, run_dir: Path):
        from agents.download.tools import LookupFailed
        mock_get_json.side_effect = LookupFailed("esearch: empty or non-JSON reply")
        result = resolve_accession("GSE110004")
        assert result["error"] == "lookup_failed"
        assert "retry" in result["message"]
        assert not SESSION.require_paths().download_metadata.exists()

    @patch("agents.download.tools._http_get_json")
    def test_esearch_error_payload_raises(self, mock_get_json):
        from agents.download.tools import LookupFailed
        mock_get_json.return_value = {"esearchresult": {"ERROR": "API rate limit exceeded"}}
        with pytest.raises(LookupFailed):
            _resolve_gse("GSE123")

    @patch("agents.download.tools._http_get_json")
    @patch("agents.download.tools._http_get")
    def test_empty_efetch_raises(self, mock_get, mock_get_json):
        from agents.download.tools import LookupFailed
        mock_get_json.return_value = {"esearchresult": {"idlist": ["1"]}}
        mock_get.return_value = ""
        with pytest.raises(LookupFailed):
            _resolve_gse("GSE123")

    @patch("agents.download.tools._http_get_json")
    @patch("agents.download.tools._http_get")
    def test_no_sra_link_is_no_runs(self, mock_get, mock_get_json, run_dir: Path):
        mock_get_json.return_value = {"esearchresult": {"idlist": ["1"]}}
        mock_get.side_effect = lambda url: (
            "1. A microarray study\nPlatform GPL570\n" if "db=gds" in url
            else "^SERIES = GSE123\n!Series_title = Arrays\n"
        )
        result = resolve_accession("GSE123")
        assert result["error"] == "no_runs_found"
        assert "reachable" in result["message"]

    def test_http_get_json_retries_empty_body(self, monkeypatch):
        import io
        from agents.download import tools
        bodies = [b"", b"", b'{"ok": 1}']
        monkeypatch.setattr(tools.urllib.request, "urlopen", lambda req, timeout: io.BytesIO(bodies.pop(0)))
        monkeypatch.setattr(tools.time, "sleep", lambda s: None)
        assert tools._http_get_json("https://example.org/x") == {"ok": 1}

    def test_http_get_json_gives_up(self, monkeypatch):
        import io
        from agents.download import tools
        monkeypatch.setattr(tools.urllib.request, "urlopen", lambda req, timeout: io.BytesIO(b""))
        monkeypatch.setattr(tools.time, "sleep", lambda s: None)
        with pytest.raises(tools.LookupFailed, match="empty or non-JSON"):
            tools._http_get_json("https://example.org/x")


# --- run subsetting ---


def _three_runs() -> None:
    _write_metadata([
        {
            "run_accession": f"SRR00{i}",
            "sample_alias": f"GSM{i}",
            "sample_title": title,
            "files": [{"filename": f"SRR00{i}.fastq.gz", "url": f"ftp://test/{i}", "md5": "x", "bytes": 1000}],
        }
        for i, title in enumerate(["Mock rep1", "Mock rep2", "CoV2 rep1"], start=1)
    ])


class TestRunSelection:
    def test_resolve_returns_sample_titles(self, run_dir: Path):
        with patch("agents.download.tools._resolve_study") as mock:
            mock.return_value = [{
                "run_accession": "SRR001", "library_layout": "PAIRED", "sample_alias": "GSM1",
                "sample_title": "Mock rep1", "files": [], "total_bytes": 0,
            }]
            result = resolve_accession("SRP001")
        assert result["runs"][0]["sample"] == "GSM1"
        assert result["runs"][0]["title"] == "Mock rep1"

    def test_parse_ena_run_keeps_sample_fields(self):
        run = _parse_ena_run({"run_accession": "SRR1", "sample_alias": "GSM1", "sample_title": "Mock rep1"})
        assert run["sample_alias"] == "GSM1"
        assert run["sample_title"] == "Mock rep1"

    def test_script_only_includes_selected_runs(self, run_dir: Path, tmp_path: Path):
        _three_runs()
        result = generate_download_script(str(tmp_path / "out"), runs=["SRR001", "SRR003"])
        assert result["n_files"] == 2
        assert result["n_runs_selected"] == 2
        assert result["n_runs_total"] == 3
        script = Path(result["script_path"]).read_text()
        assert "SRR002" not in script
        metadata = json.loads(SESSION.require_paths().download_metadata.read_text())
        assert metadata["selected_runs"] == ["SRR001", "SRR003"]

    def test_all_runs_when_omitted(self, run_dir: Path, tmp_path: Path):
        from agents.download.tools import selected_runs
        _three_runs()
        assert generate_download_script(str(tmp_path / "out"))["n_files"] == 3
        metadata = json.loads(SESSION.require_paths().download_metadata.read_text())
        assert metadata["selected_runs"] is None
        assert len(selected_runs(metadata)) == 3

    def test_selection_recorded_when_nothing_to_download(self, run_dir: Path, tmp_path: Path):
        out = tmp_path / "out"
        out.mkdir()
        (out / "SRR001.fastq.gz").write_bytes(b"data")
        md5 = _compute_md5(out / "SRR001.fastq.gz")
        _write_metadata([{"run_accession": "SRR001", "files": [
            {"filename": "SRR001.fastq.gz", "url": "ftp://t", "md5": md5, "bytes": 4}]}])
        result = generate_download_script(str(out), runs=["SRR001"])
        assert result["script_path"] is None
        metadata = json.loads(SESSION.require_paths().download_metadata.read_text())
        assert metadata["selected_runs"] == ["SRR001"]
        assert metadata["output_dir"] == str(out.resolve())

    @pytest.mark.parametrize("runs,fragment", [(["SRR999"], "SRR999"), ([], "empty")])
    def test_bad_runs(self, run_dir: Path, tmp_path: Path, runs, fragment):
        _three_runs()
        result = generate_download_script(str(tmp_path / "out"), runs=runs)
        assert result["error"] == "bad_runs"
        assert fragment in result["message"]
        assert check_existing_files(str(tmp_path), runs=runs)["error"] == "bad_runs"

    def test_check_existing_respects_selection(self, run_dir: Path, tmp_path: Path):
        _three_runs()
        result = check_existing_files(str(tmp_path), runs=["SRR002"])
        assert result["n_files"] == 1
        assert result["files"][0]["run"] == "SRR002"


class TestSuperSeries:
    SUPER_SOFT = (
        "^SERIES = GSE110004\n!Series_title = Umbrella study\n"
        "!Series_relation = SuperSeries of: GSE110000\n"
        "!Series_relation = SuperSeries of: GSE110003\n"
        "!Series_relation = BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA432544\n"
    )

    def _get(self, url):
        if "db=gds" in url:
            return "1. Umbrella study\n(Submitter supplied) This SuperSeries is composed of the SubSeries listed below.\n"
        if "acc=GSE110004" in url:
            return self.SUPER_SOFT
        if "acc=GSE110000" in url:
            return "^SERIES = GSE110000\n!Series_title = ChIP-seq of Rap1\n"
        if "acc=GSE110003" in url:
            return "^SERIES = GSE110003\n!Series_title = RNA-seq of Rap1 depletion\n"
        raise AssertionError(url)

    @patch("agents.download.tools._http_get_json")
    def test_superseries_lists_subseries(self, mock_get_json, run_dir: Path):
        mock_get_json.return_value = {"esearchresult": {"idlist": ["200110004"]}}
        with patch("agents.download.tools._http_get", side_effect=self._get):
            result = resolve_accession("GSE110004")
        assert result["error"] == "superseries"
        assert result["subseries"] == [
            {"accession": "GSE110000", "title": "ChIP-seq of Rap1"},
            {"accession": "GSE110003", "title": "RNA-seq of Rap1 depletion"},
        ]

    @patch("agents.download.tools._http_get_json")
    def test_superseries_soft_blocked_is_lookup_failure(self, mock_get_json, run_dir: Path):
        mock_get_json.return_value = {"esearchresult": {"idlist": ["200110004"]}}

        def get(url):
            return self._get(url) if "db=gds" in url else "<!doctype html><html>captcha"

        with patch("agents.download.tools._http_get", side_effect=get):
            result = resolve_accession("GSE110004")
        assert result["error"] == "lookup_failed"

    @patch("agents.download.tools._http_get_json")
    def test_subseries_titles_are_best_effort(self, mock_get_json, run_dir: Path):
        mock_get_json.return_value = {"esearchresult": {"idlist": ["200110004"]}}

        def get(url):
            return "<html>blocked" if "acc=GSE110000" in url else self._get(url)

        with patch("agents.download.tools._http_get", side_effect=get):
            result = resolve_accession("GSE110004")
        assert result["subseries"][0] == {"accession": "GSE110000"}
        assert result["subseries"][1]["title"] == "RNA-seq of Rap1 depletion"


# --- BioProject fallback and GEO labelling ---


class TestGeoFallbackAndLabelling:
    GDS_NO_SRA = "1. Recent study\nOrganism:\tHomo sapiens\nSeries\t\tAccession: GSE246386\n"
    SERIES_SOFT = (
        "^SERIES = GSE246386\n!Series_title = GFI1B\n"
        "!Series_relation = BioProject: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1032643\n"
    )
    SAMPLES_SOFT = (
        "^SAMPLE = GSM7868165\n!Sample_title = EV#1\n"
        "!Sample_relation = BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN1\n"
        "!Sample_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRX22242368\n"
        "^SAMPLE = GSM7868167\n!Sample_title = EV#3\n"
        "!Sample_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term=SRX22242370\n"
    )
    ENA = [
        {"run_accession": "SRR26539599", "experiment_accession": "SRX22242368",
         "sample_alias": "GSM7868165", "sample_title": "EV#1", "library_layout": "PAIRED"},
        {"run_accession": "SRR26539597", "experiment_accession": "SRX22242370",
         "sample_alias": "", "sample_title": "", "library_layout": "PAIRED"},
        {"run_accession": "SRR0", "experiment_accession": "SRX0",
         "sample_alias": "", "sample_title": "", "library_layout": "PAIRED"},
    ]

    def _get(self, url):
        if "db=gds" in url:
            return self.GDS_NO_SRA
        if "targ=self" in url:
            return self.SERIES_SOFT
        if "targ=gsm" in url:
            return self.SAMPLES_SOFT
        raise AssertionError(url)

    def _json(self, url):
        if "esearch" in url:
            return {"esearchresult": {"idlist": ["200246386"]}}
        assert "PRJNA1032643" in url
        return self.ENA

    def test_bioproject_fallback_and_labels(self):
        with patch("agents.download.tools._http_get", side_effect=self._get), \
             patch("agents.download.tools._http_get_json", side_effect=self._json):
            runs = {r["run_accession"]: r for r in _resolve_gse("GSE246386")}
        assert len(runs) == 3
        assert (runs["SRR26539597"]["sample_alias"], runs["SRR26539597"]["sample_title"]) == ("GSM7868167", "EV#3")
        assert runs["SRR0"]["sample_alias"] == ""  # no GEO match: left as ENA had it

    def test_labelling_blocked_is_lookup_failure(self, run_dir: Path):
        def get(url):
            return "<html>captcha" if "targ=gsm" in url else self._get(url)

        with patch("agents.download.tools._http_get", side_effect=get), \
             patch("agents.download.tools._http_get_json", side_effect=self._json):
            assert resolve_accession("GSE246386")["error"] == "lookup_failed"


class TestUnlabelledRuns:
    def test_script_refused_for_unlabelled_selection(self, run_dir: Path, tmp_path: Path):
        _write_metadata([
            {"run_accession": "SRR1", "files": [{"filename": "a.fq.gz", "url": "ftp://t", "md5": "x", "bytes": 1}]},
            {"run_accession": "SRR2", "sample_alias": "", "sample_title": "",
             "files": [{"filename": "b.fq.gz", "url": "ftp://t", "md5": "y", "bytes": 1}]},
        ])
        result = generate_download_script(str(tmp_path / "out"))
        assert result["error"] == "unlabelled_runs"
        assert "SRR2" in result["message"]
        assert not SESSION.require_paths().download_script.exists()
        # Selecting only labelled runs is fine
        assert generate_download_script(str(tmp_path / "out"), runs=["SRR1"])["n_files"] == 1

    def test_title_alone_is_not_enough(self):
        from agents.download.tools import unlabelled_runs
        runs = [{"run_accession": "SRR1", "sample_alias": "", "sample_title": "Mock rep1"},
                {"run_accession": "SRR2", "sample_alias": "GSM2", "sample_title": "CoV2 rep1"}]
        assert unlabelled_runs(runs) == ["SRR1"]

    def test_resolve_reports_unlabelled(self, run_dir: Path):
        with patch("agents.download.tools._resolve_study") as mock:
            mock.return_value = [{"run_accession": "SRR1", "library_layout": "PAIRED", "sample_alias": "",
                                  "sample_title": "", "files": [], "total_bytes": 0}]
            result = resolve_accession("SRP001")
        assert result["unlabelled_runs"] == ["SRR1"]

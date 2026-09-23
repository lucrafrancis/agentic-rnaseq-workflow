"""Tests for the sample sheet tools — offline, no API calls."""

from __future__ import annotations

from pathlib import Path

from agents.samplesheet.tools import (
    _check_samplesheet,
    _match_pairs,
    draft_samplesheet,
    match_pairs,
    read_metadata,
    sample_label_table,
    save_samplesheet,
    scan_fastqs,
    stage_fastqs,
    validate_samplesheet,
    write_report,
)
from core.session import SESSION


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
        result = _match_pairs(files)
        assert result["n_matched_pairs"] == 3
        assert result["n_unpaired"] == 0
        assert result["n_incomplete"] == 0

    def test_detects_unpaired(self):
        files = ["sampleA_R1.fastq.gz", "sampleB_R1.fastq.gz", "sampleB_R2.fastq.gz"]
        result = _match_pairs(files)
        assert result["n_matched_pairs"] == 1
        assert result["n_incomplete"] == 1

    def test_handles_underscore_numbering(self):
        files = ["sample_1.fastq.gz", "sample_2.fastq.gz"]
        result = _match_pairs(files)
        assert result["n_matched_pairs"] == 1

    def test_lanes_get_separate_pair_ids(self):
        files = ["A_S1_L001_R1.fq.gz", "A_S1_L001_R2.fq.gz", "A_S1_L002_R1.fq.gz", "A_S1_L002_R2.fq.gz"]
        matched = _match_pairs(files)["matched"]
        assert [m["pair_id"] for m in matched] == ["A_S1_L001", "A_S1_L002"]
        assert {m["sample"] for m in matched} == {"A"}

    def test_uses_scanned_files_and_stores_pairs(self, fastq_dir: Path):
        assert match_pairs()["error"] == "no_fastqs"
        scan_fastqs(str(fastq_dir))
        result = match_pairs()
        assert result["n_matched_pairs"] == 3
        assert SESSION.fastq_pairs["sampleA"]["fastq_2"] == str((fastq_dir / "sampleA_R2.fastq.gz").resolve())


def _scan_and_match(d: Path) -> None:
    scan_fastqs(str(d))
    match_pairs()


class TestDraftSamplesheet:
    def test_paths_come_from_pairs(self, fastq_dir: Path):
        _scan_and_match(fastq_dir)
        result = draft_samplesheet([{"pair_id": "sampleB", "sample": "ctrl_1", "strandedness": "reverse"}])
        assert result["rows"] == [{"sample": "ctrl_1", "pair_id": "sampleB", "strandedness": "reverse"}]
        assert result["unused_pair_ids"] == ["sampleA", "sampleC"]
        r1 = (fastq_dir / "sampleB_R1.fastq.gz").resolve()
        assert SESSION.samplesheet_draft.splitlines()[1] == f"ctrl_1,{r1},{str(r1).replace('R1', 'R2')},reverse"

    def test_strandedness_defaults_to_auto(self, fastq_dir: Path):
        _scan_and_match(fastq_dir)
        draft_samplesheet([{"pair_id": "sampleA", "sample": "a"}])
        assert SESSION.samplesheet_draft.splitlines()[1].endswith(",auto")

    def test_unknown_or_repeated_pair_id_rejects_draft(self, fastq_dir: Path):
        _scan_and_match(fastq_dir)
        result = draft_samplesheet([
            {"pair_id": "sampleA", "sample": "a"},
            {"pair_id": "sampleA", "sample": "b"},
            {"pair_id": "/data/sampleZ_R1.fastq.gz", "sample": "z"},
        ])
        assert result["error"] == "invalid_samples"
        assert len(result["errors"]) == 2
        assert SESSION.samplesheet_draft is None

    def test_requires_match_pairs(self):
        assert draft_samplesheet([{"pair_id": "x", "sample": "x"}])["error"] == "no_pairs"

    def test_single_end_unpaired_file(self, tmp_path: Path):
        d = tmp_path / "se"
        d.mkdir()
        (d / "ctrl.fastq.gz").touch()
        _scan_and_match(d)
        draft_samplesheet([{"pair_id": "ctrl.fastq.gz", "sample": "ctrl"}])
        assert SESSION.samplesheet_draft.splitlines()[1].endswith("ctrl.fastq.gz,,auto")


class TestValidateSamplesheet:
    def test_valid_sheet(self):
        sheet = (
            "sample,fastq_1,fastq_2,strandedness\n"
            "s1,/a_R1.fq.gz,/a_R2.fq.gz,auto\n"
            "s2,/b_R1.fq.gz,/b_R2.fq.gz,reverse\n"
        )
        result = _check_samplesheet(sheet)
        assert result["valid"] is True
        assert result["n_samples"] == 2
        assert result["errors"] == []

    def test_missing_columns(self):
        sheet = "sample,fastq_1\ns1,/a_R1.fq.gz\n"
        result = _check_samplesheet(sheet)
        assert result["valid"] is False
        assert any("Missing columns" in e for e in result["errors"])

    def test_duplicate_samples_warned(self):
        sheet = (
            "sample,fastq_1,fastq_2,strandedness\n"
            "s1,/a_R1.fq.gz,/a_R2.fq.gz,auto\n"
            "s1,/b_R1.fq.gz,/b_R2.fq.gz,auto\n"
        )
        result = _check_samplesheet(sheet)
        assert len(result["warnings"]) > 0

    def test_invalid_strandedness(self):
        sheet = (
            "sample,fastq_1,fastq_2,strandedness\n"
            "s1,/a_R1.fq.gz,/a_R2.fq.gz,banana\n"
        )
        result = _check_samplesheet(sheet)
        assert result["valid"] is False

    def test_empty_sheet(self):
        result = _check_samplesheet("")
        assert result["valid"] is False




class TestSaveSamplesheet:
    def test_saves_validated_draft(self, fastq_dir: Path):
        SESSION.begin_run("t")
        _scan_and_match(fastq_dir)
        draft_samplesheet([{"pair_id": "sampleA", "sample": "a"}])
        assert validate_samplesheet()["valid"] is True
        result = save_samplesheet()
        assert result["n_samples"] == 1
        assert SESSION.paths.samplesheet.read_text() == SESSION.samplesheet_draft

    def test_refuses_without_draft(self):
        SESSION.begin_run("t")
        assert save_samplesheet()["error"] == "no_draft"
        assert validate_samplesheet()["error"] == "no_draft"

    def test_refuses_unvalidated_draft(self, fastq_dir: Path):
        SESSION.begin_run("t")
        _scan_and_match(fastq_dir)
        draft_samplesheet([{"pair_id": "sampleA", "sample": "a"}])
        assert save_samplesheet()["error"] == "not_validated"
        assert not SESSION.paths.samplesheet.exists()

    def test_redraft_resets_validation(self, fastq_dir: Path):
        SESSION.begin_run("t")
        _scan_and_match(fastq_dir)
        draft_samplesheet([{"pair_id": "sampleA", "sample": "a"}])
        validate_samplesheet()
        draft_samplesheet([{"pair_id": "sampleA", "sample": "a"}, {"pair_id": "sampleB", "sample": "b"}])
        assert save_samplesheet()["error"] == "not_validated"

    def test_missing_file_fails_validation(self, fastq_dir: Path):
        SESSION.begin_run("t")
        _scan_and_match(fastq_dir)
        draft_samplesheet([{"pair_id": "sampleA", "sample": "a"}])
        (fastq_dir / "sampleA_R2.fastq.gz").unlink()
        result = validate_samplesheet()
        assert result["valid"] is False and "does not exist" in result["errors"][0]
        assert save_samplesheet()["error"] == "not_validated"

    def test_whitespace_in_sample_name_fails_validation(self, fastq_dir: Path):
        SESSION.begin_run("t")
        _scan_and_match(fastq_dir)
        draft_samplesheet([{"pair_id": "sampleA", "sample": "GFI1B rep1"}])
        assert validate_samplesheet()["valid"] is False


class TestStageFastqs:
    def test_links_and_updates_pairs(self, fastq_dir: Path):
        SESSION.begin_run("t")
        _scan_and_match(fastq_dir)
        result = stage_fastqs([{"pair_id": "sampleA", "sample": "ctrl_1"}])
        assert result["n_links"] == 2
        link = SESSION.paths.dir / "staged_fastqs" / "ctrl_1_R1.fastq.gz"
        assert link.is_symlink() and SESSION.fastq_pairs["sampleA"]["fastq_1"] == str(link)
        draft_samplesheet([{"pair_id": "sampleA", "sample": "ctrl_1"}])
        assert str(link) in SESSION.samplesheet_draft

    def test_rejects_colliding_names(self, fastq_dir: Path):
        SESSION.begin_run("t")
        _scan_and_match(fastq_dir)
        result = stage_fastqs([{"pair_id": "sampleA", "sample": "x"}, {"pair_id": "sampleB", "sample": "x"}])
        assert result["error"] == "invalid_samples"
        assert not (SESSION.paths.dir / "staged_fastqs").exists()


class TestWriteReport:
    def test_banner_when_no_samplesheet(self):
        SESSION.begin_run("t")
        result = write_report("# Report\nAll good.")
        assert "warning" in result
        assert (SESSION.paths.dir / "report.md").read_text().startswith("> ⚠ No sample sheet")

    def test_no_banner_when_saved(self, fastq_dir: Path):
        SESSION.begin_run("t")
        _scan_and_match(fastq_dir)
        draft_samplesheet([{"pair_id": "sampleA", "sample": "a"}])
        validate_samplesheet()
        save_samplesheet()
        result = write_report("# Report")
        assert "warning" not in result
        assert (SESSION.paths.dir / "report.md").read_text() == "# Report\n"


class TestSampleLabelTable:
    def test_table_from_files_metadata_and_design(self, tmp_path: Path):
        SESSION.begin_run("t")
        paths = SESSION.paths
        d = tmp_path / "fq"
        d.mkdir()
        for run in ("SRR1", "SRR2", "SRR3"):
            (d / f"{run}_1.fastq.gz").touch()
            (d / f"{run}_2.fastq.gz").touch()
        paths.sample_metadata.write_text(
            "run_accession,sample,title,library_layout\n"
            "SRR1,GSM1,EV#1,PAIRED\nSRR2,GSM2,GFI1B#1,PAIRED\nSRR3,GSM3,GFI1B#2,PAIRED\n"
        )
        paths.design.write_text("sample,condition\nEV_rep1,EV\nGFI1B_rep1,GFI1B\n")
        # EV_rep1 points at a staged symlink: the run must come from the resolved target
        link = tmp_path / "EV_rep1_R1.fastq.gz"
        link.symlink_to(d / "SRR1_1.fastq.gz")
        paths.samplesheet.write_text(
            "sample,fastq_1,fastq_2,strandedness\n"
            f"EV_rep1,{link},{d}/SRR1_2.fastq.gz,auto\n"
            f"GFI1B_rep1,{d}/SRR2_1.fastq.gz,{d}/SRR2_2.fastq.gz,auto\n"
        )
        table = sample_label_table(paths).splitlines()
        assert table[0].split() == ["sample", "condition", "run", "GEO", "sample", "GEO", "title"]
        assert table[2].split() == ["EV_rep1", "EV", "SRR1", "GSM1", "EV#1"]
        assert table[3].split() == ["GFI1B_rep1", "GFI1B", "SRR2", "GSM2", "GFI1B#1"]
        assert table[-1] == "⚠ Downloaded runs not in the sample sheet: SRR3"

    def test_none_without_metadata(self):
        SESSION.begin_run("t")
        SESSION.paths.samplesheet.write_text("sample,fastq_1,fastq_2,strandedness\n")
        assert sample_label_table(SESSION.paths) is None

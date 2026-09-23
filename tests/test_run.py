"""Tests for run.py stage orchestration: resume, skip-download, analysis-only.

Agents and approvals are monkeypatched — no API calls, no interactive input.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

import run
from core.approval import ApprovalResult
from core.session import SESSION


@pytest.fixture
def stubbed(monkeypatch: pytest.MonkeyPatch) -> dict[str, list]:
    """Stub everything past stage 1 and record which agents/approvals were called."""
    calls: dict[str, list] = {"approvals": [], "samplesheet_agent": [], "download_agent": [], "analysis": []}

    def approve(title, **kwargs):
        calls["approvals"].append(title)
        return ApprovalResult(approved=True)

    def samplesheet_agent(prompt):
        calls["samplesheet_agent"].append(prompt)
        SESSION.require_paths().samplesheet.write_text("sample,fastq_1\nA,a.fq.gz\n")

    monkeypatch.setattr(run, "present_for_approval", approve)
    monkeypatch.setattr("agents.samplesheet.loop.run_samplesheet_agent", samplesheet_agent)
    monkeypatch.setattr("agents.download.loop.run_download_agent", lambda p: calls["download_agent"].append(p))
    monkeypatch.setattr(run, "_submit_with_troubleshooting", lambda s: {"outdir": "results"})
    monkeypatch.setattr(run, "_run_analysis", lambda d, **kw: calls["analysis"].append(d))
    return calls


def _start(stages: list[str] | None = None) -> None:
    SESSION.begin_run("test_run")
    SESSION.require_paths().params_file.write_text("{}")  # skip submission agent
    for stage in stages or []:
        SESSION.mark_stage_complete(stage)


class TestRunState:
    def test_stages_persist_across_resume(self):
        SESSION.begin_run("state", source_prompt=Path("prompt.txt"))
        SESSION.mark_stage_complete("download")
        run_dir = SESSION.require_paths().dir

        SESSION.__init__()
        SESSION.resume_run(run_dir)
        assert SESSION.is_complete("download")
        assert not SESSION.is_complete("samplesheet")
        assert SESSION.source_prompt == str(Path("prompt.txt").resolve())

    def test_resume_without_state_file(self):
        SESSION.begin_run("legacy")
        run_dir = SESSION.require_paths().dir
        SESSION.require_paths().state_file.unlink()

        SESSION.__init__()
        SESSION.resume_run(run_dir)
        assert SESSION.stages_completed == []


class TestSamplesheetStage:
    def test_unapproved_samplesheet_is_reapproved_not_regenerated(self, stubbed):
        _start()
        SESSION.require_paths().samplesheet.write_text("sample,fastq_1\nX,x.fq.gz\n")

        run._run_pipeline("FASTQs in ./fastqs")

        assert stubbed["samplesheet_agent"] == []
        assert stubbed["approvals"] == ["Sample sheet"]
        assert SESSION.is_complete("samplesheet")

    def test_approved_samplesheet_is_skipped(self, stubbed):
        _start(["samplesheet"])
        SESSION.require_paths().samplesheet.write_text("sample,fastq_1\nX,x.fq.gz\n")

        run._run_pipeline("FASTQs in ./fastqs")

        assert stubbed["samplesheet_agent"] == []
        assert stubbed["approvals"] == []

    def test_rejected_samplesheet_is_not_marked(self, stubbed, monkeypatch):
        _start()
        monkeypatch.setattr(run, "present_for_approval", lambda title, **kw: ApprovalResult(approved=False))
        with pytest.raises(SystemExit):
            run._run_pipeline("FASTQs in ./fastqs")
        assert not SESSION.is_complete("samplesheet")


class TestDownloadStage:
    def test_skip_download_flag(self, stubbed):
        _start()
        run._run_pipeline("GSE123456, FASTQs in ./fastqs", skip_download=True)
        assert stubbed["download_agent"] == []
        assert SESSION.is_complete("download")

    def test_no_accession_no_download(self, stubbed):
        _start()
        run._run_pipeline("FASTQs in ./fastqs")
        assert stubbed["download_agent"] == []

    def test_completed_download_is_skipped(self, stubbed):
        _start(["download"])
        run._run_pipeline("GSE123456")
        assert stubbed["download_agent"] == []

    def test_metadata_without_completion_reruns_stage(self, stubbed):
        _start()
        SESSION.require_paths().download_metadata.write_text(json.dumps({"runs": []}))
        run._run_pipeline("GSE123456")  # metadata lacks output_dir -> agent must re-run
        assert stubbed["download_agent"] == ["GSE123456"]

    def test_interrupted_download_regenerates_without_agent(self, stubbed, monkeypatch, tmp_path):
        _start()
        paths = SESSION.require_paths()
        paths.download_script.write_text("# stale\n")
        paths.download_metadata.write_text(json.dumps({"runs": [], "output_dir": str(tmp_path)}))
        regenerated = []
        monkeypatch.setattr(
            "agents.download.tools.generate_download_script",
            lambda d: regenerated.append(d) or {"n_files": 0, "script_path": None},
        )

        run._run_pipeline("GSE123456")

        assert regenerated == [str(tmp_path)]
        assert stubbed["download_agent"] == []
        assert not paths.download_script.exists()  # stale script removed, never approved
        assert SESSION.is_complete("download")

    def test_rejected_download_exits(self, stubbed, monkeypatch, tmp_path):
        _start()
        paths = SESSION.require_paths()
        paths.download_metadata.write_text(json.dumps({"runs": [], "output_dir": str(tmp_path)}))

        def fake_generate(d):
            paths.download_script.write_text("curl ...\n")
            return {"n_files": 1}

        monkeypatch.setattr("agents.download.tools.generate_download_script", fake_generate)
        monkeypatch.setattr(run, "present_for_approval", lambda title, **kw: ApprovalResult(approved=False))
        with pytest.raises(SystemExit, match="--skip-download"):
            run._run_pipeline("GSE123456")
        assert not SESSION.is_complete("download")
        assert not SESSION.is_complete("samplesheet")


class TestAnalyzeMode:
    def test_analyze_prompt_creates_run_without_results(self, stubbed, monkeypatch, tmp_path):
        prompt = tmp_path / "myproject" / "prompt.txt"
        prompt.parent.mkdir()
        prompt.write_text("GSE123456. Count matrix at ./counts.csv\n")
        monkeypatch.setattr(sys, "argv", ["run.py", "--analyze", str(prompt)])

        run.main()

        assert stubbed["analysis"] == [None]
        assert stubbed["download_agent"] == []
        paths = SESSION.require_paths()
        assert paths.dir.name.endswith("_myproject")
        assert (paths.dir / "prompt.txt").read_text() == prompt.read_text()
        assert SESSION.source_prompt == str(prompt.resolve())

    def test_analyze_run_dir_without_results(self, stubbed, monkeypatch):
        SESSION.begin_run("counts_only")
        run_dir = SESSION.require_paths().dir
        (run_dir / "prompt.txt").write_text("Count matrix at ./counts.csv\n")
        SESSION.__init__()
        monkeypatch.setattr(sys, "argv", ["run.py", "--analyze", str(run_dir)])

        run.main()

        assert stubbed["analysis"] == [None]

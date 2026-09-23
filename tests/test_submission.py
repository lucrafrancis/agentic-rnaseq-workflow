"""Tests for the submission handler — offline, no nextflow needed."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from agents.submission.errors import Diagnosis, diagnose_log, format_diagnoses
from agents.submission.params import SubmissionParams, default_params
from agents.submission.tools import check_resources, configure_submission
from agents.submission.troubleshoot import _run_safe_command
from core.session import SESSION


class TestSubmissionParams:
    def test_default_params(self):
        params = default_params("/path/to/sheet.csv", "/path/to/outdir")
        assert params.input_samplesheet == "/path/to/sheet.csv"
        assert params.outdir == "/path/to/outdir"
        assert params.genome == "GRCh38"
        assert params.profile == "docker"

    def test_to_script(self):
        params = default_params("/sheet.csv", "/out")
        script = params.to_script("/params.yml")
        assert "nextflow run nf-core/rnaseq" in script
        assert "-params-file" in script

    def test_to_script_with_resume_and_config(self):
        params = default_params("/sheet.csv", "/out")
        script = params.to_script("/params.yml", config_file="/custom.config", resume=True)
        assert "-resume" in script
        assert "-c /custom.config" in script

    def test_to_dict(self):
        params = default_params("/sheet.csv", "/out")
        d = params.to_dict()
        assert d["input"] == "/sheet.csv"
        assert d["genome"] == "GRCh38"


class TestDiagnoseLog:
    def test_detects_genome_error(self, tmp_path: Path):
        log = tmp_path / "nextflow.log"
        log.write_text("ERROR ~ Unknown genome 'GRCh99' specified\n")
        diagnoses = diagnose_log(str(log))
        assert len(diagnoses) == 1
        assert diagnoses[0].category == "genome_not_found"

    def test_detects_missing_files(self, tmp_path: Path):
        log = tmp_path / "nextflow.log"
        log.write_text("ERROR: No such file /data/missing_R1.fastq.gz\n")
        diagnoses = diagnose_log(str(log))
        assert any(d.category == "missing_files" for d in diagnoses)

    def test_detects_oom(self, tmp_path: Path):
        log = tmp_path / "nextflow.log"
        log.write_text("Process `STAR_ALIGN` exceeded memory limit 32 GB\n")
        diagnoses = diagnose_log(str(log))
        assert any(d.category == "out_of_memory" for d in diagnoses)

    def test_clean_log_returns_empty(self, tmp_path: Path):
        log = tmp_path / "nextflow.log"
        log.write_text("Launching nf-core/rnaseq\nCompleted successfully.\n")
        diagnoses = diagnose_log(str(log))
        assert diagnoses == []

    def test_missing_log_returns_empty(self):
        diagnoses = diagnose_log("/nonexistent/nextflow.log")
        assert diagnoses == []

    def test_format_diagnoses_empty(self):
        assert "No known error" in format_diagnoses([])

    def test_format_diagnoses_with_items(self):
        d = Diagnosis(category="test", line="ERROR line", suggestion="Fix it")
        out = format_diagnoses([d])
        assert "1 issue" in out
        assert "Fix it" in out


class TestFromDict:
    def test_round_trip(self):
        params = SubmissionParams(
            input_samplesheet="/sheet.csv",
            outdir="/out",
            genome="GRCm39",
            profile="singularity",
            extra_args={"skip_alignment": True, "max_memory": "30.GB"},
        )
        rebuilt = SubmissionParams.from_dict(params.to_dict())
        assert rebuilt.genome == "GRCm39"
        assert rebuilt.profile == "singularity"
        assert rebuilt.extra_args["skip_alignment"] is True
        assert rebuilt.extra_args["max_memory"] == "30.GB"
        assert rebuilt.input_samplesheet == "/sheet.csv"

    def test_load_from_file(self, tmp_path: Path):
        data = {
            "pipeline": "nf-core/rnaseq",
            "revision": "3.26.0",
            "input": "/sheet.csv",
            "outdir": "/out",
            "genome": "GRCh38",
            "profile": "docker",
            "skip_alignment": True,
        }
        f = tmp_path / "params.json"
        f.write_text(json.dumps(data))
        params = SubmissionParams.load(f)
        assert params.extra_args["skip_alignment"] is True
        assert params.genome == "GRCh38"


class TestConfigureSubmission:
    @pytest.fixture(autouse=True)
    def machine(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """A 36 GB Mac whose Docker VM has 23.4 GB — the GSE246386 run's machine."""
        monkeypatch.setattr("agents.submission.tools._memory_gb", lambda: 36.0)
        monkeypatch.setattr("agents.submission.tools._docker_info",
                            lambda: {"available": True, "memory_gb": 23.4, "cpus": 12})

    def test_rejects_memory_above_docker(self):
        SESSION.begin_run("test_docker_mem")
        paths = SESSION.require_paths()
        result = configure_submission(genome="GRCh38", extra_params={"max_memory": "32.GB"})
        assert result["error"] == "resource_limit"
        assert "23.4 GB available to Docker" in result["message"]
        assert not paths.params_file.exists()

    def test_accepts_memory_below_docker(self):
        SESSION.begin_run("test_docker_mem_ok")
        result = configure_submission(genome="GRCh38", extra_params={"max_memory": "21.GB"})
        assert "error" not in result

    def test_docker_limit_only_applies_to_docker_profile(self):
        SESSION.begin_run("test_singularity_mem")
        result = configure_submission(genome="GRCh38", profile="singularity", extra_params={"max_memory": "32.GB"})
        assert "error" not in result

    def test_rejects_memory_above_host(self):
        SESSION.begin_run("test_host_mem")
        result = configure_submission(genome="GRCh38", profile="singularity", extra_params={"max_memory": "64 GB"})
        assert "36.0 GB of RAM" in result["message"]

    def test_rejects_unparseable_memory(self):
        SESSION.begin_run("test_bad_mem")
        result = configure_submission(genome="GRCh38", extra_params={"max_memory": "lots"})
        assert result["error"] == "resource_limit"

    def test_basic(self, tmp_path: Path):
        SESSION.begin_run("test_configure")
        paths = SESSION.require_paths()
        paths.samplesheet.write_text("sample,fastq_1,fastq_2,strandedness\n")

        result = configure_submission(genome="GRCh38")
        assert result["genome"] == "GRCh38"
        assert paths.params_file.is_file()

        saved = json.loads(paths.params_file.read_text())
        assert saved["genome"] == "GRCh38"

    def test_extra_params(self, tmp_path: Path):
        SESSION.begin_run("test_configure_extras")
        paths = SESSION.require_paths()
        paths.samplesheet.write_text("sample,fastq_1,fastq_2,strandedness\n")

        result = configure_submission(
            genome="GRCm39",
            profile="singularity",
            extra_params={"skip_alignment": True, "max_memory": "16.GB"},
        )
        assert result["genome"] == "GRCm39"
        assert result["profile"] == "singularity"
        assert result["extra_params"]["skip_alignment"] is True

        saved = json.loads(paths.params_file.read_text())
        assert saved["skip_alignment"] is True
        assert saved["max_memory"] == "16.GB"

    def test_skip_alignment_forces_pseudo_aligner(self, tmp_path: Path):
        SESSION.begin_run("test_configure_pseudo")
        SESSION.require_paths().samplesheet.write_text("sample,fastq_1,fastq_2,strandedness\n")

        configure_submission(
            genome="GRCh38",
            extra_params={"skip_alignment": True},
        )
        saved = json.loads(SESSION.require_paths().params_file.read_text())
        assert saved["skip_alignment"] is True
        assert saved["pseudo_aligner"] == "salmon"

    def test_skip_alignment_respects_explicit_pseudo_aligner(self, tmp_path: Path):
        SESSION.begin_run("test_configure_pseudo_explicit")
        SESSION.require_paths().samplesheet.write_text("sample,fastq_1,fastq_2,strandedness\n")

        configure_submission(
            genome="GRCh38",
            extra_params={"skip_alignment": True, "pseudo_aligner": "salmon"},
        )
        saved = json.loads(SESSION.require_paths().params_file.read_text())
        assert saved["pseudo_aligner"] == "salmon"

    def test_coerces_string_booleans(self, tmp_path: Path):
        SESSION.begin_run("test_configure_coerce")
        SESSION.require_paths().samplesheet.write_text("sample,fastq_1,fastq_2,strandedness\n")

        result = configure_submission(
            genome="GRCh38",
            extra_params={"skip_alignment": "true", "skip_trimming": "false"},
        )
        saved = json.loads(SESSION.require_paths().params_file.read_text())
        assert saved["skip_alignment"] is True
        assert saved["skip_trimming"] is False


class TestSchemaConsistency:
    """Check all agent schemas (samplesheet, download, submission)."""

    def _all_schemas(self):
        from agents.samplesheet.schemas import TOOL_FUNCTIONS as SS_FN, TOOL_SCHEMAS as SS_SC
        from agents.download.schemas import TOOL_FUNCTIONS as DL_FN, TOOL_SCHEMAS as DL_SC
        from agents.submission.schemas import TOOL_FUNCTIONS as SB_FN, TOOL_SCHEMAS as SB_SC
        return [(SS_SC, SS_FN), (DL_SC, DL_FN), (SB_SC, SB_FN)]

    def test_every_schema_has_a_callable_function(self):
        for schemas, functions in self._all_schemas():
            for schema in schemas:
                name = schema["name"]
                assert name in functions, f"{name} exposed to Claude but not dispatchable"
                assert callable(functions[name])

    def test_schemas_are_well_formed(self):
        for schemas, _ in self._all_schemas():
            for schema in schemas:
                assert {"name", "description", "input_schema"} <= schema.keys()
                assert schema["input_schema"]["type"] == "object"

    def test_no_duplicate_tool_names(self):
        for schemas, _ in self._all_schemas():
            names = [s["name"] for s in schemas]
            assert len(names) == len(set(names))


class TestLimitParsing:
    def test_units(self):
        from agents.submission.tools import _limit_gb
        assert _limit_gb("32.GB") == 32
        assert _limit_gb("32 GB") == 32
        assert _limit_gb("1.5.GB") == 1.5
        assert _limit_gb("512.MB") == 0.5
        assert _limit_gb("1.TB") == 1024
        assert _limit_gb("32") is None


class TestCheckResources:
    def test_returns_structured_info(self):
        result = check_resources()
        assert "cpus" in result
        assert "memory_gb" in result
        assert "disk" in result
        assert "docker" in result
        assert result["cpus"] is None or result["cpus"] > 0
        assert result["memory_gb"] is None or result["memory_gb"] > 0


class TestRunSafeCommand:
    def test_allowed_binary(self):
        result = _run_safe_command("uname -s")
        assert result["returncode"] == 0
        assert result["stdout"]

    def test_rejects_shell_metacharacters(self):
        result = _run_safe_command("uname -s; rm -rf /")
        assert result["error"] == "blocked"

    def test_rejects_pipe(self):
        result = _run_safe_command("cat /etc/passwd | grep root")
        assert result["error"] == "blocked"

    def test_rejects_subshell(self):
        result = _run_safe_command("$(whoami)")
        assert result["error"] == "blocked"

    def test_rejects_disallowed_binary(self):
        result = _run_safe_command("rm -rf /")
        assert result["error"] == "blocked"
        assert "not in the allowed" in result["message"]

    def test_rejects_path_traversal(self):
        result = _run_safe_command("/usr/bin/rm -rf /")
        assert result["error"] == "blocked"

    def test_empty_command(self):
        result = _run_safe_command("")
        assert result["error"] == "blocked"


class TestWriteNfConfig:
    def test_quotes_string_values(self, tmp_path: Path):
        params = SubmissionParams(
            extra_args={"max_memory": "30.GB", "max_cpus": 8, "max_time": "12.h"},
        )
        config_path = tmp_path / "custom.config"
        params.write_nf_config(config_path)
        content = config_path.read_text()
        assert "'30.GB'" in content
        assert "'12.h'" in content
        assert "cpus: 8" in content

    def test_integer_values_unquoted(self, tmp_path: Path):
        params = SubmissionParams(extra_args={"max_cpus": 4})
        config_path = tmp_path / "custom.config"
        params.write_nf_config(config_path)
        content = config_path.read_text()
        assert "cpus: 4" in content
        assert "'" not in content


class TestTroubleshootResourceCheck:
    def test_rejects_fix_above_docker_memory(self, monkeypatch: pytest.MonkeyPatch):
        from types import SimpleNamespace

        from agents.submission import troubleshoot
        from tests.test_loop import FakeClient, _response

        monkeypatch.setattr("agents.submission.tools._memory_gb", lambda: 36.0)
        monkeypatch.setattr("agents.submission.tools._docker_info", lambda: {"available": True, "memory_gb": 23.4})
        SESSION.begin_run("test_ts_mem")

        def fix(tool_id, memory):
            return SimpleNamespace(type="tool_use", id=tool_id, name="propose_fix", input={
                "root_cause": "OOM", "fix_description": "more memory", "category": "parameter_change",
                "parameter_changes": {"max_memory": memory},
            })

        client = FakeClient([_response([fix("t1", "32.GB")], "tool_use"),
                             _response([fix("t2", "21.GB")], "tool_use")])
        monkeypatch.setattr(troubleshoot.anthropic, "Anthropic", lambda: client)
        monkeypatch.setattr("builtins.input", lambda *a: "a")

        params = SubmissionParams(input_samplesheet="s.csv", outdir="out", genome="GRCh38")
        proposal = troubleshoot.diagnose_and_propose({"log_path": ""}, params, 1, [])

        assert proposal["parameter_changes"] == {"max_memory": "21.GB"}
        results = [c for m in client.calls[1]["messages"] if isinstance(m["content"], list)
                   for c in m["content"] if isinstance(c, dict) and c.get("tool_use_id") == "t1"]
        rejection = results[0]
        assert rejection["is_error"] and "available to Docker" in rejection["content"]



def _publish_reference(results: Path, gtfs=("genes.filtered.gtf",), index_files=("info.json", "versionInfo.json")) -> None:
    """What nf-core 3.26.0 publishes to results/genome/ with save_reference (Salmon-only)."""
    genome = results / "genome"
    (genome / "index" / "salmon").mkdir(parents=True)
    for g in gtfs:
        (genome / g).write_text("chr1\tgtf\n")
    (genome / "genome.transcripts.fa").write_text(">tx1\nACGT\n")
    (genome / "genome.fa.fai").write_text("chr1\t4\n")
    for f in index_files:
        (genome / "index" / "salmon" / f).write_text("{}")
    (results / "pipeline_info").mkdir()
    (results / "pipeline_info" / "nf_core_rnaseq_software_mqc_versions.yml").write_text(
        "SALMON_INDEX:\n  salmon: 1.10.3\nWorkflow:\n  nf-core/rnaseq: v3.26.0\n")


class TestReferenceCache:
    SALMON_ONLY = {"skip_alignment": True, "max_memory": "21.GB"}

    @pytest.fixture(autouse=True)
    def machine(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setattr("agents.submission.tools._memory_gb", lambda: 36.0)
        monkeypatch.setattr("agents.submission.tools._docker_info", lambda: {"available": True, "memory_gb": 23.4})

    def test_no_cache_builds_and_saves(self):
        SESSION.begin_run("ref_build")
        result = configure_submission(genome="GRCh38", extra_params=dict(self.SALMON_ONLY))
        assert result["extra_params"]["save_reference"] is True
        params = SubmissionParams.load(SESSION.paths.params_file)
        assert params.reference is None
        params.write_nf_params(SESSION.paths.nf_params)
        assert 'genome: "GRCh38"' in SESSION.paths.nf_params.read_text()

    def test_alignment_runs_untouched(self):
        SESSION.begin_run("ref_star")
        result = configure_submission(genome="GRCh38", extra_params={"max_memory": "21.GB"})
        assert "save_reference" not in result["extra_params"]

    def test_cache_then_reuse(self, tmp_path: Path):
        from agents.submission.reference import cache_reference, find_reference

        results = tmp_path / "results"
        _publish_reference(results)
        cached = cache_reference("GRCh38", "3.26.0", results, "run_1")
        assert cached["status"] == "cached"
        ref = find_reference("GRCh38", "3.26.0")
        info = json.loads((ref / "reference.json").read_text())
        assert info["salmon_version"] == "1.10.3" and info["built_by_run"] == "run_1"
        assert info["copied_from"]["gtf"] == "genes.filtered.gtf"
        assert not list(ref.parent.glob(".*tmp*"))
        assert cache_reference("GRCh38", "3.26.0", results, "run_2")["status"] == "already_cached"

        SESSION.begin_run("ref_reuse")
        result = configure_submission(genome="GRCh38", extra_params=dict(self.SALMON_ONLY))
        assert "save_reference" not in result["extra_params"]
        assert "reusing cached" in result["reference"]
        params = SubmissionParams.load(SESSION.paths.params_file)
        assert params.reference == str(ref)
        params.write_nf_params(SESSION.paths.nf_params)
        yml = SESSION.paths.nf_params.read_text()
        assert "genome:" not in yml
        assert f'salmon_index: "{ref / "salmon_index"}"' in yml
        assert f'gtf: "{ref / "genes.gtf"}"' in yml and f'transcript_fasta: "{ref / "transcripts.fa"}"' in yml

    def test_other_revision_not_reused(self, tmp_path: Path):
        from agents.submission.reference import cache_reference, find_reference

        _publish_reference(tmp_path / "results")
        cache_reference("GRCh38", "3.25.0", tmp_path / "results", "old")
        assert find_reference("GRCh38", "3.26.0") is None

    def test_user_reference_files_respected(self):
        SESSION.begin_run("ref_user")
        result = configure_submission(genome="GRCh38", extra_params={**self.SALMON_ONLY, "salmon_index": "/my/index"})
        assert "save_reference" not in result["extra_params"]

    def test_refuses_incomplete_index(self, tmp_path: Path):
        from agents.submission.reference import cache_reference, reference_dir

        _publish_reference(tmp_path / "results", index_files=("info.json",))
        result = cache_reference("GRCh38", "3.26.0", tmp_path / "results", "r")
        assert result["status"] == "not_cached" and "versionInfo.json" in result["problems"][0]
        assert not reference_dir("GRCh38", "3.26.0").exists()

    def test_refuses_ambiguous_gtf(self, tmp_path: Path):
        from agents.submission.reference import cache_reference

        _publish_reference(tmp_path / "results", gtfs=("a.gtf", "b.gtf"))
        assert cache_reference("GRCh38", "3.26.0", tmp_path / "results", "r")["status"] == "not_cached"

    def test_run_caches_after_success(self, tmp_path: Path):
        import run
        from agents.submission.reference import find_reference

        SESSION.begin_run("ref_run")
        results = tmp_path / "results"
        _publish_reference(results)
        params = SubmissionParams(input_samplesheet="s.csv", outdir=str(results), genome="GRCh38",
                                  extra_args={"skip_alignment": True, "pseudo_aligner": "salmon", "save_reference": True})
        run._cache_reference(params)
        assert find_reference("GRCh38", params.revision) is not None

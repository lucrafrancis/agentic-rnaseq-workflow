"""Tests for the submission handler — offline, no nextflow needed."""

from __future__ import annotations

import json
from pathlib import Path

from agents.submission.errors import Diagnosis, diagnose_log, format_diagnoses
from agents.submission.params import SubmissionParams, default_params
from agents.submission.tools import configure_submission
from core.session import SESSION


class TestSubmissionParams:
    def test_default_params(self):
        params = default_params("/path/to/sheet.csv", "/path/to/outdir")
        assert params.input_samplesheet == "/path/to/sheet.csv"
        assert params.outdir == "/path/to/outdir"
        assert params.genome == "GRCh38"
        assert params.profile == "docker"

    def test_to_command_string(self):
        params = default_params("/sheet.csv", "/out")
        cmd = params.to_command_string()
        assert "nextflow run nf-core/rnaseq" in cmd
        assert "--input /sheet.csv" in cmd
        assert "--outdir /out" in cmd

    def test_extra_args(self):
        params = SubmissionParams(
            input_samplesheet="/sheet.csv",
            outdir="/out",
            extra_args={"aligner": "star_salmon", "skip_qc": True},
        )
        args = params.to_nextflow_args()
        assert "--aligner" in args
        assert "star_salmon" in args
        assert "--skip_qc" in args

    def test_bool_false_excluded(self):
        params = SubmissionParams(
            input_samplesheet="/sheet.csv",
            outdir="/out",
            extra_args={"skip_qc": False},
        )
        args = params.to_nextflow_args()
        assert "--skip_qc" not in args

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
    def test_basic(self, tmp_path: Path):
        SESSION.begin_run("test_configure")
        paths = SESSION.require_paths()
        paths.samplesheet.write_text("sample,fastq_1,fastq_2,strandedness\n")

        result = configure_submission(genome="GRCh38")
        assert "command" in result
        assert "GRCh38" in result["command"]
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

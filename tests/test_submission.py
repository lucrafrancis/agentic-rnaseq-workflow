"""Tests for the submission handler — offline, no nextflow needed."""

from __future__ import annotations

from pathlib import Path

from agents.submission.errors import Diagnosis, diagnose_log, format_diagnoses
from agents.submission.params import SubmissionParams, default_params
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


class TestSchemaConsistency:
    def test_every_schema_has_a_callable_function(self):
        from agents.samplesheet.schemas import TOOL_FUNCTIONS, TOOL_SCHEMAS

        for schema in TOOL_SCHEMAS:
            name = schema["name"]
            assert name in TOOL_FUNCTIONS, f"{name} exposed to Claude but not dispatchable"
            assert callable(TOOL_FUNCTIONS[name])

    def test_schemas_are_well_formed(self):
        from agents.samplesheet.schemas import TOOL_SCHEMAS

        for schema in TOOL_SCHEMAS:
            assert {"name", "description", "input_schema"} <= schema.keys()
            assert schema["input_schema"]["type"] == "object"

    def test_no_duplicate_tool_names(self):
        from agents.samplesheet.schemas import TOOL_SCHEMAS

        names = [s["name"] for s in TOOL_SCHEMAS]
        assert len(names) == len(set(names))

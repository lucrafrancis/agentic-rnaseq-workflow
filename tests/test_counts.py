"""Tests for the GEO count-matrix tools — offline, _fetch is mocked.

Fixtures mirror real GEO layouts: R write.table output (GSE157852), comma-separated
with a BOM (GSE164073), featureCounts with comment + annotation columns, and
NCBI-generated counts keyed by Entrez ID with GSM columns.
"""

from __future__ import annotations

import gzip
import json

import pandas as pd
import pytest

from agents.download import counts
from core import config
from core.session import SESSION

SOFT = """\
^SAMPLE = GSM1
!Sample_title = Mock rep1
!Sample_organism_ch1 = Homo sapiens
!Sample_characteristics_ch1 = tissue: cornea
!Sample_characteristics_ch1 = infection: mock
!Sample_library_strategy = RNA-Seq
^SAMPLE = GSM2
!Sample_title = Mock rep2
!Sample_organism_ch1 = Homo sapiens
!Sample_characteristics_ch1 = tissue: limbus
!Sample_characteristics_ch1 = infection: mock
!Sample_library_strategy = RNA-Seq
^SAMPLE = GSM3
!Sample_title = CoV2 rep1
!Sample_organism_ch1 = Homo sapiens
!Sample_characteristics_ch1 = tissue: cornea
!Sample_characteristics_ch1 = infection: SARS-CoV-2
!Sample_library_strategy = RNA-Seq
^SAMPLE = GSM4
!Sample_title = CoV2 rep2
!Sample_organism_ch1 = Homo sapiens
!Sample_characteristics_ch1 = tissue: limbus
!Sample_characteristics_ch1 = infection: SARS-CoV-2
!Sample_library_strategy = RNA-Seq
"""

FTP_INDEX = """<pre>Name                             Last modified      Size  <hr><a href="/geo/series/GSE1nnn/GSE1234/">Parent Directory</a>   -
<a href="GSE1234_counts.txt.gz">GSE1234_counts.txt.gz</a>                2020-09-11 14:58  416K
<a href="GSE1234_bom.csv.gz">GSE1234_bom.csv.gz</a>                2020-09-11 14:58  416K
<a href="GSE1234_fc.txt">GSE1234_fc.txt</a>                2020-09-11 14:58  416K
<a href="GSE1234_FPKM.xlsx">GSE1234_FPKM.xlsx</a>                2020-09-11 14:58  1.0M
<a href="GSE1234_RAW.tar">GSE1234_RAW.tar</a>                2020-09-11 14:58  1.0M
<a href="filelist.txt">filelist.txt</a>                     2026-02-03 00:42  557
<hr></pre>"""

FILELIST = (
    "#Archive/File\tName\tTime\tSize\tType\n"
    "Archive\tGSE1234_RAW.tar\t02/03/2026 00:42:34\t1013760\tTAR\n"
    "File\tGSM1_counts.txt.gz\t10/20/2023 08:02:03\t166519\tTXT\n"
)

DOWNLOAD_PAGE = (
    '<a href="/geo/download/?type=rnaseq_counts&amp;acc=GSE1234&amp;format=file&amp;'
    'file=GSE1234_raw_counts_GRCh38.p13_NCBI.tsv.gz">raw</a>'
    '<a href="/geo/download/?type=rnaseq_counts&amp;acc=GSE1234&amp;format=file&amp;'
    'file=GSE1234_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz">tpm</a>'
)

R_STYLE = "S_mock_1 S_mock_2 S_cov_1 S_cov_2\nGAPDH 100 120 90 95\nACTB 50 55 60 58\nGAPDH 1 1 1 1\n"
BOM_CSV = "﻿Gene,MW1_mock,MW2_mock,MW3_cov,MW4_cov\nA1BG,91,131,86,77\nA1BG-AS1,292,284,271,232\n"
FEATURECOUNTS = (
    "# Program:featureCounts v2.0.1\n"
    "Geneid\tChr\tStart\tEnd\tStrand\tLength\tMock rep1\tMock rep2\tCoV2 rep1\tCoV2 rep2\n"
    "ENSG00000111640.15\tchr12\t1\t10\t+\t1000\t10\t20\t30\t40\n"
    "ENSG00000075624.17\tchr7\t1\t10\t-\t2000\t5\t6\t7\t8\n"
)
NCBI = "GeneID\tGSM1\tGSM2\tGSM3\tGSM4\n2597\t100\t120\t90\t95\n60\t50\t55\t60\t58\n"
GENE_INFO = (
    "#tax_id\tGeneID\tSymbol\tLocusTag\tSynonyms\tdbXrefs\n"
    "9606\t2597\tGAPDH\t-\t-\tHGNC:4141|Ensembl:ENSG00000111640\n"
    "9606\t60\tACTB\t-\t-\tHGNC:132|Ensembl:ENSG00000075624\n"
)


@pytest.fixture
def geo(monkeypatch: pytest.MonkeyPatch, tmp_path):
    """Mock NCBI responses keyed by URL substring; return the list of fetched URLs."""
    monkeypatch.setattr(config, "REFERENCE_DIR", tmp_path / "ref")
    responses = {
        "targ=gsm": SOFT.encode(),
        "suppl/filelist.txt": FILELIST.encode(),
        "suppl/GSE1234_counts.txt.gz": gzip.compress(R_STYLE.encode()),
        "suppl/GSE1234_bom.csv.gz": gzip.compress(BOM_CSV.encode("utf-8")),
        "suppl/GSE1234_fc.txt": FEATURECOUNTS.encode(),
        "suppl/": FTP_INDEX.encode(),
        "geo/download/?acc=": DOWNLOAD_PAGE.encode(),
        "file=GSE1234_raw_counts": gzip.compress(NCBI.encode()),
        "Homo_sapiens.gene_info.gz": gzip.compress(GENE_INFO.encode()),
    }
    fetched: list[str] = []

    def fake_fetch(url, timeout=60):
        fetched.append(url)
        for key, body in responses.items():
            if key in url:
                return body
        raise AssertionError(f"unexpected URL {url}")

    monkeypatch.setattr(counts, "_fetch", fake_fetch)
    SESSION.begin_run("counts_test")
    return fetched


def _listed(geo) -> dict:
    return counts.list_geo_count_sources("GSE1234")


class TestListSources:
    def test_lists_files_samples_and_ncbi(self, geo):
        r = _listed(geo)
        assert r["n_samples"] == 4
        assert r["organisms"] == ["Homo sapiens"]
        assert r["samples"]["GSM1"]["characteristics"] == {"tissue": "cornea", "infection": "mock"}
        files = {f["name"]: f for f in r["files"]}
        assert "filelist.txt" not in files
        assert files["GSE1234_counts.txt.gz"]["supported"]
        assert not files["GSE1234_FPKM.xlsx"]["supported"]
        assert files["GSE1234_RAW.tar"]["contents"][0]["name"] == "GSM1_counts.txt.gz"
        assert r["ncbi_generated_available"]
        assert "GSE1234_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz" not in files  # raw only
        assert all("url" not in f for f in r["files"])  # URLs stay on disk

    def test_series_dir(self):
        assert counts._series_dir("GSE164073").endswith("/GSE164nnn/GSE164073")
        assert counts._series_dir("GSE12").endswith("/GSEnnn/GSE12")

    def test_rejects_non_gse(self, geo):
        assert counts.list_geo_count_sources("SRP123")["error"] == "bad_accession"


class TestPreview:
    def test_r_style_header(self, geo):
        _listed(geo)
        p = counts.preview_geo_file("GSE1234_counts.txt.gz")
        assert p["format"]["separator"] == "whitespace"
        assert p["format"]["unnamed_first_column"]
        assert p["columns"][0]["name"] == "row_names"
        assert p["columns"][1] == {"name": "S_mock_1", "numeric": True, "integer": True}

    def test_bom_is_stripped(self, geo):
        _listed(geo)
        p = counts.preview_geo_file("GSE1234_bom.csv.gz")
        assert p["columns"][0]["name"] == "Gene"

    def test_featurecounts_annotation_and_title_matches(self, geo):
        _listed(geo)
        p = counts.preview_geo_file("GSE1234_fc.txt")
        cols = {c["name"]: c for c in p["columns"]}
        assert p["format"]["comment_lines_skipped"] == 1
        assert cols["Length"]["looks_like_annotation"]
        assert cols["Mock rep1"]["matches_gsm"] == "GSM1"

    def test_unknown_and_unsupported(self, geo):
        _listed(geo)
        assert counts.preview_geo_file("nope.txt")["error"] == "unknown_file"
        assert counts.preview_geo_file("GSE1234_RAW.tar")["error"] == "unsupported"

    def test_requires_listing(self, geo):
        assert counts.preview_geo_file("GSE1234_counts.txt.gz")["error"] == "no_sources"

    def test_captcha_page_is_refused(self, geo, monkeypatch):
        _listed(geo)
        monkeypatch.setattr(counts, "_fetch", lambda url, timeout=60: b"<!doctype html><html>captcha")
        r = counts.preview_geo_file("GSE1234_bom.csv.gz")
        assert r["error"] == "fetch_failed"
        assert "rate limiting" in r["message"]


class TestFetchCounts:
    MAP = {"S_mock_1": "GSM1", "S_mock_2": "GSM2", "S_cov_1": "GSM3", "S_cov_2": "GSM4"}

    def test_writes_counts_renamed_to_titles(self, geo):
        _listed(geo)
        r = counts.fetch_geo_counts("GSE1234_counts.txt.gz", "row_names", self.MAP)
        assert r["value_type"] == "raw_integer_counts"
        assert r["n_duplicates_summed"] == 1
        df = pd.read_csv(SESSION.paths.counts_matrix, sep="\t", index_col=0)
        assert list(df.columns) == ["gene_name", "Mock_rep1", "Mock_rep2", "CoV2_rep1", "CoV2_rep2"]
        assert df.loc["GAPDH", "Mock_rep1"] == 101  # duplicates summed
        meta = json.loads(SESSION.paths.counts_metadata.read_text())
        assert meta["samples"][0] == {"sample": "Mock_rep1", "gsm": "GSM1", "title": "Mock rep1", "file_column": "S_mock_1"}
        assert meta["md5"]

    def test_featurecounts_ensembl_gets_symbols(self, geo):
        _listed(geo)
        smap = {"Mock rep1": "GSM1", "Mock rep2": "GSM2", "CoV2 rep1": "GSM3", "CoV2 rep2": "GSM4"}
        r = counts.fetch_geo_counts("GSE1234_fc.txt", "Geneid", smap)
        assert r["gene_id_type"] == "ensembl"
        assert set(r["dropped_columns"]) == {"Chr", "Start", "End", "Strand", "Length"}
        df = pd.read_csv(SESSION.paths.counts_matrix, sep="\t", index_col=0)
        assert df.loc["ENSG00000111640", "gene_name"] == "GAPDH"  # version stripped, symbol mapped

    def test_ncbi_entrez_gets_symbols(self, geo):
        _listed(geo)
        r = counts.fetch_geo_counts(
            "GSE1234_raw_counts_GRCh38.p13_NCBI.tsv.gz", "GeneID", {g: g for g in ("GSM1", "GSM2", "GSM3", "GSM4")},
        )
        assert r["gene_id_type"] == "entrez"
        df = pd.read_csv(SESSION.paths.counts_matrix, sep="\t", index_col=0)
        assert df.loc[2597, "gene_name"] == "GAPDH"

    def test_partial_map_reports_unmapped(self, geo):
        _listed(geo)
        r = counts.fetch_geo_counts("GSE1234_counts.txt.gz", "row_names", {"S_mock_1": "GSM1", "S_cov_1": "GSM3"})
        assert {u["gsm"] for u in r["unmapped_gsms"]} == {"GSM2", "GSM4"}
        assert set(r["dropped_columns"]) == {"S_mock_2", "S_cov_2"}

    @pytest.mark.parametrize("gene_col,smap,error", [
        ("nope", MAP, "bad_gene_id_column"),
        ("row_names", {"missing": "GSM1"}, "unknown_columns"),
        ("row_names", {"S_mock_1": "GSM9"}, "unknown_gsm"),
        ("row_names", {"S_mock_1": "GSM1", "S_mock_2": "GSM1"}, "duplicate_gsm"),
        ("row_names", {"row_names": "GSM1"}, "non_numeric"),
        ("row_names", {}, "empty_sample_map"),
    ])
    def test_validation(self, geo, gene_col, smap, error):
        _listed(geo)
        assert counts.fetch_geo_counts("GSE1234_counts.txt.gz", gene_col, smap)["error"] == error

    def test_value_types(self):
        assert counts._value_type(pd.DataFrame({"a": [1.0, 2.0], "b": [3.0, 4.0]}))[0] == "raw_integer_counts"
        assert counts._value_type(pd.DataFrame({"a": [1.5, 7.2], "b": [3.1, 4.4]}))[0] == "log_scale"
        tpm = pd.DataFrame({"a": [500000.5, 499999.5], "b": [250000.25, 749999.75]})
        assert counts._value_type(tpm)[0] == "tpm_like"


class TestSaveDesign:
    def _fetch(self, geo):
        _listed(geo)
        counts.fetch_geo_counts("GSE1234_counts.txt.gz", "row_names", TestFetchCounts.MAP)

    def test_condition_field_with_covariates(self, geo):
        self._fetch(geo)
        r = counts.save_geo_design(condition_field="infection")
        assert r["conditions"] == {"mock": 2, "SARS-CoV-2": 2}
        assert r["covariates"] == ["tissue"]
        assert r["low_replication"] == ["mock", "SARS-CoV-2"]
        design = pd.read_csv(SESSION.paths.design)
        assert list(design.columns) == ["sample", "condition", "gsm", "title", "tissue"]
        assert design.loc[0, "sample"] == "Mock_rep1"

    def test_explicit_conditions(self, geo):
        self._fetch(geo)
        r = counts.save_geo_design(conditions={"GSM1": "ctrl", "GSM2": "ctrl", "GSM3": "cov", "GSM4": "cov"})
        assert r["conditions"] == {"ctrl": 2, "cov": 2}
        assert set(r["covariates"]) == {"tissue", "infection"}

    @pytest.mark.parametrize("kwargs,error", [
        ({}, "bad_arguments"),
        ({"condition_field": "x", "conditions": {}}, "bad_arguments"),
        ({"condition_field": "genotype"}, "bad_field"),
        ({"conditions": {"GSM1": "a"}}, "missing_conditions"),
        ({"conditions": {"GSM1": "a", "GSM2": "a", "GSM3": "a", "GSM4": "a"}}, "single_condition"),
    ])
    def test_validation(self, geo, kwargs, error):
        self._fetch(geo)
        assert counts.save_geo_design(**kwargs)["error"] == error

    def test_requires_counts(self, geo):
        _listed(geo)
        assert counts.save_geo_design(condition_field="infection")["error"] == "no_counts"

    def test_design_loads_in_analysis(self, geo):
        from agents.analysis.tools import load_counts
        self._fetch(geo)
        counts.save_geo_design(condition_field="infection")
        r = load_counts(str(SESSION.paths.counts_matrix), str(SESSION.paths.design))
        assert r["n_samples"] == 4
        assert r["has_gene_names"]
        assert r["conditions"] == ["SARS-CoV-2", "mock"]

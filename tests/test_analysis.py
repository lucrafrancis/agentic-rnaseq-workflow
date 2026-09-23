"""Tests for the analysis agent tools.

All tests run offline — no API calls, no network. Tool functions are called directly
with synthetic data from conftest fixtures.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from agents.analysis.tools import (
    compute_qc,
    fetch_abstract,
    fetch_geo_metadata,
    filter_low_counts,
    generate_figures,
    generate_report,
    get_top_genes,
    inspect_counts,
    load_counts,
    run_deseq2,
    run_enrichment,
    scan_results,
    set_design,
    summarize_findings,
    write_report,
)
from core.session import SESSION


class TestFetchGeoMetadata:
    def test_rejects_non_gse(self):
        result = fetch_geo_metadata("SRP123456")
        assert result["error"] == "bad_accession"

    def test_parses_response(self, monkeypatch: pytest.MonkeyPatch):
        search_json = json.dumps({
            "esearchresult": {"idlist": ["200245856"]}
        })
        summary_json = json.dumps({
            "result": {
                "200245856": {
                    "title": "VPA treatment study",
                    "summary": "RNA-seq of VPA-treated cells",
                    "taxon": "Homo sapiens",
                    "n_samples": 6,
                    "samples": [
                        {"accession": "GSM1", "title": "CTRL Rep1"},
                        {"accession": "GSM2", "title": "VPA Rep1"},
                    ],
                    "pubmedids": ["12345678"],
                    "gpl": "28038",
                }
            }
        })
        call_count = [0]

        def mock_get(url):
            call_count[0] += 1
            if "esearch" in url:
                return search_json
            return summary_json

        monkeypatch.setattr("agents.analysis.tools._ncbi_get", mock_get)
        result = fetch_geo_metadata("GSE245856")
        assert result["accession"] == "GSE245856"
        assert result["organism"] == "Homo sapiens"
        assert result["title"] == "VPA treatment study"
        assert "12345678" in result["pubmed_ids"]
        assert len(result["samples"]) == 2

    def test_handles_network_error(self, monkeypatch: pytest.MonkeyPatch):
        import urllib.error
        def mock_get(url):
            raise urllib.error.URLError("network down")
        monkeypatch.setattr("agents.analysis.tools._ncbi_get", mock_get)
        result = fetch_geo_metadata("GSE245856")
        assert result["error"] == "fetch_failed"


class TestFetchAbstract:
    def test_rejects_non_numeric(self):
        result = fetch_abstract("not_a_number")
        assert result["error"] == "bad_pmid"

    def test_parses_xml(self, monkeypatch: pytest.MonkeyPatch):
        xml_response = """<?xml version="1.0"?>
        <PubmedArticleSet>
        <PubmedArticle><MedlineCitation>
        <Article>
            <ArticleTitle>VPA induces neuronal genes</ArticleTitle>
            <Abstract><AbstractText>VPA is an HDAC inhibitor that affects gene expression.</AbstractText></Abstract>
            <AuthorList>
                <Author><LastName>Smith</LastName><Initials>AB</Initials></Author>
                <Author><LastName>Jones</LastName><Initials>CD</Initials></Author>
            </AuthorList>
            <Journal>
                <Title>Nature Genetics</Title>
                <JournalIssue><PubDate><Year>2024</Year></PubDate></JournalIssue>
            </Journal>
        </Article>
        </MedlineCitation></PubmedArticle>
        </PubmedArticleSet>"""

        monkeypatch.setattr("agents.analysis.tools._ncbi_get", lambda url: xml_response)
        result = fetch_abstract("12345678")
        assert result["pmid"] == "12345678"
        assert "VPA induces" in result["title"]
        assert result["authors"] == ["Smith AB", "Jones CD"]
        assert result["journal"] == "Nature Genetics"
        assert result["year"] == "2024"
        assert "HDAC inhibitor" in result["abstract"]

    def test_handles_network_error(self, monkeypatch: pytest.MonkeyPatch):
        import urllib.error
        def mock_get(url):
            raise urllib.error.URLError("timeout")
        monkeypatch.setattr("agents.analysis.tools._ncbi_get", mock_get)
        result = fetch_abstract("12345678")
        assert result["error"] == "fetch_failed"


class TestScanResults:
    def test_finds_count_matrix(self, nfcore_results_dir: Path):
        result = scan_results(str(nfcore_results_dir))
        assert result["counts_found"] is True
        assert "salmon.merged.gene_counts.tsv" in result["counts_path"]

    def test_finds_tpm(self, nfcore_results_dir: Path):
        result = scan_results(str(nfcore_results_dir))
        assert result["tpm_found"] is True

    def test_finds_multiqc(self, nfcore_results_dir: Path):
        result = scan_results(str(nfcore_results_dir))
        assert result["multiqc_found"] is True

    def test_returns_error_for_bad_dir(self):
        result = scan_results("/nonexistent/path")
        assert "error" in result

    def test_detects_design_csv(self, nfcore_results_dir: Path):
        design_path = nfcore_results_dir.parent / "design.csv"
        design_path.write_text("sample,condition\nWT_REP1,WT\n")
        result = scan_results(str(nfcore_results_dir))
        assert result["design_found"] is True


class TestLoadCounts:
    def test_loads_matrix(self, count_matrix_tsv: Path):
        result = load_counts(str(count_matrix_tsv))
        assert "error" not in result
        assert result["n_genes"] == 20
        assert result["n_samples"] == 4
        assert "WT_REP1" in result["samples"]
        assert result["has_gene_names"] is True

    def test_loads_with_design(self, count_matrix_tsv: Path, design_csv: Path):
        result = load_counts(str(count_matrix_tsv), design_path=str(design_csv))
        assert result["design_loaded"] is True
        assert set(result["conditions"]) == {"KO", "WT"}

    def test_returns_error_for_missing_file(self):
        result = load_counts("/nonexistent/file.tsv")
        assert result["error"] == "not_a_file"

    def test_library_sizes_are_positive(self, count_matrix_tsv: Path):
        result = load_counts(str(count_matrix_tsv))
        for size in result["library_sizes"].values():
            assert size > 0

    def test_summary_is_json_serializable(self, count_matrix_tsv: Path):
        result = load_counts(str(count_matrix_tsv))
        json.dumps(result)


class TestSetDesign:
    def test_requires_loaded_counts(self):
        result = set_design([{"sample": "A", "condition": "X"}])
        assert result["error"] == "counts_not_loaded"

    def test_sets_design(self, count_matrix_tsv: Path):
        load_counts(str(count_matrix_tsv))
        SESSION.begin_run("test")
        result = set_design([
            {"sample": "WT_REP1", "condition": "WT"},
            {"sample": "WT_REP2", "condition": "WT"},
            {"sample": "KO_REP1", "condition": "KO"},
            {"sample": "KO_REP2", "condition": "KO"},
        ])
        assert "error" not in result
        assert set(result["conditions"]) == {"KO", "WT"}

    def test_rejects_missing_samples(self, count_matrix_tsv: Path):
        load_counts(str(count_matrix_tsv))
        result = set_design([{"sample": "NONEXISTENT", "condition": "X"}])
        assert result["error"] == "samples_missing_from_design"

    def test_rejects_empty(self, count_matrix_tsv: Path):
        load_counts(str(count_matrix_tsv))
        result = set_design([])
        assert result["error"] == "empty_design"


class TestComputeQC:
    def test_requires_loaded_counts(self):
        result = compute_qc()
        assert result["error"] == "counts_not_loaded"

    def test_computes_metrics(self, count_matrix_tsv: Path):
        load_counts(str(count_matrix_tsv))
        result = compute_qc()
        assert "error" not in result
        assert len(result["library_sizes"]) == 4
        assert len(result["pca"]) == 4
        assert len(result["variance_explained"]) == 2

    def test_parses_multiqc(self, nfcore_results_dir: Path):
        scan_results(str(nfcore_results_dir))
        counts_path = nfcore_results_dir / "star_salmon" / "salmon.merged.gene_counts.tsv"
        load_counts(str(counts_path))
        result = compute_qc()
        assert result["multiqc_summary"] is not None

    def test_summary_is_json_serializable(self, count_matrix_tsv: Path):
        load_counts(str(count_matrix_tsv))
        result = compute_qc()
        json.dumps(result)


class TestFilterLowCounts:
    def test_requires_loaded_counts(self):
        result = filter_low_counts()
        assert result["error"] == "counts_not_loaded"

    def test_filters_genes(self, count_matrix_tsv: Path):
        load_counts(str(count_matrix_tsv))
        result = filter_low_counts(min_count=10, min_samples=2)
        assert "error" not in result
        assert result["genes_before"] == 20
        assert result["genes_after"] <= 20
        assert result["genes_removed"] == result["genes_before"] - result["genes_after"]


class TestRunDeseq2:
    def test_requires_counts(self):
        result = run_deseq2(["condition", "KO", "WT"])
        assert result["error"] == "counts_not_loaded"

    def test_requires_design(self, count_matrix_tsv: Path):
        load_counts(str(count_matrix_tsv))
        result = run_deseq2(["condition", "KO", "WT"])
        assert result["error"] == "design_not_set"

    def test_rejects_bad_contrast(self, count_matrix_tsv: Path, design_csv: Path):
        load_counts(str(count_matrix_tsv), design_path=str(design_csv))
        result = run_deseq2(["condition", "KO"])
        assert result["error"] == "bad_contrast"

    def test_runs_de(self, count_matrix_tsv: Path, design_csv: Path):
        SESSION.begin_run("test")
        load_counts(str(count_matrix_tsv), design_path=str(design_csv))
        filter_low_counts(min_count=1, min_samples=1)
        result = run_deseq2(["condition", "KO", "WT"])
        assert "error" not in result
        assert "n_tested" in result
        assert "n_significant" in result
        assert "top_up" in result
        assert "top_down" in result

    def test_summary_is_json_serializable(self, count_matrix_tsv: Path, design_csv: Path):
        SESSION.begin_run("test")
        load_counts(str(count_matrix_tsv), design_path=str(design_csv))
        filter_low_counts(min_count=1, min_samples=1)
        result = run_deseq2(["condition", "KO", "WT"])
        json.dumps(result)

    def _design_with(self, tmp_path: Path, extra: dict[str, list[str]]) -> Path:
        path = tmp_path / "design_cov.csv"
        rows = {"sample": ["WT_REP1", "WT_REP2", "KO_REP1", "KO_REP2"], "condition": ["WT", "WT", "KO", "KO"], **extra}
        import pandas as pd
        pd.DataFrame(rows).to_csv(path, index=False)
        return path

    def test_paired_covariate(self, tmp_path: Path, count_matrix_tsv: Path):
        SESSION.begin_run("test")
        design = self._design_with(tmp_path, {"donor": ["A", "B", "A", "B"]})
        load_counts(str(count_matrix_tsv), design_path=str(design))
        filter_low_counts(min_count=1, min_samples=1)
        result = run_deseq2(["condition", "KO", "WT"], covariates=["donor"])
        assert "error" not in result
        assert result["design"] == "~donor + condition"
        assert summarize_findings()["de_summary"]["design_formula"] == "~donor + condition"

    @pytest.mark.parametrize("extra,contrast,covariates,error", [
        ({"batch": ["1", "1", "2", "2"]}, ["condition", "KO", "WT"], ["batch"], "confounded_design"),
        ({"batch": ["1", "1", "1", "1"]}, ["condition", "KO", "WT"], ["batch"], "constant_covariate"),
        ({}, ["condition", "KO", "WT"], ["donor"], "bad_factor"),
        ({}, ["condition", "KO", "Control"], None, "bad_level"),
        ({}, ["condition", "KO", "WT"], ["condition"], "bad_covariates"),
        ({"cell line": ["a", "b", "a", "b"]}, ["condition", "KO", "WT"], ["cell line"], "bad_column_name"),
    ])
    def test_design_validation(self, tmp_path, count_matrix_tsv, extra, contrast, covariates, error):
        design = self._design_with(tmp_path, extra)
        load_counts(str(count_matrix_tsv), design_path=str(design))
        assert run_deseq2(contrast, covariates=covariates)["error"] == error

    def test_load_counts_reports_design_columns(self, tmp_path, count_matrix_tsv):
        design = self._design_with(tmp_path, {"donor": ["A", "B", "A", "B"]})
        result = load_counts(str(count_matrix_tsv), design_path=str(design))
        assert result["design_columns"] == ["condition", "donor"]


class TestGetTopGenes:
    def test_requires_deseq_results(self):
        result = get_top_genes()
        assert result["error"] == "no_deseq_results"

    def test_returns_genes(self, count_matrix_tsv: Path, design_csv: Path):
        SESSION.begin_run("test")
        load_counts(str(count_matrix_tsv), design_path=str(design_csv))
        filter_low_counts(min_count=1, min_samples=1)
        run_deseq2(["condition", "KO", "WT"])
        result = get_top_genes(n=5)
        assert "error" not in result
        assert len(result["genes"]) <= 5
        for gene in result["genes"]:
            assert "gene_id" in gene
            assert "log2FC" in gene
            assert "padj" in gene


class TestRunEnrichment:
    def test_rejects_empty_list(self):
        result = run_enrichment([])
        assert result["error"] == "empty_gene_list"

    def test_runs_with_mock(self, monkeypatch: pytest.MonkeyPatch):
        import pandas as pd
        import gseapy

        mock_df = pd.DataFrame({
            "Gene_set": ["GO_Biological_Process_2023"],
            "Term": ["regulation of transcription (GO:0006355)"],
            "P-value": [0.001],
            "Adjusted P-value": [0.01],
            "Overlap": ["5/100"],
            "Genes": ["Gene1;Gene2;Gene3"],
        })

        class MockResult:
            results = mock_df

        monkeypatch.setattr(gseapy, "enrichr", lambda **kw: MockResult())
        result = run_enrichment(["Gene1", "Gene2", "Gene3"])
        assert "error" not in result
        assert "enrichment" in result

    def test_label_accumulates(self, monkeypatch: pytest.MonkeyPatch):
        import pandas as pd
        import gseapy

        mock_df = pd.DataFrame({
            "Gene_set": ["GO_Biological_Process_2023"],
            "Term": ["some term"],
            "P-value": [0.001],
            "Adjusted P-value": [0.01],
            "Overlap": ["3/50"],
            "Genes": ["A;B;C"],
        })

        class MockResult:
            results = mock_df

        monkeypatch.setattr(gseapy, "enrichr", lambda **kw: MockResult())
        run_enrichment(["Gene1"], label="upregulated")
        run_enrichment(["Gene2"], label="downregulated")
        assert "upregulated" in SESSION.enrichment_results
        assert "downregulated" in SESSION.enrichment_results
        assert "results" in SESSION.enrichment_results["upregulated"]


class TestSummarizeFindings:
    def test_returns_empty_when_nothing_loaded(self):
        result = summarize_findings()
        assert isinstance(result, dict)
        assert "software_versions" in result
        assert "python" in result["software_versions"]
        assert "data_source" in result

    def test_includes_de_summary(self, count_matrix_tsv: Path, design_csv: Path):
        SESSION.begin_run("test")
        load_counts(str(count_matrix_tsv), design_path=str(design_csv))
        filter_low_counts(min_count=1, min_samples=1)
        run_deseq2(["condition", "KO", "WT"])
        result = summarize_findings()
        assert "de_summary" in result
        assert "n_significant" in result["de_summary"]
        assert "numpy" in result["software_versions"]
        assert "pandas" in result["software_versions"]


def _strong_de_dataset(tmp_path: Path, n_up: int = 40, n_down: int = 30) -> tuple[Path, Path]:
    """3 vs 3 counts with many clearly up/down genes, plus identifier and covariate columns."""
    import numpy as np
    import pandas as pd
    rng = np.random.default_rng(1)
    n_genes = 300
    base = rng.integers(200, 2000, size=n_genes).astype(float)
    fold = np.ones(n_genes)
    fold[:n_up] = 8.0
    fold[n_up:n_up + n_down] = 0.125
    samples = ["WT_1", "WT_2", "WT_3", "KO_1", "KO_2", "KO_3"]
    data = {s: rng.poisson(base * (fold if s.startswith("KO") else 1)) for s in samples}
    counts = pd.DataFrame(data, index=[f"G{i:03d}" for i in range(n_genes)])
    counts.insert(0, "gene_name", [f"Gene{i}" for i in range(n_genes)])
    counts_path = tmp_path / "counts.tsv"
    counts.rename_axis("gene_id").to_csv(counts_path, sep="\t")
    design_path = tmp_path / "design.csv"
    pd.DataFrame({
        "sample": samples, "condition": ["WT"] * 3 + ["KO"] * 3,
        "gsm": [f"GSM{i}" for i in range(6)], "title": samples,
        "batch": ["a", "b", "c", "a", "b", "c"], "cell_line": ["X"] * 6,
    }).to_csv(design_path, index=False)
    return counts_path, design_path


def _analysed(tmp_path: Path, **kwargs):
    SESSION.begin_run("test")
    counts, design = _strong_de_dataset(tmp_path, **kwargs)
    load_counts(str(counts), design_path=str(design))
    filter_low_counts(min_count=10, min_samples=2)
    run_deseq2(["condition", "KO", "WT"])


class TestGenerateFigures:
    def test_requires_counts(self):
        assert generate_figures()["error"] == "counts_not_loaded"

    def test_figures_have_paths_and_captions(self, tmp_path: Path):
        _analysed(tmp_path)
        figs = generate_figures()["figures"]
        for name in ("library_sizes", "pca", "sample_correlation", "pc_association", "volcano", "ma_plot", "de_heatmap"):
            assert f"figures/{name}.png" in figs
        assert all(caption for caption in figs.values())
        assert all((SESSION.paths.analysis_dir / p).is_file() for p in figs)

    def test_heatmap_is_25_up_and_25_down(self, tmp_path: Path):
        _analysed(tmp_path)
        caption = generate_figures()["figures"]["figures/de_heatmap.png"]
        assert caption.startswith("Top 25 up-regulated (above the line) and 25 down-regulated")

    def test_heatmap_uses_fewer_when_not_enough(self, tmp_path: Path):
        _analysed(tmp_path, n_up=40, n_down=10)
        caption = generate_figures()["figures"]["figures/de_heatmap.png"]
        assert "Top 25 up-regulated" in caption and "and 10 down-regulated" in caption

    def test_pca_skips_identifier_and_constant_columns(self, tmp_path: Path):
        _analysed(tmp_path)
        result = generate_figures()
        assert result["pca_color_options"] == ["batch"]
        assert "figures/pca_batch.png" in result["figures"]
        assert not any(p in result["figures"] for p in ("figures/pca_gsm.png", "figures/pca_title.png",
                                                         "figures/pca_cell_line.png"))
        assert "gsm" not in result["figures"]["figures/pc_association.png"]

    def test_pca_color_by_choice(self, tmp_path: Path):
        _analysed(tmp_path)
        assert "figures/pca_batch.png" not in generate_figures(pca_color_by=[])["figures"]
        assert generate_figures(pca_color_by=["gsm"])["error"] == "bad_pca_columns"


class TestWriteReport:
    def test_requires_figures(self, tmp_path: Path):
        _analysed(tmp_path)
        assert write_report("# R")["error"] == "no_figures"

    def test_captions_and_lossless_path_fix(self, tmp_path: Path):
        _analysed(tmp_path)
        generate_figures()
        md = "# R\n\n![Volcano](volcano.png)\n\n![Heat](figures/de_heatmap.png)\n\nText."
        result = write_report(md)
        content = Path(result["report_path"]).read_text()
        assert "![Volcano](figures/volcano.png)" in content  # prefix added
        assert "*File: `figures/volcano.png` — Volcano plot:" in content
        assert "*File: `figures/de_heatmap.png` — Top 25 up-regulated" in content
        assert result["figures_referenced"] == 2
        assert "figures/pca.png" in result["unreferenced_figures"]

    def test_unknown_figure_rejected_nothing_written(self, tmp_path: Path):
        _analysed(tmp_path)
        generate_figures()
        result = write_report("# R\n\n![Made up](figures/invented_plot.png)")
        assert result["error"] == "unknown_figures"
        assert "invented_plot.png" in result["message"]
        assert not SESSION.paths.analysis_report.exists()

    def test_llm_text_kept_and_one_standard_disclaimer_appended(self, tmp_path: Path):
        _analysed(tmp_path)
        generate_figures()
        md = "# R\n\nDisclaimer: only 3 replicates — exploratory.\n\n---\n"
        content = Path(write_report(md)["report_path"]).read_text()
        assert "Disclaimer: only 3 replicates — exploratory." in content  # never removed
        assert content.count("**Disclaimer:** This report was generated by an AI system.") == 1
        assert "\n---\n\n---" not in content

    def test_generate_report_back_compat(self, tmp_path: Path):
        _analysed(tmp_path)
        result = generate_report("# R\n\n![PCA](figures/pca.png)")
        assert Path(result["report_path"]).is_file()


class TestProvenance:
    def test_geo_counts(self, tmp_path: Path):
        SESSION.begin_run("test")
        counts, design = _strong_de_dataset(tmp_path)
        paths = SESSION.require_paths()
        paths.counts_matrix.write_text(counts.read_text())
        paths.counts_metadata.write_text(json.dumps({
            "accession": "GSE1", "source": "author", "filename": "GSE1_counts.txt.gz", "md5": "abc",
            "value_type": "raw_integer_counts", "gene_id_type": "symbol", "n_duplicates_summed": 0,
        }))
        load_counts(str(paths.counts_matrix), design_path=str(design))
        filter_low_counts(min_count=10, min_samples=2)
        result = summarize_findings()
        assert result["data_source"] == "GEO count matrix"
        assert result["data_provenance"]["source"] == "authors' supplementary file"
        assert result["data_provenance"]["value_type"] == "raw_integer_counts"
        assert result["filtering"]["min_count"] == 10

    def test_nfcore_salmon_only(self, nfcore_results_dir: Path):
        import shutil
        SESSION.begin_run("test")
        shutil.move(nfcore_results_dir / "star_salmon", nfcore_results_dir / "salmon")
        scan_results(str(nfcore_results_dir))
        load_counts(str(nfcore_results_dir / "salmon" / "salmon.merged.gene_counts.tsv"))
        result = summarize_findings()
        assert result["data_source"] == "nf-core/rnaseq"
        assert "alignment skipped" in result["data_provenance"]["quantification"]

    def test_user_provided(self, count_matrix_tsv: Path):
        SESSION.begin_run("test")
        load_counts(str(count_matrix_tsv))
        assert summarize_findings()["data_source"] == "user-provided"


class TestScanResultsLooseFiles:
    """scan_results should detect user-provided loose count files."""

    def test_finds_loose_counts_csv(self, tmp_path: Path):
        d = tmp_path / "data"
        d.mkdir()
        (d / "counts.csv").write_text("gene,sampleA,sampleB\nGENE1,100,200\n")
        result = scan_results(str(d))
        assert result["counts_found"] is True
        assert result["source"] == "user-provided"
        assert "counts.csv" in result["counts_path"]

    def test_finds_design_file(self, tmp_path: Path):
        d = tmp_path / "data"
        d.mkdir()
        (d / "counts.csv").write_text("gene,sampleA,sampleB\nGENE1,100,200\n")
        (d / "design.csv").write_text("sample,condition\nsampleA,ctrl\nsampleB,treat\n")
        result = scan_results(str(d))
        assert result["design_found"] is True

    def test_lists_tabular_files(self, tmp_path: Path):
        d = tmp_path / "data"
        d.mkdir()
        (d / "counts.tsv").write_text("gene\tsA\nG1\t10\n")
        (d / "metadata.csv").write_text("sample,condition\nsA,ctrl\n")
        (d / "readme.txt").write_text("notes")
        result = scan_results(str(d))
        assert len(result["tabular_files"]) == 3

    def test_nfcore_takes_priority(self, nfcore_results_dir: Path):
        result = scan_results(str(nfcore_results_dir))
        assert result["source"] == "nf-core"


class TestLoadCountsFormats:
    """load_counts should handle multiple input formats."""

    def test_generic_csv(self, tmp_path: Path):
        path = tmp_path / "counts.csv"
        path.write_text("gene,sampleA,sampleB\nGENE1,100,200\nGENE2,50,75\n")
        result = load_counts(str(path))
        assert "error" not in result
        assert result["detected_format"] == "generic"
        assert result["n_genes"] == 2
        assert result["n_samples"] == 2

    def test_featurecounts_format(self, tmp_path: Path):
        path = tmp_path / "counts.tsv"
        path.write_text(
            "Geneid\tChr\tStart\tEnd\tStrand\tLength\tsA\tsB\n"
            "GENE1\tchr1\t100\t200\t+\t100\t500\t600\n"
            "GENE2\tchr1\t300\t400\t-\t100\t50\t75\n"
        )
        result = load_counts(str(path))
        assert "error" not in result
        assert result["detected_format"] == "featureCounts"
        assert result["n_samples"] == 2
        assert "sA" in result["samples"]

    def test_auto_detects_tsv(self, tmp_path: Path):
        path = tmp_path / "data.tsv"
        path.write_text("gene\tsA\tsB\nG1\t10\t20\nG2\t30\t40\n")
        result = load_counts(str(path))
        assert "error" not in result
        assert result["n_samples"] == 2

    def test_loads_design_tsv(self, tmp_path: Path):
        counts = tmp_path / "counts.csv"
        counts.write_text("gene,sA,sB\nG1,10,20\n")
        design = tmp_path / "design.tsv"
        design.write_text("sample\tcondition\nsA\tctrl\nsB\ttreat\n")
        result = load_counts(str(counts), design_path=str(design))
        assert result["design_loaded"] is True
        assert set(result["conditions"]) == {"ctrl", "treat"}

    def test_rejects_non_numeric(self, tmp_path: Path):
        path = tmp_path / "bad.csv"
        path.write_text("gene,sA,sB\nG1,hello,world\n")
        result = load_counts(str(path))
        assert result["error"] == "non_numeric_columns"


class TestInspectCounts:
    def test_requires_loaded_counts(self):
        result = inspect_counts()
        assert result["error"] == "counts_not_loaded"

    def test_raw_counts_hints(self, count_matrix_tsv: Path):
        load_counts(str(count_matrix_tsv))
        result = inspect_counts()
        assert "error" not in result
        assert result["global"]["fraction_non_integer"] < 0.01
        assert any("raw counts" in h for h in result["hints"])

    def test_log_transformed_hints(self, tmp_path: Path):
        import numpy as np
        path = tmp_path / "log_counts.csv"
        np.random.seed(42)
        vals = np.random.uniform(0, 15, (100, 4))
        lines = ["gene,sA,sB,sC,sD"]
        for i, row in enumerate(vals):
            lines.append(f"G{i},{row[0]:.4f},{row[1]:.4f},{row[2]:.4f},{row[3]:.4f}")
        path.write_text("\n".join(lines) + "\n")
        load_counts(str(path))
        result = inspect_counts()
        assert result["global"]["fraction_non_integer"] > 0.5
        assert any("log-transformed" in h for h in result["hints"])

    def test_per_sample_stats(self, count_matrix_tsv: Path):
        load_counts(str(count_matrix_tsv))
        result = inspect_counts()
        assert "WT_REP1" in result["per_sample"]
        assert "mean" in result["per_sample"]["WT_REP1"]

    def test_json_serializable(self, count_matrix_tsv: Path):
        load_counts(str(count_matrix_tsv))
        result = inspect_counts()
        json.dumps(result)


class TestSchemas:
    def test_every_schema_has_a_callable(self):
        from agents.analysis.schemas import TOOL_FUNCTIONS, TOOL_SCHEMAS
        for schema in TOOL_SCHEMAS:
            name = schema["name"]
            assert name in TOOL_FUNCTIONS, f"{name} exposed to Claude but not dispatchable"
            assert callable(TOOL_FUNCTIONS[name])

    def test_schemas_are_well_formed(self):
        from agents.analysis.schemas import TOOL_SCHEMAS
        for schema in TOOL_SCHEMAS:
            assert {"name", "description", "input_schema"} <= schema.keys()
            assert schema["input_schema"]["type"] == "object"

    def test_no_duplicate_tool_names(self):
        from agents.analysis.schemas import TOOL_SCHEMAS
        names = [s["name"] for s in TOOL_SCHEMAS]
        assert len(names) == len(set(names))


class TestReplayScript:
    def _log(self, path: Path, entries: list[tuple[str, dict, dict]]) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("a") as f:
            for tool, args, summary in entries:
                f.write(json.dumps({"tool": tool, "args": args, "summary": summary}) + "\n")
        return path

    def test_skips_failed_and_read_only_calls(self, tmp_path: Path):
        from agents.analysis.replay import write_replay_script

        log = self._log(tmp_path / "run" / "analysis" / "tool_calls.jsonl", [
            ("fetch_geo_metadata", {"accession": "GSE1"}, {"title": "x"}),
            ("load_counts", {"counts_path": "missing.tsv"}, {"error": "not_a_file"}),
            ("load_counts", {"counts_path": "counts.tsv"}, {"n_genes": 20}),
            ("get_top_genes", {"n": 10}, {"genes": []}),
            ("run_deseq2", {"contrast": ["condition", "KO", "WT"]}, {"n_up": 5}),
        ])
        script = write_replay_script(log, log.parent / "replay.py")
        text = script.read_text()
        compile(text, str(script), "exec")
        assert "missing.tsv" not in text
        assert "fetch_geo_metadata" not in text
        assert "get_top_genes" not in text
        assert "counts_path='counts.tsv'" in text
        assert "contrast=['condition', 'KO', 'WT']" in text

    def test_start_line_skips_earlier_sessions(self, tmp_path: Path):
        from agents.analysis.replay import write_replay_script

        log = self._log(tmp_path / "run" / "analysis" / "tool_calls.jsonl", [
            ("filter_low_counts", {"min_count": 99}, {"genes_after": 1}),
            ("filter_low_counts", {"min_count": 5}, {"genes_after": 15}),
        ])
        text = write_replay_script(log, log.parent / "replay.py", start_line=1).read_text()
        assert "min_count=5" in text
        assert "min_count=99" not in text

    def test_nothing_to_replay(self, tmp_path: Path):
        from agents.analysis.replay import write_replay_script

        log = self._log(tmp_path / "run" / "analysis" / "tool_calls.jsonl", [
            ("inspect_counts", {}, {"likely_raw_counts": True}),
        ])
        assert write_replay_script(log, log.parent / "replay.py") is None
        assert write_replay_script(tmp_path / "absent.jsonl", tmp_path / "replay.py") is None

    def test_replay_reproduces_de_results(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
        count_matrix_tsv: Path, design_csv: Path,
    ):
        import pandas as pd
        from agents.analysis.replay import write_replay_script
        from core import config

        SESSION.begin_run("original")
        paths = SESSION.require_paths()
        calls = [
            ("load_counts", {"counts_path": str(count_matrix_tsv), "design_path": str(design_csv)}),
            ("filter_low_counts", {"min_count": 10, "min_samples": 2}),
            ("run_deseq2", {"contrast": ["condition", "KO", "WT"]}),
        ]
        tools_mod = __import__("agents.analysis.tools", fromlist=["x"])
        self._log(paths.analysis_tool_log, [
            (name, args, getattr(tools_mod, name)(**args)) for name, args in calls
        ])
        original = pd.read_csv(paths.de_results, index_col=0)

        script = write_replay_script(paths.analysis_tool_log, paths.analysis_dir / "replay.py")
        monkeypatch.chdir(tmp_path)  # the script chdirs to the repo root; restore afterwards
        exec(compile(script.read_text(), str(script), "exec"), {
            "__file__": str(config.ROOT / "runs" / "original" / "analysis" / "replay.py"),
            "__name__": "__main__",
        })

        assert SESSION.require_paths().name.endswith("original_replay")
        replayed = pd.read_csv(SESSION.require_paths().de_results, index_col=0)
        pd.testing.assert_frame_equal(original, replayed)

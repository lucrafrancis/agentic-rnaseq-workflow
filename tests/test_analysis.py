"""Tests for the analysis agent tools.

All tests run offline — no API calls, no network. Tool functions are called directly
with synthetic data from conftest fixtures.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from agents.analysis import tools
from agents.analysis.tools import (
    compute_qc,
    read_multiqc,
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
        assert result["cite_as"] == "{{cite:12345678}}"
        assert SESSION.references["12345678"]["authors"] == ["Smith AB", "Jones CD"]

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
        assert len(result["multiqc_files"]) == 1

    def test_finds_multiqc_in_any_layout(self, nfcore_results_dir: Path):
        # nf-core/rnaseq 3.26 salmon-only writes multiqc/multiqc_report_data/ (no aligner subfolder)
        flat = nfcore_results_dir / "multiqc" / "multiqc_report_data"
        flat.mkdir(parents=True)
        (flat / "multiqc_general_stats.txt").write_text("Sample\tsalmon-percent_mapped\nWT_REP1\t85.0\n")
        result = scan_results(str(nfcore_results_dir))
        assert len(result["multiqc_files"]) == 2

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

    def test_reads_multiqc(self, nfcore_results_dir: Path):
        mqc = nfcore_results_dir / "multiqc" / "star_salmon" / "multiqc_report_data" / "multiqc_general_stats.txt"
        mqc.write_text(mqc.read_text() + "WT_REP1 Read 1\t\t\n")
        scan_results(str(nfcore_results_dir))
        load_counts(str(nfcore_results_dir / "star_salmon" / "salmon.merged.gene_counts.tsv"))
        result = read_multiqc()
        assert result["n_samples_matched"] == 4
        assert result["columns"]["star-mapped_percent"]["KO_REP2"] == 92.0
        facts = summarize_findings()["facts"]
        assert facts["multiqc.star_mapped_percent.min"] == "89.00"
        assert facts["multiqc.star_mapped_percent.max"] == "92.00"

    def test_multiqc_asks_to_choose_between_files(self, nfcore_results_dir: Path):
        flat = nfcore_results_dir / "multiqc" / "multiqc_report_data"
        flat.mkdir(parents=True)
        (flat / "multiqc_general_stats.txt").write_text("Sample\tsalmon-percent_mapped\nWT_REP1\t85.0\n")
        scan_results(str(nfcore_results_dir))
        assert read_multiqc()["error"] == "choose_file"
        assert read_multiqc("/elsewhere/multiqc_general_stats.txt")["error"] == "unknown_file"
        assert read_multiqc(str(flat / "multiqc_general_stats.txt"))["n_samples_matched"] == 1

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


def _mock_enrichr(monkeypatch, captured: list | None = None):
    import pandas as pd
    import gseapy

    def enrichr(**kw):
        if captured is not None:
            captured.append(kw["gene_list"])

        class Result:
            results = pd.DataFrame({
                "Gene_set": ["GO_Biological_Process_2023", "KEGG_2021_Human"],
                "Term": ["regulation of transcription (GO:0006355)", "Pathway X"],
                "P-value": [0.001, 0.002], "Adjusted P-value": [0.01, 0.02],
                "Overlap": ["5/100", "3/40"], "Genes": ["Gene1;Gene2", "Gene3"],
            })
        return Result()

    monkeypatch.setattr(gseapy, "enrichr", enrichr)


class TestRunEnrichment:
    def test_requires_deseq(self):
        assert run_enrichment("up")["error"] == "no_deseq_results"

    @pytest.mark.parametrize("kwargs,error", [
        ({"direction": "sideways"}, "bad_direction"),
        ({"direction": "up", "organism": "zebrafish"}, "bad_organism"),
        ({"direction": "up", "max_genes": 0}, "bad_max_genes"),
    ])
    def test_validation(self, tmp_path, kwargs, error):
        _analysed(tmp_path)
        assert run_enrichment(**kwargs)["error"] == error

    def test_selects_genes_itself(self, tmp_path, monkeypatch):
        from agents.analysis.tools import _ranked
        _analysed(tmp_path)
        captured: list = []
        _mock_enrichr(monkeypatch, captured)
        result = run_enrichment("down")
        expected = [f"Gene{int(g[1:])}" for g in _ranked("down", 0.05, 1.0).index]
        assert captured[0] == expected  # symbols, in rank order, chosen by the tool
        assert result["selection"]["n_passing"] == len(expected)
        assert result["label"] == "downregulated"
        saved = (SESSION.paths.analysis_dir / "enrichment_input_downregulated.txt").read_text().split()
        assert saved == expected
        assert (SESSION.paths.analysis_dir / "enrichment_downregulated.csv").is_file()

    def test_cap_and_thresholds(self, tmp_path, monkeypatch):
        _analysed(tmp_path)
        captured: list = []
        _mock_enrichr(monkeypatch, captured)
        result = run_enrichment("up", max_genes=5)
        assert len(captured[0]) == 5
        assert result["selection"]["capped"] is True
        assert run_enrichment("up", lfc_min=50)["error"] == "no_genes"

    def test_both_directions_accumulate(self, tmp_path, monkeypatch):
        _analysed(tmp_path)
        _mock_enrichr(monkeypatch)
        run_enrichment("up")
        run_enrichment("down")
        assert set(SESSION.enrichment_results) == {"upregulated", "downregulated"}


class TestTableFormatting:
    def test_underflow_padj(self):
        from agents.analysis.tools import _fmt_p
        assert _fmt_p(0.0) == "< 1e-300"
        assert _fmt_p(2.5e-8) == "2.50e-08"
        assert _fmt_p(0.0123) == "0.012"

    def test_gene_id_column_only_when_ids_are_not_symbols(self):
        import pandas as pd
        from agents.analysis.tools import _tables
        SESSION.counts_df = None
        SESSION.deseq_results = pd.DataFrame(
            {"log2FoldChange": [2.0, -2.0], "padj": [0.001, 0.002], "baseMean": [10.0, 20.0]},
            index=["CTGF", "KDR"])
        SESSION.deseq_design, SESSION.deseq_contrast = "~condition", ["condition", "B", "A"]
        assert _tables()["top_up"].splitlines()[0] == "| Gene | log2FC | padj | baseMean |"
        SESSION.gene_names = pd.Series({"CTGF": "CTGF", "KDR": "KDR"})
        SESSION.deseq_results.index = ["ENSG1", "ENSG2"]
        SESSION.gene_names = pd.Series({"ENSG1": "CTGF", "ENSG2": "KDR"})
        assert _tables()["top_up"].splitlines()[0] == "| Gene | Gene ID | log2FC | padj | baseMean |"


class TestRanking:
    def test_ties_broken_by_abs_lfc_then_id(self):
        import pandas as pd
        from agents.analysis.tools import _ranked
        SESSION.deseq_results = pd.DataFrame({
            "log2FoldChange": [2.0, 5.0, 5.0, 3.0, -4.0, 0.5],
            "padj": [0.0, 0.0, 0.0, 1e-5, 0.0, 1e-9],
            "baseMean": [1.0] * 6,
        }, index=["B", "C", "A", "D", "E", "F"])
        assert list(_ranked("up").index) == ["A", "C", "B", "F", "D"]
        assert list(_ranked("up", lfc_min=1.0).index) == ["A", "C", "B", "D"]
        assert list(_ranked("down").index) == ["E"]


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

    def test_condition_accepted_in_pca_color_by(self, tmp_path: Path):
        _analysed(tmp_path)
        result = generate_figures(pca_color_by=["condition", "batch"])
        assert "error" not in result
        assert "figures/pca_batch.png" in result["figures"]
        assert "figures/pca_condition.png" not in result["figures"]  # condition is pca.png

    def test_pca_color_by_choice(self, tmp_path: Path):
        _analysed(tmp_path)
        assert "figures/pca_batch.png" not in generate_figures(pca_color_by=[])["figures"]
        assert generate_figures(pca_color_by=["gsm"])["error"] == "bad_pca_columns"


TABLES = "{{table:de_summary}}\n\n{{table:top_up}}\n\n{{table:top_down}}\n\n{{table:methods}}\n\n"


def _figs(skip: tuple[str, ...] = ()) -> str:
    """Markdown placing every generated figure (write_report requires all of them)."""
    return "".join(f"![{Path(f['path']).stem}]({f['path']})\n\n"
                   for f in SESSION.figures if Path(f["path"]).name not in skip)


class TestWriteReport:
    def test_requires_figures(self, tmp_path: Path):
        _analysed(tmp_path)
        assert write_report("# R")["error"] == "no_figures"

    def test_captions_and_lossless_path_fix(self, tmp_path: Path):
        _analysed(tmp_path)
        generate_figures()
        md = ("# R\n\n![Volcano](volcano.png)\n\n![Heat](figures/de_heatmap.png)\n\nText.\n\n" + TABLES
              + _figs(skip=("volcano.png", "de_heatmap.png")))
        result = write_report(md)
        content = Path(result["report_path"]).read_text()
        assert "![Volcano](figures/volcano.png)" in content  # prefix added
        assert "*File: `figures/volcano.png` — Volcano plot:" in content
        assert "*File: `figures/de_heatmap.png` — Top 25 up-regulated" in content
        assert result["figures_referenced"] == len(SESSION.figures)
        assert result["unreferenced_figures"] == []

    def test_unknown_figure_rejected_nothing_written(self, tmp_path: Path):
        _analysed(tmp_path)
        generate_figures()
        result = write_report("# R\n\n![Made up](figures/invented_plot.png)\n\n" + TABLES)
        assert result["error"] == "report_problems"
        assert "invented_plot.png" in result["message"]
        assert not SESSION.paths.analysis_report.exists()

    def test_llm_text_kept_and_one_standard_disclaimer_appended(self, tmp_path: Path):
        _analysed(tmp_path)
        generate_figures()
        md = "# R\n\n" + TABLES + _figs() + "Disclaimer: only three replicates — exploratory.\n\n---\n"
        content = Path(write_report(md)["report_path"]).read_text()
        assert "Disclaimer: only three replicates — exploratory." in content  # never removed
        assert content.count("**Disclaimer:** This report was generated by an AI system.") == 1
        assert "\n---\n\n---" not in content

    def test_generate_report_back_compat(self, tmp_path: Path):
        _analysed(tmp_path)
        generate_figures()  # only to learn the figure names; generate_report redraws them
        result = generate_report("# R\n\n" + _figs() + TABLES)
        assert Path(result["report_path"]).is_file()


class TestPlaceholders:
    def _ready(self, tmp_path, monkeypatch):
        _analysed(tmp_path)
        _mock_enrichr(monkeypatch)
        run_enrichment("up")
        run_enrichment("down")
        generate_figures()

    def _full(self, body: str) -> str:
        return (body + "\n\n" + TABLES + "{{table:enrichment_up}}\n\n{{table:enrichment_down}}\n\n" + _figs())

    def test_facts_render_from_results(self, tmp_path, monkeypatch):
        self._ready(tmp_path, monkeypatch)
        facts = summarize_findings()["facts"]
        md = self._full("Of {{de.n_tested}} genes, {{de.n_significant}} ({{de.pct_significant}}) were "
                        "significant; design {{de.design}}; filter {{filter.min_count}}; "
                        "PyDESeq2 {{versions.pydeseq2}}; source {{provenance.data_source}}.")
        content = Path(write_report(md)["report_path"]).read_text()
        for key in ("de.n_tested", "de.n_significant", "de.pct_significant", "de.design"):
            assert facts[key] in content
        assert "{{" not in content

    def test_percent_after_percentage_fact_not_doubled(self, tmp_path, monkeypatch):
        self._ready(tmp_path, monkeypatch)
        pct = summarize_findings()["facts"]["de.pct_up"]
        md = self._full("A: {{de.pct_up}}% up. B: {{de.n_up}}% kept.")
        content = Path(write_report(md)["report_path"]).read_text()
        assert f"A: {pct} up." in content and "%%" not in content
        assert "B: " + summarize_findings()["facts"]["de.n_up"] + "% kept." in content  # non-% fact untouched

    def test_full_stop_after_sentence_fact_not_doubled(self, tmp_path, monkeypatch):
        self._ready(tmp_path, monkeypatch)
        summary = summarize_findings()["facts"]["qc.pca_summary"]
        assert summary.endswith(".")
        content = Path(write_report(self._full("PCA: {{qc.pca_summary}}. Next."))["report_path"]).read_text()
        assert f"PCA: {summary} Next." in content

    def test_fact_prefix_is_accepted(self, tmp_path, monkeypatch):
        self._ready(tmp_path, monkeypatch)
        design = summarize_findings()["facts"]["de.design"]
        content = Path(write_report(self._full("Design {{fact.de.design}}."))["report_path"]).read_text()
        assert f"Design {design}." in content

    def test_table_inside_a_sentence_rejected(self, tmp_path, monkeypatch):
        self._ready(tmp_path, monkeypatch)
        result = write_report(self._full("- Genes fell (see {{table:top_down}})."))
        assert "{{table:top_down}} must be on a line of its own" in result["message"]

    def test_version_facts_name_the_software_once(self, tmp_path, monkeypatch):
        self._ready(tmp_path, monkeypatch)
        facts = summarize_findings()["facts"]
        assert facts["versions.pydeseq2"].startswith("PyDESeq2 ")
        version = facts["versions.pydeseq2"].split()[-1]
        md = self._full("A: tested with {{versions.pydeseq2}}. B: PyDESeq2 {{versions.pydeseq2}}. "
                        "C: pydeseq2 ({{versions.pydeseq2}}).")
        content = Path(write_report(md)["report_path"]).read_text()
        assert f"A: tested with PyDESeq2 {version}. B: PyDESeq2 {version}. C: pydeseq2 ({version})." in content

    def test_fact_values_are_correct(self, tmp_path, monkeypatch):
        self._ready(tmp_path, monkeypatch)
        facts = summarize_findings()["facts"]
        res = SESSION.deseq_results.dropna(subset=["padj"])
        sig = res[res["padj"] < 0.05]
        assert facts["de.n_tested"] == f"{len(res):,}"
        assert facts["de.n_up"] == f"{int((sig['log2FoldChange'] > 0).sum()):,}"
        assert facts["de.pct_significant"] == f"{100 * len(sig) / len(res):.1f}%"
        assert facts["de.contrast"] == "KO vs WT"
        assert facts["samples.n_ko"] == "3"
        assert facts["enrichment.up.lfc_min"] == "1.0"

    def test_new_qc_filter_design_and_term_facts(self, tmp_path, monkeypatch):
        from agents.analysis.tools import _pca
        self._ready(tmp_path, monkeypatch)
        facts = summarize_findings()["facts"]
        pca = _pca()
        assert facts["qc.pca_pc1_pct"] == f"{pca['var_exp'][0]:.1%}"
        caption = generate_figures()["figures"]["figures/pca.png"]
        assert facts["qc.pca_pc1_pct"] in caption and facts["qc.pca_pc2_pct"] in caption
        for key in ("qc.genes_detected_min", "qc.genes_detected_max", "qc.library_size_median_millions",
                    "qc.sample_correlation_min", "filter.pct_removed"):
            assert facts[key]
        assert facts["design.cell_line"] == "X"            # constant column -> fact
        assert facts["design.batch_values"] == "a, b, c"    # varying column -> its levels
        assert not any(k.startswith(("design.gsm", "design.title")) for k in facts)
        assert facts["enrichment.up.go_bp.1"] == (
            "regulation of transcription (GO:0006355) (5/100 genes, padj 0.010)")
        assert facts["enrichment.down.kegg.1"].startswith("Pathway X (3/40 genes")

    def test_pct_rounding_is_correct(self):
        from agents.analysis.tools import _fmt_pct
        assert _fmt_pct(13009, 27946) == "46.6%"  # the value Haiku typed as 46.5%

    def test_library_stats_are_pre_filter_everywhere(self, tmp_path):
        SESSION.begin_run("test")
        counts, design = _strong_de_dataset(tmp_path)
        load_counts(str(counts), design_path=str(design))
        compute_qc()
        raw_lib = SESSION.counts_df.sum(axis=0)
        filter_low_counts(min_count=500, min_samples=6)   # aggressive: changes totals
        facts = summarize_findings()["facts"]
        assert facts["qc.library_size_max_millions"] == f"{raw_lib.max() / 1e6:.1f}"
        figs = generate_figures()["figures"]
        assert "before low-count filtering" in figs["figures/library_sizes.png"]
        assert f"{raw_lib.max() / 1e6:.1f} million" in figs["figures/library_sizes.png"]

    def test_tables_come_from_results(self, tmp_path, monkeypatch):
        from agents.analysis.tools import _ranked
        self._ready(tmp_path, monkeypatch)
        content = Path(write_report(self._full("x"))["report_path"]).read_text()
        top_gene = _ranked("up").index[0]
        assert f"| Gene{int(top_gene[1:])} | {top_gene} |" in content
        assert "| Genes tested |" in content
        assert "regulation of transcription (GO:0006355)" in content

    def test_gene_and_cite(self, tmp_path, monkeypatch):
        from agents.analysis.tools import _ranked
        self._ready(tmp_path, monkeypatch)
        SESSION.references["33010822"] = {"authors": ["Jacob F", "Pather SR"], "year": "2020"}
        gene = f"Gene{int(_ranked('up').index[0][1:])}"
        content = Path(write_report(self._full(f"{{{{gene:{gene}}}}} rose {{{{cite:33010822}}}}."))
                       ["report_path"]).read_text()
        assert f"{gene} (log2FC " in content and ", padj " in content
        assert "(Jacob et al., 2020; PMID: 33010822)" in content

    def test_typed_numbers_rejected(self, tmp_path, monkeypatch):
        self._ready(tmp_path, monkeypatch)
        result = write_report(self._full("Of 14,937 genes, 16.7% were significant (padj < 0.05)."))
        assert result["error"] == "report_problems"
        assert result["message"].count("number typed directly") == 3
        for number in ("14,937", "16.7%", "0.05"):
            assert number in result["message"]

    def test_names_and_design_labels_with_digits_allowed(self, tmp_path, monkeypatch):
        self._ready(tmp_path, monkeypatch)
        SESSION.design_df["timepoint"] = "24 hours"
        SESSION.design_df["infection"] = "SARS-CoV-2, MOI = 1.0"
        md = self._full("## 1. Summary\n\n1. Gene12 and Gene3 in GSE164073 at 24 hours (MOI = 1.0); "
                        "SARS-CoV-2, PC1, log2FC, GO:0006355 and S100A8/9 are names.")
        assert "report_path" in write_report(md)

    def test_gene_named_with_term_must_be_one_of_its_genes(self, tmp_path, monkeypatch):
        self._ready(tmp_path, monkeypatch)
        result = write_report(self._full("{{enrichment.up.go_bp.1}} includes Gene1 and Gene7."))
        assert result["error"] == "report_problems"
        assert "'Gene7' is named with 'regulation of transcription (GO:0006355)'" in result["message"]
        assert "'Gene1'" not in result["message"]
        result = write_report(self._full("{{term:up:regulation of transcription (GO:0006355)}} includes Gene7."))
        assert "'Gene7' is named with" in result["message"]

    def test_gene_placeholder_and_other_sentences_exempt_from_term_check(self, tmp_path, monkeypatch):
        self._ready(tmp_path, monkeypatch)
        md = self._full("{{enrichment.up.go_bp.1}} includes Gene1 and Gene2, and {{gene:Gene7}}. "
                        "Separately, Gene7 rose.")
        assert "report_path" in write_report(md)

    def test_all_problems_reported_at_once(self, tmp_path, monkeypatch):
        self._ready(tmp_path, monkeypatch)
        result = write_report("{{de.n_madeup}} {{gene:NOTAGENE}} {{cite:123}} ![x](figures/nope.png)")
        assert result["error"] == "report_problems"
        msg = result["message"]
        for bit in ("de.n_madeup", "gene:NOTAGENE", "cite:123", "nope.png",
                    "table:de_summary", "table:enrichment_down"):
            assert bit in msg
        assert "de.n_significant" in result["available_facts"]
        assert not SESSION.paths.analysis_report.exists()

    def test_written_with_markers_after_two_rejections(self, tmp_path, monkeypatch):
        self._ready(tmp_path, monkeypatch)
        bad = "Result: {{de.n_madeup}}. Kept text."
        assert "error" in write_report(bad)
        assert "error" in write_report(bad)
        result = write_report(bad)
        content = Path(result["report_path"]).read_text()
        assert "⚠[unknown: {{de.n_madeup}}]" in content
        assert "Kept text." in content
        assert "## ⚠ Required tables the report did not place" in content
        assert "| Genes tested |" in content
        assert "warning" in result

    def test_bare_image_reference_fixed_losslessly(self, tmp_path, monkeypatch):
        self._ready(tmp_path, monkeypatch)
        md = self._full("x").replace("![volcano](figures/volcano.png)", "![figures/volcano.png]")
        assert "![figures/volcano.png]" in md
        content = Path(write_report(md)["report_path"]).read_text()
        assert "![volcano](figures/volcano.png)\n\n*File: `figures/volcano.png` — Volcano plot:" in content

    def test_bare_reference_to_unknown_figure_rejected(self, tmp_path, monkeypatch):
        self._ready(tmp_path, monkeypatch)
        result = write_report(self._full("![figures/made_up.png]"))
        assert result["error"] == "report_problems"
        assert "malformed image reference ![figures/made_up.png]" in result["message"]

    def test_every_figure_required(self, tmp_path, monkeypatch):
        self._ready(tmp_path, monkeypatch)
        md = self._full("x").replace("![pca](figures/pca.png)", "")
        result = write_report(md)
        assert result["error"] == "report_problems"
        assert "figure not placed: ![...](figures/pca.png)" in result["message"]
        assert "figures/pca.png" in result["required_figures"]

    def test_missing_figures_appended_when_forced(self, tmp_path, monkeypatch):
        self._ready(tmp_path, monkeypatch)
        md = self._full("x").replace("![pca](figures/pca.png)", "")
        write_report(md)
        write_report(md)
        content = Path(write_report(md)["report_path"]).read_text()
        section = content.split("## ⚠ Figures the report did not place")[1]
        assert "![pca](figures/pca.png)" in section and "*File: `figures/pca.png`" in section

    def test_rejection_counter_resets_with_new_figures(self, tmp_path, monkeypatch):
        self._ready(tmp_path, monkeypatch)
        write_report("{{nope}}")
        write_report("{{nope}}")
        generate_figures()
        assert write_report("{{nope}}")["error"] == "report_problems"

    def test_no_de_only_methods_required(self, count_matrix_tsv):
        SESSION.begin_run("test")
        load_counts(str(count_matrix_tsv))
        generate_figures()
        assert "{{table:methods}}" in write_report("QC only.\n\n" + _figs())["message"]
        result = write_report("QC only: {{samples.n}} samples.\n\n{{table:qc}}\n\n{{table:methods}}\n\n" + _figs())
        assert "error" not in result
        content = Path(result["report_path"]).read_text()
        assert "4 samples" in content
        assert "### Differential expression" not in content

    def test_methods_section_is_code_only(self, tmp_path, monkeypatch):
        self._ready(tmp_path, monkeypatch)
        body = ("## Results\n\n" + TABLES.replace("{{table:methods}}", "") + "{{table:enrichment_up}}\n\n"
                "{{table:enrichment_down}}\n\n" + _figs())
        bad = write_report(body + "## 9. Methods\n\nThe Wald test was applied to log counts.\n\n"
                           "{{table:methods}}\n\n## References\n\nNone.")
        assert bad["error"] == "report_problems" and "'## 9. Methods' must contain only" in bad["message"]
        good = write_report(body + "## 9. Methods\n\n{{table:methods}}\n\n## References\n\nNone.")
        content = Path(good["report_path"]).read_text()
        assert "negative binomial" in content and "Wald test" in content
        assert "Genes were kept when at least 2 samples had ≥ 10 counts" in content
        assert "enrichment_input_upregulated.txt" in content and "Enrichr's default background" in content
        assert "| PyDESeq2 |" in content


class TestMethods:
    def test_nfcore_reference_and_rounding(self, nfcore_results_dir: Path):
        SESSION.begin_run("test")
        info = nfcore_results_dir / "pipeline_info"
        info.mkdir()
        (info / "params_2026-01-01_00-00-00.json").write_text(json.dumps({
            "genome": "GRCh38",
            "gtf": "s3://ngi-igenomes/igenomes//Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf"}))
        scan_results(str(nfcore_results_dir))
        load_counts(str(nfcore_results_dir / "star_salmon" / "salmon.merged.gene_counts.tsv"))
        SESSION.deseq_rounded_pct = 37.5
        facts = summarize_findings()["facts"]
        assert facts["reference.annotation_provider"] == "NCBI"
        assert facts["reference.annotation"] == "iGenomes Homo sapiens NCBI GRCh38"
        methods = tools._methods()
        assert "Reference genome: `GRCh38`" in methods and "iGenomes Homo sapiens NCBI GRCh38" in methods

    def test_cached_reference_reports_original_annotation(self, nfcore_results_dir: Path, tmp_path: Path):
        SESSION.begin_run("test")
        cache = tmp_path / "reference" / "GRCh38" / "nf-core-rnaseq-3.26.0"
        cache.mkdir(parents=True)
        (cache / "reference.json").write_text(json.dumps({
            "gtf_source": "s3://ngi-igenomes/igenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf"}))
        info = nfcore_results_dir / "pipeline_info"
        info.mkdir()
        (info / "params_x.json").write_text(json.dumps({"gtf": str(cache / "genes.gtf")}))
        scan_results(str(nfcore_results_dir))
        ref = tools._reference_provenance()
        assert ref["annotation_provider"] == "NCBI" and ref["cached_reference"] == str(cache)


class TestPcaGroups:
    def test_spread_and_misclustered(self, tmp_path: Path):
        _analysed(tmp_path)
        facts = summarize_findings()["facts"]
        for key in ("qc.pca_spread.wt", "qc.pca_spread.ko", "qc.pca_centroid_distance",
                    "qc.pca_condition_r2_pc1", "qc.pca_misclustered", "qc.pca_summary"):
            assert key in facts
        assert facts["qc.pca_misclustered"] == "none"
        caption = {f["path"]: f["caption"] for f in generate_figures() and SESSION.figures}["figures/pca.png"]
        assert facts["qc.pca_summary"] in caption

    def test_flags_sample_closer_to_other_group(self, tmp_path: Path):
        _analysed(tmp_path)
        wt = [s for s in SESSION.counts_df.columns if SESSION.design_df.loc[s, "condition"] == "WT"]
        SESSION.design_df.loc[wt[0], "condition"] = "KO"  # a mislabelled sample sits with the other group
        groups = tools._pca_groups()
        assert wt[0] in groups["misclustered"]

    def test_groups_compared_within_covariate_levels(self, tmp_path: Path):
        """Tissue dominates PC1 and the infection effect differs by tissue (like GSE164073):
        pooled condition centroids nearly coincide, so pooled nearest-centroid flags samples
        that are in fact well separated within their tissue."""
        import numpy as np
        import pandas as pd
        rng = np.random.default_rng(2)
        n_genes = 300
        base = rng.integers(200, 2000, size=n_genes).astype(float)
        samples, rows = [], []
        for tissue in ("cornea", "sclera"):
            for cond in ("mock", "CoV2"):
                for rep in (1, 2, 3):
                    fold = np.ones(n_genes)
                    if tissue == "sclera":
                        fold[100:200] = 20.0  # tissue effect
                    if cond == "CoV2":
                        fold[:40] = 8.0 if tissue == "cornea" else 0.125  # tissue-specific response
                    samples.append(f"{tissue}_{cond}_{rep}")
                    rows.append({"sample": samples[-1], "condition": cond, "tissue": tissue})
                    rows[-1]["counts"] = rng.poisson(base * fold)
        pd.DataFrame({r["sample"]: r["counts"] for r in rows}, index=[f"G{i:03d}" for i in range(n_genes)]) \
            .rename_axis("gene_id").to_csv(tmp_path / "counts.tsv", sep="\t")
        pd.DataFrame([{k: v for k, v in r.items() if k != "counts"} for r in rows]).to_csv(
            tmp_path / "design.csv", index=False)
        SESSION.begin_run("test")
        load_counts(str(tmp_path / "counts.tsv"), design_path=str(tmp_path / "design.csv"))
        filter_low_counts(min_count=10, min_samples=2)

        run_deseq2(["condition", "CoV2", "mock"])
        assert tools._pca_groups()["misclustered"]  # pooled comparison: the misleading result

        run_deseq2(["condition", "CoV2", "mock"], covariates=["tissue"])
        groups = tools._pca_groups()
        assert groups["misclustered"] == {}
        assert groups["covariate_r2"]["tissue"]["PC1"] > 0.9
        facts = summarize_findings()["facts"]
        assert facts["qc.pca_misclustered"] == "none"
        assert facts["qc.pca_r2_pc1.tissue"] == f"{groups['covariate_r2']['tissue']['PC1']:.0%}"
        assert "within each tissue level" in facts["qc.pca_summary"]
        assert "tissue explains" in facts["qc.pca_summary"]


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

    def test_geo_counts_loaded_from_another_run(self, tmp_path: Path):
        """A replay writes to a new run directory but loads the original run's counts."""
        original = tmp_path / "original"
        original.mkdir()
        counts, design = _strong_de_dataset(original)
        (original / "counts_metadata.json").write_text(json.dumps({
            "accession": "GSE1", "source": "author", "filename": "GSE1_counts.txt.gz", "md5": "abc",
            "value_type": "raw_integer_counts", "gene_id_type": "symbol", "n_duplicates_summed": 0,
        }))
        SESSION.begin_run("replay")
        load_counts(str(counts), design_path=str(design))
        facts = summarize_findings()["facts"]
        assert facts["provenance.data_source"] == "GEO count matrix"
        assert facts["provenance.file"] == "GSE1_counts.txt.gz"

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

    def test_repo_paths_written_relative(self, tmp_path: Path):
        from agents.analysis.replay import write_replay_script
        from core import config

        log = self._log(tmp_path / "run" / "analysis" / "tool_calls.jsonl", [
            ("load_counts", {"counts_path": str(config.ROOT / "runs" / "r" / "counts.tsv"),
                             "design_path": "/elsewhere/design.csv"}, {"n_genes": 20}),
        ])
        text = write_replay_script(log, log.parent / "replay.py").read_text()
        assert "counts_path='runs/r/counts.tsv'" in text
        assert "design_path='/elsewhere/design.csv'" in text  # outside the repo: unchanged

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


def _report_parts() -> str:
    """Everything write_report requires, for tests that only vary one section."""
    required = "".join(f"{{{{table:{t}}}}}\n\n" for t in tools._required_tables())
    return required + _figs()


class TestEvidenceTools:
    def _ready(self, tmp_path, monkeypatch):
        _analysed(tmp_path)
        _mock_enrichr(monkeypatch)
        run_enrichment("up")
        run_enrichment("down")

    def test_query_genes_by_symbol_and_family(self, tmp_path, monkeypatch):
        self._ready(tmp_path, monkeypatch)
        result = tools.query_genes(symbols=["gene5", "Gene45", "NOPE"], prefix="Gene2")
        status = {g["gene"]: g["status"] for g in result["genes"]}
        assert status == {"Gene5": "up", "Gene45": "down"}
        assert result["not_found"] == {"NOPE": "not in the data (check the symbol)"}
        family = result["family"]
        assert family["n_genes"] == 111 and family["n_up"] == 11 and family["n_down"] == 0
        assert family["cite_as"] == "{{genes:Gene2*}}"
        assert tools._render_gene_family("Gene2*") == "Gene2* genes: 11 of 111 significant (11 up, 0 down)"

    def test_search_enrichment_beyond_top_terms(self, tmp_path, monkeypatch):
        self._ready(tmp_path, monkeypatch)
        result = tools.search_enrichment("pathway", direction="up")
        (match,) = result["matches"]
        assert match["term"] == "Pathway X" and match["rank_in_library"] == 1
        assert match["genes"] == ["Gene3"] and match["cite_as"] == "{{term:up:Pathway X}}"
        assert tools._render_term("up:pathway x") == "Pathway X (3/40 genes, padj 0.020)"
        assert tools._render_term("up:Not a term") is None

    def test_paper_genes_table(self, tmp_path, monkeypatch):
        self._ready(tmp_path, monkeypatch)
        SESSION.references["123"] = {"title": "T", "authors": ["Doe J"], "year": "2024",
                                     "abstract": "KO raised Gene5 and lowered Gene45; HE cells kept Gene250. gene7 no."}
        rows = {r["gene"]: r["result"] for r in summarize_findings()["paper_genes"]}
        assert rows == {"Gene5": "up", "Gene45": "down", "Gene250": "not significant"}
        generate_figures()
        assert "paper_genes" in tools._required_tables()
        content = Path(write_report("# R\n\n" + _report_parts())["report_path"]).read_text()
        assert "| Gene5 | 123 |" in content

    def test_interpretation_must_be_short_and_anchored(self, tmp_path, monkeypatch):
        self._ready(tmp_path, monkeypatch)
        generate_figures()
        bare = write_report(_report_parts() + "## Biological Interpretation\n\n"
                            "{{gene:Gene5}} rises.\n\nThis is a hallmark of everything.\n")
        assert "points at no result" in bare["message"] and "hallmark" in bare["message"]
        long = write_report(_report_parts() + "## Biological Interpretation\n\n{{gene:Gene5}} " + "word " * 260)
        assert "keep it to 250" in long["message"]
        good = write_report(_report_parts() + "## Biological Interpretation\n\n"
                            "{{gene:Gene5}} rises.\n\n- {{genes:Gene2*}}\n- {{term:up:Pathway X}}\n")
        content = Path(good["report_path"]).read_text()
        assert "Gene2* genes: 11 of 111 significant" in content and "Pathway X (3/40 genes" in content

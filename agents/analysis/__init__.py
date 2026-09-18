"""Post-pipeline analysis agent: differential expression, enrichment, and reporting.

Takes nf-core/rnaseq outputs (count matrices, QC data) and runs downstream analysis:
QC checks, differential expression via PyDESeq2, gene set enrichment via gseapy,
and generates a Markdown report with figures.
"""

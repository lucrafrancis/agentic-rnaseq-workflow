"""JSON tool definitions sent to Claude, plus the name -> function dispatch table.

Kept in its own module so the loop stays generic.
"""

from __future__ import annotations

from typing import Any, Callable

from agents.submission import tools

TOOL_SCHEMAS: list[dict[str, Any]] = [
    {
        "name": "check_resources",
        "description": (
            "Check the machine's available resources: RAM (GB), CPU count, "
            "disk space, and Docker availability/memory. Call this before "
            "configure_submission to set appropriate resource limits."
        ),
        "input_schema": {
            "type": "object",
            "properties": {},
        },
    },
    {
        "name": "configure_submission",
        "description": (
            "Configure nextflow submission parameters for nf-core/rnaseq. "
            "Set genome, profile, and any extra nf-core/rnaseq parameters based on "
            "the user's prompt. The samplesheet path and output directory are set "
            "automatically from the session — do not provide them. "
            "Saves params.json to the run directory and returns the full nextflow command."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "genome": {
                    "type": "string",
                    "description": (
                        "Reference genome identifier for nf-core igenomes. "
                        "Infer from the organism in the prompt. "
                        "Common values: GRCh38 (human), GRCh37 (human legacy), "
                        "GRCm39 (mouse), GRCm38 (mouse legacy), BDGP6 (Drosophila), "
                        "WBcel235 (C. elegans), R64-1-1 (yeast)."
                    ),
                },
                "profile": {
                    "type": "string",
                    "description": (
                        "Nextflow execution profile. Default: docker. "
                        "Use 'singularity' if the user mentions Singularity or HPC."
                    ),
                },
                "extra_params": {
                    "type": "object",
                    "description": (
                        "Any additional nf-core/rnaseq parameters as key-value pairs. "
                        "Boolean params use true/false. The LLM is free to set any valid "
                        "nf-core/rnaseq parameter here. Common examples:\n"
                        "  skip_alignment: true — salmon pseudo-alignment only (faster, less RAM)\n"
                        "  max_memory: '30.GB' — cap memory (set below available to leave headroom)\n"
                        "  max_cpus: 4 — cap CPU usage\n"
                        "  skip_trimming: true — skip adapter trimming\n"
                        "  aligner: 'star_salmon' — aligner choice (star_salmon, star_rsem, hisat2)\n"
                        "  pseudo_aligner: 'salmon' — pseudo-aligner\n"
                        "  min_mapped_reads: 5 — minimum % mapped reads\n"
                        "Not limited to these — any valid nf-core/rnaseq param is accepted."
                    ),
                    "additionalProperties": True,
                },
            },
            "required": ["genome"],
        },
    },
]

TOOL_FUNCTIONS: dict[str, Callable[..., dict[str, Any]]] = {
    "check_resources": tools.check_resources,
    "configure_submission": tools.configure_submission,
}

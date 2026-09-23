"""The agent loop for data download resolution and script generation."""

from __future__ import annotations

from agents.download.prompts import COUNTS_SYSTEM_PROMPT, SYSTEM_PROMPT
from agents.download.schemas import COUNTS_TOOL_FUNCTIONS, COUNTS_TOOL_SCHEMAS, TOOL_FUNCTIONS, TOOL_SCHEMAS
from core.loop import run_agent_loop


def run_download_agent(user_prompt: str) -> list[dict]:
    """Resolve accessions and generate a download script if needed."""
    return run_agent_loop(
        system_prompt=SYSTEM_PROMPT,
        tool_schemas=TOOL_SCHEMAS,
        tool_functions=TOOL_FUNCTIONS,
        user_prompt=user_prompt,
        label="download",
        emoji="\U0001f310",
    )


def run_counts_agent(user_prompt: str) -> list[dict]:
    """Fetch a processed count matrix + design from GEO (analysis-only runs)."""
    return run_agent_loop(
        system_prompt=COUNTS_SYSTEM_PROMPT,
        tool_schemas=COUNTS_TOOL_SCHEMAS,
        tool_functions=COUNTS_TOOL_FUNCTIONS,
        user_prompt=user_prompt,
        label="counts",
        emoji="\U0001f310",
        hide_args=frozenset({"sample_map", "conditions"}),
    )

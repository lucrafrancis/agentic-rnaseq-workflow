"""The agent loop for data download resolution and script generation."""

from __future__ import annotations

from agents.download.prompts import SYSTEM_PROMPT
from agents.download.schemas import TOOL_FUNCTIONS, TOOL_SCHEMAS
from core.loop import run_agent_loop


def run_download_agent(user_prompt: str) -> list[dict]:
    """Resolve accessions and generate a download script if needed."""
    return run_agent_loop(
        system_prompt=SYSTEM_PROMPT,
        tool_schemas=TOOL_SCHEMAS,
        tool_functions=TOOL_FUNCTIONS,
        user_prompt=user_prompt,
        emoji="\U0001f310",
    )

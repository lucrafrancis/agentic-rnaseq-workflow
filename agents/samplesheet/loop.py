"""The agent loop for sample sheet generation."""

from __future__ import annotations

from agents.samplesheet.prompts import SYSTEM_PROMPT
from agents.samplesheet.schemas import TOOL_FUNCTIONS, TOOL_SCHEMAS
from core.loop import run_agent_loop


def run_samplesheet_agent(user_prompt: str) -> list[dict]:
    """Drive sample sheet generation to completion. Returns the final message transcript."""
    return run_agent_loop(
        system_prompt=SYSTEM_PROMPT,
        tool_schemas=TOOL_SCHEMAS,
        tool_functions=TOOL_FUNCTIONS,
        user_prompt=user_prompt,
        label="samplesheet",
    )

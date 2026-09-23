"""The agent loop for submission configuration."""

from __future__ import annotations

from agents.submission.prompts import SYSTEM_PROMPT
from agents.submission.schemas import TOOL_FUNCTIONS, TOOL_SCHEMAS
from core.loop import run_agent_loop


def run_submission_agent(user_prompt: str) -> list[dict]:
    """Configure submission parameters based on the user's prompt."""
    return run_agent_loop(
        system_prompt=SYSTEM_PROMPT,
        tool_schemas=TOOL_SCHEMAS,
        tool_functions=TOOL_FUNCTIONS,
        user_prompt=user_prompt,
        label="submission",
        emoji="⚙️",
        hide_args=frozenset({"extra_params"}),
    )

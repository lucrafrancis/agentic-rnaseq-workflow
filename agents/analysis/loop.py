"""The agent loop for downstream analysis."""

from __future__ import annotations

from agents.analysis.prompts import SYSTEM_PROMPT
from agents.analysis.schemas import TOOL_FUNCTIONS, TOOL_SCHEMAS
from core import config
from core.loop import run_agent_loop
from core.session import SESSION


def run_analysis_agent(user_prompt: str | list) -> list[dict]:
    """Drive the analysis to completion. Returns the full message transcript.

    user_prompt can be a string or a list of content blocks (for PDF inclusion).
    """
    return run_agent_loop(
        system_prompt=SYSTEM_PROMPT,
        tool_schemas=TOOL_SCHEMAS,
        tool_functions=TOOL_FUNCTIONS,
        user_prompt=user_prompt,
        label="analysis",
        log_path=SESSION.require_paths().analysis_tool_log,
        hide_args=frozenset({"report_markdown", "gene_list", "rows"}),
        model=config.ANALYSIS_MODEL,
    )

"""The agent loop for sample sheet generation.

Same pattern as the scrna workflow: send conversation + tool schemas to Claude, dispatch
tool calls, append results, repeat until the agent finishes or MAX_TURNS is hit.
"""

from __future__ import annotations

import json

import anthropic

from core import config
from agents.samplesheet.prompts import SYSTEM_PROMPT
from agents.samplesheet.schemas import TOOL_FUNCTIONS, TOOL_SCHEMAS
from core.session import SESSION


def _run_tool(name: str, args: dict) -> dict:
    try:
        return TOOL_FUNCTIONS[name](**args)
    except Exception as exc:
        return {"error": type(exc).__name__, "message": str(exc)}


def _log_tool_call(name: str, args: dict, summary: dict) -> None:
    log_path = SESSION.require_paths().tool_log
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as f:
        f.write(json.dumps({"tool": name, "args": args, "summary": summary}) + "\n")


def run_samplesheet_agent(user_prompt: str) -> dict:
    """Drive sample sheet generation to completion. Returns the final message transcript."""
    client = anthropic.Anthropic()
    messages: list[dict] = [{"role": "user", "content": user_prompt}]

    for _turn in range(config.MAX_TURNS):
        response = client.messages.create(
            model=config.MODEL,
            max_tokens=config.MAX_TOKENS,
            system=SYSTEM_PROMPT,
            tools=TOOL_SCHEMAS,
            messages=messages,
        )
        messages.append({"role": "assistant", "content": response.content})

        for block in response.content:
            if block.type == "text" and block.text.strip():
                print(f"\n\U0001f9e0 {block.text.strip()}")

        if response.stop_reason != "tool_use":
            if response.stop_reason != "end_turn":
                print(f"\n⚠️  Stopped early: stop_reason={response.stop_reason!r}")
            break

        tool_results = []
        for block in response.content:
            if block.type != "tool_use":
                continue
            summary = _run_tool(block.name, block.input)
            _log_tool_call(block.name, block.input, summary)
            arg_str = ", ".join(f"{k}={v}" for k, v in block.input.items() if k != "csv_content")
            print(f"\U0001f527 {block.name}({arg_str})")
            tool_results.append({
                "type": "tool_result",
                "tool_use_id": block.id,
                "content": json.dumps(summary),
            })
        messages.append({"role": "user", "content": tool_results})

    return messages

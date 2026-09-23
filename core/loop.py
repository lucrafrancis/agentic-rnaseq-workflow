"""Generic agent loop shared by all agents.

Each agent supplies its system prompt, tool schemas, and tool functions.
This module handles the conversation loop, tool dispatch, and logging.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable

import anthropic

from core import config
from core.session import SESSION


def run_agent_loop(
    *,
    system_prompt: str,
    tool_schemas: list[dict[str, Any]],
    tool_functions: dict[str, Callable[..., dict[str, Any]]],
    user_prompt: str | list,
    emoji: str = "\U0001f9e0",
    log_path: Path | None = None,
    hide_args: frozenset[str] = frozenset(),
    label: str = "agent",
) -> list[dict]:
    """Run an agent loop to completion.

    Parameters
    ----------
    system_prompt:
        The system prompt for this agent.
    tool_schemas:
        JSON tool definitions sent to Claude.
    tool_functions:
        Name -> callable dispatch table.
    user_prompt:
        The initial user message (string or content blocks).
    emoji:
        Prefix for assistant text output.
    log_path:
        Where to write tool call logs. Defaults to the run's tool_calls.jsonl.
    hide_args:
        Tool argument keys to omit from the console display.
    label:
        Name recorded in the run's usage.jsonl.
    """
    if log_path is None:
        log_path = SESSION.require_paths().tool_log

    client = anthropic.Anthropic()
    messages: list[dict] = [{"role": "user", "content": user_prompt}]
    usage = new_usage()

    for _turn in range(config.MAX_TURNS):
        response = client.messages.create(
            model=config.MODEL,
            max_tokens=config.MAX_TOKENS,
            system=system_prompt,
            tools=tool_schemas,
            messages=messages,
            # Automatic caching: the breakpoint follows the growing history, so each turn
            # re-reads tools + system + prior turns from cache. (Haiku 4.5 only caches
            # prefixes >= 4096 tokens, so early turns may not cache — that's expected.)
            cache_control={"type": "ephemeral"},
        )
        add_usage(usage, response.usage)
        messages.append({"role": "assistant", "content": response.content})

        for block in response.content:
            if block.type == "text" and block.text.strip():
                print(f"\n{emoji} {block.text.strip()}")

        if response.stop_reason != "tool_use":
            if response.stop_reason != "end_turn":
                print(f"\n⚠️  Stopped early: stop_reason={response.stop_reason!r}")
            break

        tool_results = []
        for block in response.content:
            if block.type != "tool_use":
                continue
            summary = _run_tool(tool_functions, block.name, block.input)
            _log_tool_call(log_path, block.name, block.input, summary)
            arg_str = ", ".join(
                f"{k}={v}" for k, v in block.input.items() if k not in hide_args
            )
            print(f"\U0001f527 {block.name}({arg_str})")
            tool_results.append({
                "type": "tool_result",
                "tool_use_id": block.id,
                "content": json.dumps(summary),
            })
        messages.append({"role": "user", "content": tool_results})
    else:
        print(f"\n⚠️  Agent used all {config.MAX_TURNS} turns without finishing.")

    log_usage(label, config.MODEL, usage)
    return messages


def new_usage() -> dict[str, int]:
    return {"requests": 0, "input_tokens": 0, "cache_creation_input_tokens": 0,
            "cache_read_input_tokens": 0, "output_tokens": 0}


def add_usage(totals: dict[str, int], usage: Any) -> None:
    totals["requests"] += 1
    for key in ("input_tokens", "cache_creation_input_tokens", "cache_read_input_tokens", "output_tokens"):
        totals[key] += getattr(usage, key, None) or 0


def log_usage(label: str, model: str, totals: dict[str, int]) -> None:
    """Print a one-line token summary and append it to the run's usage.jsonl."""
    total_in = totals["input_tokens"] + totals["cache_creation_input_tokens"] + totals["cache_read_input_tokens"]
    cached_pct = 100 * totals["cache_read_input_tokens"] / total_in if total_in else 0
    print(f"\n📊 {label}: {totals['requests']} requests, {total_in:,} input tokens "
          f"({cached_pct:.0f}% from cache), {totals['output_tokens']:,} output tokens")
    if SESSION.paths is None:
        return
    with SESSION.paths.usage_log.open("a") as f:
        f.write(json.dumps({"agent": label, "model": model, **totals}) + "\n")


def _run_tool(
    functions: dict[str, Callable[..., dict[str, Any]]],
    name: str,
    args: dict,
) -> dict:
    try:
        return functions[name](**args)
    except Exception as exc:
        return {"error": type(exc).__name__, "message": str(exc)}


def _log_tool_call(log_path: Path, name: str, args: dict, summary: dict) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a") as f:
        f.write(json.dumps({"tool": name, "args": args, "summary": summary}) + "\n")

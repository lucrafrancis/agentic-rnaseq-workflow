"""Tests for the shared agent loop — a fake client stands in for the API."""

from __future__ import annotations

import json
from types import SimpleNamespace

from core import loop
from core.session import SESSION


def _response(content, stop_reason, cached=0):
    usage = SimpleNamespace(input_tokens=100, cache_creation_input_tokens=50,
                            cache_read_input_tokens=cached, output_tokens=10)
    return SimpleNamespace(content=content, stop_reason=stop_reason, usage=usage)


class FakeClient:
    def __init__(self, responses):
        self.responses = list(responses)
        self.calls: list[dict] = []
        self.messages = self

    def create(self, **kwargs):
        self.calls.append(kwargs)
        return self.responses.pop(0)


def test_caches_and_logs_usage(monkeypatch):
    SESSION.begin_run("loop")
    tool_use = SimpleNamespace(type="tool_use", id="t1", name="echo", input={"x": 1})
    done = SimpleNamespace(type="text", text="done")
    client = FakeClient([_response([tool_use], "tool_use"), _response([done], "end_turn", cached=400)])
    monkeypatch.setattr(loop.anthropic, "Anthropic", lambda: client)

    loop.run_agent_loop(
        system_prompt="sys", tool_schemas=[], tool_functions={"echo": lambda x: {"x": x}},
        user_prompt="go", label="test",
    )

    assert all(c["cache_control"] == {"type": "ephemeral"} for c in client.calls)
    record = json.loads(SESSION.paths.usage_log.read_text())
    price = loop.config.PRICE_PER_MTOK[loop.config.MODEL]
    cost = (200 * price["input"] + 100 * price["input"] * 1.25 + 400 * price["input"] * 0.1
            + 20 * price["output"]) / 1e6
    assert record == {"agent": "test", "model": loop.config.MODEL, "requests": 2, "input_tokens": 200,
                      "cache_creation_input_tokens": 100, "cache_read_input_tokens": 400, "output_tokens": 20,
                      "estimated_cost_usd": round(cost, 4)}


def test_unpriced_model_logged_without_cost():
    assert loop.estimate_cost("claude-unknown", loop.new_usage()) is None


def test_model_override_used_and_logged(monkeypatch):
    SESSION.begin_run("loop")
    client = FakeClient([_response([SimpleNamespace(type="text", text="done")], "end_turn")])
    monkeypatch.setattr(loop.anthropic, "Anthropic", lambda: client)

    loop.run_agent_loop(system_prompt="sys", tool_schemas=[], tool_functions={}, user_prompt="go",
                        label="test", model="claude-other")

    assert client.calls[0]["model"] == "claude-other"
    assert json.loads(SESSION.paths.usage_log.read_text())["model"] == "claude-other"


def test_tool_log_paths_are_repo_relative(tmp_path):
    root = loop.config.ROOT
    log = tmp_path / "tool_calls.jsonl"
    loop._log_tool_call(log, "load_counts", {"counts_path": f"{root}/runs/r/counts.tsv", "other": "/tmp/x.tsv"},
                        {"files": [f"{root}/runs/r/a.png"], "message": f"wrote {root}/runs/r/b.csv", "n": 3})
    entry = json.loads(log.read_text())
    assert entry["args"] == {"counts_path": "runs/r/counts.tsv", "other": "/tmp/x.tsv"}
    assert entry["summary"] == {"files": ["runs/r/a.png"], "message": "wrote runs/r/b.csv", "n": 3}

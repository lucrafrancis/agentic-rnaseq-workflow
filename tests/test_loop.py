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
    assert record == {"agent": "test", "model": loop.config.MODEL, "requests": 2, "input_tokens": 200,
                      "cache_creation_input_tokens": 100, "cache_read_input_tokens": 400, "output_tokens": 20}

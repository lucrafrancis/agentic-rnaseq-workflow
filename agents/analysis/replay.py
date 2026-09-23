"""Turn the analysis tool log into a standalone replay script.

The agent's decisions (contrast, filter thresholds, gene lists, report text) are captured
as tool arguments in analysis/tool_calls.jsonl. Replaying the successful calls in order,
against the same code, reproduces the analysis without the LLM.

  uv run python -m agents.analysis.replay <run_dir>   — generate for an existing run
"""

from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from pprint import pformat

from core import config

# Tools that only inform the agent — no effect on outputs. Fetches also need network.
SKIP_TOOLS = frozenset({
    "fetch_geo_metadata", "fetch_abstract", "inspect_counts", "get_top_genes", "summarize_findings",
})

# Tools whose results depend on an external service and may differ on replay.
NETWORK_TOOLS = frozenset({"run_enrichment"})


def _git_commit() -> str:
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=config.ROOT, capture_output=True, text=True, check=True,
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "status", "--porcelain", "--", "agents", "core"],
            cwd=config.ROOT, capture_output=True, text=True, check=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"
    return f"{commit}-dirty" if dirty else commit


def _format_call(name: str, args: dict) -> str:
    if not args:
        return f'_run("{name}")'
    lines = [f'_run(\n    "{name}",']
    for key, value in args.items():
        formatted = pformat(value, width=88, compact=True, sort_dicts=False).replace("\n", "\n    ")
        lines.append(f"    {key}={formatted},")
    lines.append(")")
    return "\n".join(lines)


def write_replay_script(log_path: Path, out_path: Path, *, start_line: int = 0) -> Path | None:
    """Write a replay script from the log's successful, output-affecting calls.

    start_line skips earlier sessions appended to the same log (e.g. a previous --analyze).
    Returns the script path, or None if there was nothing to replay.
    """
    if not log_path.is_file():
        return None
    entries = [json.loads(line) for line in log_path.read_text().splitlines()[start_line:] if line.strip()]
    calls = [
        e for e in entries
        if e["tool"] not in SKIP_TOOLS and "error" not in e.get("summary", {})
    ]
    if not calls:
        return None

    run_name = out_path.parent.parent.name
    header = f'''"""Replay of the analysis agent's tool calls for run {run_name}.

Generated: {datetime.now().isoformat(timespec="seconds")}
Git commit: {_git_commit()}
Source log: {log_path.name} ({len(calls)} of {len(entries)} calls; failed and read-only calls omitted)

Re-runs the same tool functions with the same arguments, without the LLM. Output goes
to a new run directory (runs/<date>_{run_name}_replay/). Package versions are pinned by
uv.lock at the commit above — check it out first if the code has changed since.

Caveat: run_enrichment queries Enrichr live; its gene-set libraries change over time,
so enrichment results may differ from the original run.

  uv run python {out_path.relative_to(config.ROOT) if out_path.is_relative_to(config.ROOT) else out_path}
"""

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]  # runs/<run>/analysis/replay.py -> repo root
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)  # logged paths may be relative to the repo root

from agents.analysis import tools
from core.session import SESSION

NETWORK_TOOLS = {sorted(NETWORK_TOOLS)!r}


def _run(name, **kwargs):
    print(f"-> {{name}}")
    result = getattr(tools, name)(**kwargs)
    if "error" in result:
        if name in NETWORK_TOOLS:
            print(f"   WARNING: {{result['error']}}: {{result.get('message', '')}}")
            return result
        sys.exit(f"   FAILED: {{result['error']}}: {{result.get('message', '')}}")
    return result


SESSION.begin_run({run_name + "_replay"!r})
print(f"Replay output: {{SESSION.paths.dir}}")

'''
    body = "\n\n".join(_format_call(e["tool"], e["args"]) for e in calls)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(header + body + '\n\nprint(f"Done. Report: {SESSION.paths.analysis_report}")\n')
    return out_path


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: python -m agents.analysis.replay <run_dir>")
    paths = config.RunPaths(Path(sys.argv[1]).resolve().name)
    script = write_replay_script(paths.analysis_tool_log, paths.analysis_dir / "replay.py")
    print(f"Wrote {script}" if script else f"Nothing to replay in {paths.analysis_tool_log}")

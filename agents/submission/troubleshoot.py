"""LLM-assisted troubleshooting and post-run warning review.

Two modes:
  - Failure: interactive conversation with the LLM — it reads the full log,
    explains the issue, and proposes a fix. The user can chat, ask questions,
    and redirect before approving.
  - Success: scans the log for warnings and produces a concise summary.

All events are logged to troubleshooting.jsonl for auditability.
The LLM is only active during diagnosis — no tokens are burned while the pipeline runs.
"""

from __future__ import annotations

import json
import logging
import re
import shlex
import subprocess
from datetime import datetime
from pathlib import Path

import anthropic

from agents.submission.params import SubmissionParams
from core import config
from core.loop import add_usage, log_usage, new_usage
from core.session import SESSION

logger = logging.getLogger(__name__)

MAX_RETRIES = 3

_ALLOWED_BINARIES = frozenset({
    "cat", "df", "docker", "du", "env", "file", "free", "head",
    "hostname", "id", "java", "ls", "lsb_release", "lscpu", "lsof",
    "mount", "nextflow", "nproc", "ps", "sysctl", "tail", "top",
    "uname", "wc", "which", "whoami",
})

_SHELL_METACHARACTERS = re.compile(r"[;|&`$(){}<>!\\]")


def _run_safe_command(command: str) -> dict:
    """Run a command after validating it against the binary allowlist.

    Rejects shell metacharacters and commands whose binary isn't allowed.
    """
    command = command.strip()
    if _SHELL_METACHARACTERS.search(command):
        msg = f"Rejected: shell metacharacters not allowed in command: {command!r}"
        logger.warning(msg)
        return {"error": "blocked", "message": msg}

    try:
        parts = shlex.split(command)
    except ValueError as exc:
        return {"error": "blocked", "message": f"Could not parse command: {exc}"}

    if not parts:
        return {"error": "blocked", "message": "Empty command."}

    binary = Path(parts[0]).name
    if binary not in _ALLOWED_BINARIES:
        msg = f"Rejected: '{binary}' is not in the allowed binary list."
        logger.warning(msg)
        return {"error": "blocked", "message": msg}

    result = subprocess.run(
        parts, capture_output=True, text=True, timeout=30,
    )
    return {
        "stdout": result.stdout.strip(),
        "stderr": result.stderr.strip(),
        "returncode": result.returncode,
    }

_SYSTEM_PROMPT = """\
You are an expert nf-core/rnaseq troubleshooter. A pipeline run has failed and you
need to help the user diagnose and fix it.

Read the full nextflow log, identify what went wrong, and explain it clearly to the
user. Then call the propose_fix tool with a concrete fix.

The user may want to discuss — ask questions, provide context about their environment,
or suggest alternatives. Respond helpfully and adapt your proposal as needed. When the
user is satisfied, they'll approve your fix.

## Diagnosing the environment

You have a check_environment tool to inspect the machine (read-only). Use it to check
things the log alone can't tell you — Docker memory limits, disk space, running
processes, etc. For OOM kills (exit 137), always check `docker info` to see if
Docker's memory limit is the bottleneck, not just system RAM.

IMPORTANT: Do NOT use check_environment to run pipelines, modify files, or execute
nextflow commands. It is strictly for inspection. Retries with -resume are handled
automatically by the framework when you propose a fix and the user approves it.

## nf-core/rnaseq knowledge

Common fixes for known failure modes:
- **Memory exceeded (RSEM/STAR)**: RSEM and STAR genome indexing need lots of RAM.
  If the machine doesn't have enough, use salmon-only: `--skip_alignment` uses salmon
  pseudo-alignment without STAR/RSEM, drastically reducing memory. Alternatively,
  add `-resume` to reuse completed steps and set max memory:
  `--max_memory '30.GB'` (set below available RAM to leave headroom).
- **Docker OOM (exit 137)**: the process was killed by Docker, not the OS. Check
  Docker's memory allocation with `docker info | grep Memory`. The fix is a user
  action: increase Docker Desktop's memory in Settings → Resources.
- **Docker/Singularity not running**: user needs to start Docker Desktop or the
  Docker daemon.
- **Container pull failure**: network issue or rate-limited. Retry, or use
  `--singularity_pull_docker_container` to switch to Singularity.
- **Genome not found**: check `--genome` value against nf-core igenomes. Common
  values: GRCh38, GRCh37, GRCm39, GRCm38.
- **Missing input files**: FASTQ paths must be absolute. Check the samplesheet.
- **Disk space**: nextflow work directories can be large. Suggest `nextflow clean`.
- **Strandedness issues**: if salmon quant fails on strandedness, suggest `auto`.

## Fix categories

When calling propose_fix, classify the fix as:
- "parameter_change" — a nextflow parameter to add or change (always include -resume)
- "user_action" — something the user must do outside the pipeline (start Docker, etc.)
- "config_change" — a nextflow.config edit
- "unfixable" — beyond automated troubleshooting

If this is attempt 2+, previous fixes didn't work — try something different.
"""

_CHECK_ENV_TOOL = {
    "name": "check_environment",
    "description": (
        "Run a read-only shell command to inspect the machine and environment. "
        "Use this to check Docker memory limits, disk space, running processes, "
        "or anything else the log alone can't tell you. "
        "Do NOT use this to run pipelines, modify files, or execute nextflow."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "command": {
                "type": "string",
                "description": "The shell command to run.",
            },
        },
        "required": ["command"],
    },
}

_PROPOSE_FIX_TOOL = {
    "name": "propose_fix",
    "description": (
        "Propose a concrete fix for the pipeline failure. Call this after explaining "
        "the issue to the user. They will see your proposal and can approve, discuss "
        "further, or ask you to try a different approach."
    ),
    "input_schema": {
        "type": "object",
        "properties": {
            "root_cause": {
                "type": "string",
                "description": "One-line summary of what went wrong.",
            },
            "category": {
                "type": "string",
                "enum": ["parameter_change", "user_action", "config_change", "unfixable"],
                "description": "Type of fix needed.",
            },
            "fix_description": {
                "type": "string",
                "description": "What to change and why.",
            },
            "parameter_changes": {
                "type": "object",
                "description": "Nextflow parameters to add or change (for parameter_change).",
            },
            "user_instructions": {
                "type": "string",
                "description": "What the user needs to do (for user_action or config_change).",
            },
        },
        "required": ["root_cause", "category", "fix_description"],
    },
}


def _log_event(event: dict) -> None:
    paths = SESSION.require_paths()
    with paths.troubleshooting_log.open("a") as f:
        event["timestamp"] = datetime.now().isoformat()
        f.write(json.dumps(event) + "\n")


_MAX_LOG_BYTES = 512 * 1024

def _read_log(log_path: str) -> str:
    path = Path(log_path)
    if not path.is_file():
        return "(log file not found)"
    content = path.read_text()
    if len(content) > _MAX_LOG_BYTES:
        return content[:_MAX_LOG_BYTES // 2] + "\n\n... [truncated] ...\n\n" + content[-_MAX_LOG_BYTES // 2:]
    return content


def _print_proposal(proposal: dict) -> None:
    """Display a fix proposal to the user."""
    print(f"\n  Root cause: {proposal.get('root_cause', 'unknown')}")
    print(f"  Category:   {proposal.get('category', 'unknown')}")
    print(f"  Fix:        {proposal.get('fix_description', 'none proposed')}")

    if proposal.get("parameter_changes"):
        print(f"\n  Parameter changes:")
        for k, v in proposal["parameter_changes"].items():
            print(f"    --{k} {v}")

    if proposal.get("user_instructions"):
        print(f"\n  Action needed:")
        print(f"    {proposal['user_instructions']}")


def diagnose_and_propose(
    outcome: dict,
    params: SubmissionParams,
    attempt: int,
    history: list[dict],
) -> dict | None:
    """Interactive troubleshooting conversation with the LLM.

    The LLM reads the full nextflow log, explains the issue, and proposes a fix
    via the propose_fix tool. The user can chat — ask questions, provide context,
    or redirect — before approving. Returns the approved proposal dict, or None
    if the user quits.
    """
    client = anthropic.Anthropic()
    log_content = _read_log(outcome.get("log_path", ""))

    history_text = ""
    if history:
        history_text = "\n\nPrevious attempts:\n" + "\n".join(
            f"  Attempt {h['attempt']}: {h['root_cause']} → fix: {h['fix_description']} → outcome: {h['outcome']}"
            for h in history
        )

    initial_message = (
        f"Attempt {attempt} of {MAX_RETRIES}.\n\n"
        f"Submission parameters:\n{json.dumps(params.to_dict(), indent=2)}\n\n"
        f"Nextflow log:\n{log_content}"
        + history_text
    )

    messages: list[dict] = [{"role": "user", "content": initial_message}]
    proposal: dict | None = None

    print(f"\n{'='*60}")
    print(f"  Troubleshooting — attempt {attempt}/{MAX_RETRIES}")
    print(f"{'='*60}")

    usage = new_usage()
    try:
        _MAX_DIAGNOSTIC_TURNS = 20
        for _diag_turn in range(_MAX_DIAGNOSTIC_TURNS):
            response = client.messages.create(
                model=config.MODEL_SONNET,
                max_tokens=2048,
                system=_SYSTEM_PROMPT,
                tools=[_CHECK_ENV_TOOL, _PROPOSE_FIX_TOOL],
                messages=messages,
                cache_control={"type": "ephemeral"},  # full log is re-sent every turn
            )
            add_usage(usage, response.usage)
            messages.append({"role": "assistant", "content": response.content})

            for block in response.content:
                if hasattr(block, "text") and block.text.strip():
                    print(f"\n🔍 {block.text.strip()}")

            if response.stop_reason == "tool_use":
                tool_results = []
                has_proposal = False
                for block in response.content:
                    if block.type != "tool_use":
                        continue
                    if block.name == "check_environment":
                        result = _run_safe_command(block.input["command"])
                        print(f"\n🔧 check_environment({block.input['command']})")
                        tool_results.append({
                            "type": "tool_result",
                            "tool_use_id": block.id,
                            "content": json.dumps(result),
                        })
                    elif block.name == "propose_fix":
                        has_proposal = True
                        proposal = block.input
                        _log_event({"event": "diagnosis", "attempt": attempt, "proposal": proposal})
                        _print_proposal(proposal)
                        tool_results.append({
                            "type": "tool_result",
                            "tool_use_id": block.id,
                            "content": "Proposal shown to user. Waiting for their response.",
                        })
                messages.append({"role": "user", "content": tool_results})
                if not has_proposal:
                    continue

            try:
                if proposal:
                    raw = input("\n[a]pply fix / [q]uit / or type to discuss: ").strip()
                else:
                    raw = input("\n[q]uit / or type to discuss: ").strip()
            except EOFError:
                return None

            if raw.lower() in ("a", "apply") and proposal:
                return proposal
            if raw.lower() in ("q", "quit"):
                return None
            if not raw:
                continue

            messages.append({"role": "user", "content": raw})
            proposal = None
        else:
            print(f"\n⚠️  Diagnostic conversation hit {_MAX_DIAGNOSTIC_TURNS} turns without resolution.")
            return None
    finally:
        log_usage("troubleshoot", config.MODEL_SONNET, usage)


def apply_parameter_fix(params: SubmissionParams, changes: dict) -> None:
    """Apply parameter changes to the submission params."""
    for key, value in changes.items():
        if hasattr(params, key):
            setattr(params, key, value)
        else:
            params.extra_args[key] = value

    paths = SESSION.require_paths()
    paths.params_file.write_text(json.dumps(params.to_dict(), indent=2) + "\n")
    _log_event({"event": "params_updated", "changes": changes})


def wait_for_user_action() -> bool:
    """Wait for the user to complete a manual action before retrying."""
    print("\nComplete the action described above, then press Enter to retry (or 'q' to quit).")
    try:
        raw = input("> ").strip().lower()
    except EOFError:
        return False
    return raw != "q"


_WARNING_REVIEW_PROMPT = """\
You are reviewing the nextflow log from a successful nf-core/rnaseq run. Extract any
WARN lines or notable issues the user should be aware of. Ignore routine informational
messages.

For each warning, provide:
- What it means in plain language
- Whether it affects downstream analysis
- What the user could do about it (if anything)

Be concise — one or two sentences per warning. If there are no warnings worth
flagging, respond with exactly: "No warnings to report."
"""


def review_warnings(log_path: str) -> str | None:
    """Scan a successful run's log for warnings and return a concise LLM summary.

    Returns the summary string, or None if there are no warnings.
    """
    log_content = _read_log(log_path)
    if "(log file not found)" in log_content:
        return None

    warn_lines = [line for line in log_content.splitlines() if "WARN" in line]
    if not warn_lines:
        return None

    client = anthropic.Anthropic()
    response = client.messages.create(
        model=config.MODEL_SONNET,
        max_tokens=1024,
        system=_WARNING_REVIEW_PROMPT,
        messages=[{"role": "user", "content": f"Nextflow log warnings:\n\n" + "\n".join(warn_lines)}],
    )

    summary = response.content[0].text.strip()
    _log_event({"event": "warning_review", "n_warnings": len(warn_lines), "summary": summary})

    if summary == "No warnings to report.":
        return None
    return summary

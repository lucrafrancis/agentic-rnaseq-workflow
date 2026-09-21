"""LLM-assisted troubleshooting and post-run warning review.

Two modes:
  - Failure: feeds the error log to the LLM, which proposes a fix. User approves
    before anything changes.
  - Success: scans the log for warnings and produces a concise summary.

All events are logged to troubleshooting.jsonl for auditability.
The LLM is only active during diagnosis — no tokens are burned while the pipeline runs.
"""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import anthropic

from agents.submission.errors import Diagnosis
from agents.submission.params import SubmissionParams
from core import config
from core.session import SESSION

MAX_RETRIES = 3

_SYSTEM_PROMPT = """\
You are an expert nf-core/rnaseq troubleshooter. A pipeline run has failed and you
need to diagnose the issue and propose a concrete fix.

## nf-core/rnaseq knowledge

Common fixes for known failure modes:
- **Memory exceeded (RSEM/STAR)**: RSEM and STAR genome indexing need lots of RAM.
  If the machine doesn't have enough, use salmon-only: `--skip_alignment` uses salmon
  pseudo-alignment without STAR/RSEM, drastically reducing memory. Alternatively,
  add `-resume` to reuse completed steps and set max memory:
  `--max_memory '30.GB'` (set below available RAM to leave headroom).
- **Docker/Singularity not running**: user needs to start Docker Desktop or the
  Docker daemon.
- **Container pull failure**: network issue or rate-limited. Retry, or use
  `--singularity_pull_docker_container` to switch to Singularity.
- **Genome not found**: check `--genome` value against nf-core igenomes. Common
  values: GRCh38, GRCh37, GRCm39, GRCm38.
- **Missing input files**: FASTQ paths must be absolute. Check the samplesheet.
- **Disk space**: nextflow work directories can be large. Suggest `nextflow clean`.
- **Strandedness issues**: if salmon quant fails on strandedness, suggest `auto`.

## Instructions

1. Read the error log and identify the root cause.
2. Classify the fix:
   - "parameter_change" — a nextflow/nf-core parameter to add or change
   - "user_action" — something the user must do outside the pipeline
   - "config_change" — a nextflow.config change (memory limits, retry strategy)
   - "unfixable" — beyond automated troubleshooting
3. For parameter_change: include the `-resume` flag so completed steps aren't re-run.
4. Be specific — name the exact parameter and value.
5. If this is attempt 2+, the previous fix didn't work. Try something different.

Respond as JSON (no markdown fences):
{
  "root_cause": "one-line summary",
  "category": "parameter_change|user_action|config_change|unfixable",
  "fix_description": "what to change and why",
  "parameter_changes": {"key": "value"},
  "user_instructions": "what the user needs to do"
}

Do not wrap the JSON in markdown code fences. Return only the JSON object.
"""


def _log_event(event: dict) -> None:
    paths = SESSION.require_paths()
    with paths.troubleshooting_log.open("a") as f:
        event["timestamp"] = datetime.now().isoformat()
        f.write(json.dumps(event) + "\n")


def _get_log_tail(log_path: str, n_lines: int = 200) -> str:
    path = Path(log_path)
    if not path.is_file():
        return "(log file not found)"
    lines = path.read_text().splitlines()
    return "\n".join(lines[-n_lines:])


def _parse_json_response(raw: str) -> dict:
    """Extract JSON from an LLM response, handling markdown fences."""
    try:
        return json.loads(raw)
    except json.JSONDecodeError:
        pass
    # Strip markdown code fences
    import re
    match = re.search(r"```(?:json)?\s*\n(.*?)\n```", raw, re.DOTALL)
    if match:
        try:
            return json.loads(match.group(1))
        except json.JSONDecodeError:
            pass
    # Try finding first { to last }
    start = raw.find("{")
    end = raw.rfind("}")
    if start != -1 and end != -1:
        try:
            return json.loads(raw[start : end + 1])
        except json.JSONDecodeError:
            pass
    return {"category": "unfixable", "root_cause": "Could not parse LLM response", "raw_response": raw}


def diagnose_and_propose(
    outcome: dict,
    params: SubmissionParams,
    diagnoses: list[Diagnosis],
    attempt: int,
    history: list[dict],
) -> dict:
    """Ask the LLM to diagnose a failure and propose a fix.

    Returns the parsed proposal dict, or an error dict if parsing fails.
    """
    log_tail = _get_log_tail(outcome.get("log_path", ""))

    history_text = ""
    if history:
        history_text = "\n\nPrevious attempts:\n" + "\n".join(
            f"  Attempt {h['attempt']}: {h['root_cause']} → fix: {h['fix_description']} → outcome: {h['outcome']}"
            for h in history
        )

    user_message = (
        f"Attempt {attempt} of {MAX_RETRIES}.\n\n"
        f"Submission parameters:\n{json.dumps(params.to_dict(), indent=2)}\n\n"
        f"Pattern-matched diagnoses:\n"
        + ("\n".join(f"  - [{d.category}] {d.line}\n    {d.suggestion}" for d in diagnoses) or "  (none matched)")
        + f"\n\nNextflow log (last 200 lines):\n{log_tail}"
        + history_text
    )

    client = anthropic.Anthropic()
    response = client.messages.create(
        model=config.MODEL_STRONG,
        max_tokens=2048,
        system=_SYSTEM_PROMPT,
        messages=[{"role": "user", "content": user_message}],
    )

    raw = response.content[0].text.strip()
    proposal = _parse_json_response(raw)

    _log_event({"event": "diagnosis", "attempt": attempt, "proposal": proposal})
    return proposal


def present_proposal(proposal: dict, attempt: int) -> str:
    """Show the LLM's proposal to the user and get their decision.

    Returns: "approve", "skip", or "abort".
    """
    print(f"\n{'='*60}")
    print(f"  Troubleshooting — attempt {attempt}/{MAX_RETRIES}")
    print(f"{'='*60}\n")
    print(f"Root cause: {proposal.get('root_cause', 'unknown')}")
    print(f"Category:   {proposal.get('category', 'unknown')}")
    print(f"Fix:        {proposal.get('fix_description', 'none proposed')}")

    if proposal.get("parameter_changes"):
        print(f"\nParameter changes:")
        for k, v in proposal["parameter_changes"].items():
            print(f"  --{k} {v}")

    if proposal.get("user_instructions"):
        print(f"\nAction required from you:")
        print(f"  {proposal['user_instructions']}")

    print()
    while True:
        try:
            raw = input("[a]pply fix & retry / [s]kip troubleshooting / [q]uit pipeline: ")
        except EOFError:
            return "abort"
        choice = raw.strip().lower()
        if choice in ("a", "apply"):
            return "approve"
        if choice in ("s", "skip"):
            return "skip"
        if choice in ("q", "quit"):
            return "abort"
        print("Please enter 'a', 's', or 'q'.")


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


def wait_for_user_action(proposal: dict) -> bool:
    """Wait for the user to complete a manual action. Returns True if they're ready to retry."""
    print(f"\nComplete the action above, then press Enter to retry (or 'q' to quit).")
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
    log_tail = _get_log_tail(log_path, n_lines=500)
    if "(log file not found)" in log_tail:
        return None

    warn_lines = [line for line in log_tail.splitlines() if "WARN" in line]
    if not warn_lines:
        return None

    client = anthropic.Anthropic()
    response = client.messages.create(
        model=config.MODEL_STRONG,
        max_tokens=1024,
        system=_WARNING_REVIEW_PROMPT,
        messages=[{"role": "user", "content": f"Nextflow log warnings:\n\n" + "\n".join(warn_lines)}],
    )

    summary = response.content[0].text.strip()
    _log_event({"event": "warning_review", "n_warnings": len(warn_lines), "summary": summary})

    if summary == "No warnings to report.":
        return None
    return summary

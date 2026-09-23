"""Data download tools.

THE CONTRACT:
  Every tool returns a structured, JSON-serializable summary dict. The LLM decides
  the next step entirely from it. URLs and MD5 checksums are saved to disk and never
  relayed through the LLM — downstream tools read them back directly.
"""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

Summary = dict[str, Any]

_MAX_RETRIES = 3
_BACKOFF_BASE = 2.0

_ENA_FIELDS = (
    "run_accession,fastq_ftp,fastq_aspera,fastq_md5,fastq_bytes,library_layout,sample_alias,sample_title"
)
_ENA_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"
_NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
_NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


class LookupFailed(RuntimeError):
    """NCBI/ENA could not be queried (network error, rate limiting, malformed reply).

    Distinct from "queried fine, nothing found" — resolvers return [] only for that.
    """


def _http_get(url: str, timeout: int = 30) -> str:
    req = urllib.request.Request(
        url, headers={"User-Agent": "agentic-rnaseq-workflow/0.1"}
    )
    for attempt in range(_MAX_RETRIES):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read().decode("utf-8")
        except urllib.error.HTTPError as e:
            if e.code == 429 and attempt < _MAX_RETRIES - 1:
                time.sleep(_BACKOFF_BASE ** (attempt + 1))
                continue
            raise
    return ""


def _http_get_json(url: str, timeout: int = 30) -> Any:
    """GET and parse JSON. Retries 429s and empty/non-JSON bodies (NCBI returns these
    when throttling); raises LookupFailed once retries are exhausted."""
    req = urllib.request.Request(
        url,
        headers={
            "User-Agent": "agentic-rnaseq-workflow/0.1",
            "Accept": "application/json",
        },
    )
    last_error = "no response"
    for attempt in range(_MAX_RETRIES):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.loads(resp.read().decode("utf-8"))
        except urllib.error.HTTPError as e:
            last_error = f"HTTP {e.code}"
            if e.code != 429:
                break
        except json.JSONDecodeError:
            last_error = "empty or non-JSON reply"
        except urllib.error.URLError as e:
            last_error = str(e.reason)
            break
        if attempt < _MAX_RETRIES - 1:
            time.sleep(_BACKOFF_BASE ** (attempt + 1))
    raise LookupFailed(f"{url.split('?')[0]}: {last_error}")


def _compute_md5(filepath: Path) -> str:
    h = hashlib.md5()
    with filepath.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _parse_ena_run(row: dict) -> dict:
    ftp_raw = row.get("fastq_ftp", "")
    aspera_raw = row.get("fastq_aspera", "")
    md5_raw = row.get("fastq_md5", "")
    bytes_raw = row.get("fastq_bytes", "")

    ftp_urls = [u.strip() for u in ftp_raw.split(";") if u.strip()]
    aspera_urls = [u.strip() for u in aspera_raw.split(";") if u.strip()]
    md5s = [m.strip() for m in md5_raw.split(";") if m.strip()]
    sizes = [s.strip() for s in bytes_raw.split(";") if s.strip()]

    files = []
    total_bytes = 0
    for i, ftp in enumerate(ftp_urls):
        filename = ftp.split("/")[-1]
        size = int(sizes[i]) if i < len(sizes) and sizes[i] else 0
        md5 = md5s[i] if i < len(md5s) else ""
        total_bytes += size
        if ftp.startswith("ftp://"):
            url = ftp.replace("ftp://", "http://", 1)
        elif ftp.startswith("http://") or ftp.startswith("https://"):
            url = ftp
        else:
            url = f"http://{ftp}"
        entry = {"filename": filename, "url": url, "md5": md5, "bytes": size}
        if i < len(aspera_urls) and aspera_urls[i]:
            entry["aspera"] = aspera_urls[i]
        files.append(entry)

    return {
        "run_accession": row.get("run_accession", ""),
        "library_layout": row.get("library_layout", ""),
        "sample_alias": row.get("sample_alias", ""),
        "sample_title": row.get("sample_title", ""),
        "files": files,
        "total_bytes": total_bytes,
    }


# --- Internal resolution helpers ---


def _resolve_gse(gse: str) -> list[dict]:
    """GSE -> GDS record -> BioProject/SRP -> ENA FASTQ URLs.

    Returns [] when the series has no linked SRA project; raises LookupFailed when
    NCBI/ENA can't be queried.
    """
    data = _http_get_json(f"{_NCBI_ESEARCH}?db=gds&term={gse}[ACCN]&retmode=json")
    result = data.get("esearchresult")
    if result is None or "ERROR" in result:
        raise LookupFailed(f"NCBI esearch returned no result for {gse}: {data}")
    uid_list = result.get("idlist", [])
    if not uid_list:
        return []

    try:
        record = _http_get(f"{_NCBI_EFETCH}?db=gds&id={uid_list[0]}&retmode=text")
    except urllib.error.URLError as exc:
        raise LookupFailed(f"NCBI efetch failed for {gse}: {exc}") from exc
    if not record.strip():
        raise LookupFailed(f"NCBI efetch returned an empty record for {gse}")

    match = re.search(r"acc=(PRJNA\d+|SRP\d+|ERP\d+|DRP\d+)", record)
    if match:
        return _resolve_study(match.group(1))

    return []


def _resolve_study(study: str) -> list[dict]:
    """Query ENA for all runs in an SRA study / BioProject."""
    url = (
        f"{_ENA_API}?accession={study}&result=read_run"
        f"&fields={_ENA_FIELDS}&format=json"
    )
    data = _http_get_json(url)
    if not isinstance(data, list):
        raise LookupFailed(f"Unexpected ENA reply for {study}: {str(data)[:200]}")
    return [_parse_ena_run(r) for r in data if r.get("run_accession")]


def _resolve_runs(srr_ids: list[str]) -> list[dict]:
    """Query ENA for individual SRR accessions."""
    runs = []
    for srr in srr_ids:
        url = (
            f"{_ENA_API}?accession={srr}&result=read_run"
            f"&fields={_ENA_FIELDS}&format=json"
        )
        data = _http_get_json(url)
        if isinstance(data, list):
            runs.extend(_parse_ena_run(r) for r in data if r.get("run_accession"))
    return runs


# --- Agent tools (called by the LLM via the agent loop) ---


def resolve_accession(accession: str) -> Summary:
    """Resolve a GEO/SRA accession to download metadata.

    Queries NCBI and ENA APIs. Saves full metadata (URLs, MD5s) to
    download_metadata.json — downstream tools read it directly so
    URLs and checksums never pass through the LLM.
    """
    from core.session import SESSION

    accession = accession.strip()
    upper = accession.upper()

    if not upper.startswith(("GSE", "SRP", "ERP", "DRP", "PRJNA", "SRR")):
        return {
            "error": "unknown_accession_type",
            "message": f"Cannot resolve '{accession}'. Expected GSE, SRP, PRJNA, or SRR.",
        }
    try:
        if upper.startswith("GSE"):
            runs = _resolve_gse(accession)
        elif upper.startswith("SRR"):
            runs = _resolve_runs([accession])
        else:
            runs = _resolve_study(accession)
    except LookupFailed as exc:
        return {
            "error": "lookup_failed",
            "message": (
                f"Could not query NCBI/ENA for '{accession}' ({exc}). This is usually "
                "temporary (rate limiting or network) and says nothing about whether the "
                "data exists. Tell the user to wait a few minutes and retry."
            ),
        }

    if not runs:
        return {
            "error": "no_runs_found",
            "message": (
                f"NCBI/ENA were reachable but '{accession}' has no downloadable runs. The "
                "series may have no linked SRA data (e.g. microarray, or only processed "
                "files), be embargoed, require dbGaP authorization, or the accession may "
                "be invalid. Ask the user to verify the accession and access requirements."
            ),
        }

    paths = SESSION.require_paths()
    paths.download_metadata.write_text(
        json.dumps({"accession": accession, "runs": runs}, indent=2) + "\n"
    )

    total_bytes = sum(r.get("total_bytes", 0) for r in runs)
    total_files = sum(len(r.get("files", [])) for r in runs)

    return {
        "accession": accession,
        "n_runs": len(runs),
        "n_files": total_files,
        "total_gb": round(total_bytes / (1024**3), 2) if total_bytes else 0,
        "runs": [
            {
                "run_accession": r["run_accession"],
                "sample": r.get("sample_alias", ""),
                "title": r.get("sample_title", ""),
                "library_layout": r.get("library_layout", ""),
                "n_files": len(r.get("files", [])),
                "gb": round(r.get("total_bytes", 0) / (1024**3), 2),
            }
            for r in runs
        ],
        "metadata_saved": str(paths.download_metadata),
    }


def _select_runs(metadata: dict, runs: list[str] | None) -> tuple[list[dict], str | None]:
    """Resolve a list of run accessions against the metadata. None means all runs."""
    if runs is None:
        return metadata["runs"], None
    if not runs:
        return [], "runs is empty — omit it to use all runs."
    by_acc = {r["run_accession"]: r for r in metadata["runs"]}
    unknown = [r for r in runs if r not in by_acc]
    if unknown:
        return [], f"Not runs of {metadata.get('accession', 'this accession')}: {unknown[:10]}"
    return [by_acc[r] for r in dict.fromkeys(runs)], None


def selected_runs(metadata: dict) -> list[dict]:
    """The runs chosen by generate_download_script (all runs if none were chosen)."""
    return _select_runs(metadata, metadata.get("selected_runs"))[0]


def check_existing_files(search_dir: str, runs: list[str] | None = None) -> Summary:
    """Check if expected FASTQ files already exist and validate MD5 checksums.

    Reads expected files from download_metadata.json (saved by resolve_accession).
    Computes MD5 for any files found and compares to expected checksums.
    runs optionally restricts the check to those run accessions.
    """
    from core.session import SESSION

    paths = SESSION.require_paths()
    if not paths.download_metadata.is_file():
        return {"error": "no_metadata", "message": "Run resolve_accession first."}

    metadata = json.loads(paths.download_metadata.read_text())
    runs, err = _select_runs(metadata, runs)
    if err:
        return {"error": "bad_runs", "message": err}
    dirpath = Path(search_dir)

    results = []
    n_valid = 0
    n_missing = 0
    n_corrupted = 0

    for run in runs:
        for finfo in run.get("files", []):
            filename = finfo["filename"]
            expected_md5 = finfo.get("md5", "")
            filepath = dirpath / filename

            if not filepath.is_file():
                results.append(
                    {"filename": filename, "run": run["run_accession"], "status": "missing"}
                )
                n_missing += 1
            elif expected_md5:
                actual_md5 = _compute_md5(filepath)
                if actual_md5 == expected_md5:
                    results.append(
                        {"filename": filename, "run": run["run_accession"], "status": "valid"}
                    )
                    n_valid += 1
                else:
                    results.append(
                        {"filename": filename, "run": run["run_accession"], "status": "corrupted"}
                    )
                    n_corrupted += 1
            else:
                results.append(
                    {"filename": filename, "run": run["run_accession"], "status": "exists_unverified"}
                )
                n_valid += 1

    return {
        "search_dir": str(dirpath),
        "n_files": len(results),
        "n_valid": n_valid,
        "n_missing": n_missing,
        "n_corrupted": n_corrupted,
        "all_valid": n_missing == 0 and n_corrupted == 0,
        "files": results,
    }


def generate_download_script(output_dir: str, runs: list[str] | None = None) -> Summary:
    """Generate a shell script to download FASTQ files from ENA.

    Reads download metadata from disk. runs optionally restricts the download to
    those run accessions (default: all). Checks output_dir for files that already
    exist and pass MD5 — those are skipped. Writes download.sh and records the
    output directory and selected runs in download_metadata.json.
    """
    from core.session import SESSION

    paths = SESSION.require_paths()
    if not paths.download_metadata.is_file():
        return {"error": "no_metadata", "message": "Run resolve_accession first."}

    metadata = json.loads(paths.download_metadata.read_text())
    selection = runs
    runs, err = _select_runs(metadata, runs)
    if err:
        return {"error": "bad_runs", "message": err}
    dirpath = Path(output_dir).resolve()
    # Record the selection up front so resume and validation see it even when
    # nothing needs downloading.
    metadata["selected_runs"] = [r["run_accession"] for r in runs] if selection is not None else None
    metadata["output_dir"] = str(dirpath)
    paths.download_metadata.write_text(json.dumps(metadata, indent=2) + "\n")

    skip = set()
    if dirpath.is_dir():
        for run in runs:
            for finfo in run.get("files", []):
                fp = dirpath / finfo["filename"]
                if fp.is_file() and finfo.get("md5") and _compute_md5(fp) == finfo["md5"]:
                    skip.add(finfo["filename"])

    downloads = []
    for run in runs:
        for finfo in run.get("files", []):
            if finfo["filename"] not in skip:
                downloads.append(finfo)

    if not downloads:
        return {
            "n_files": 0,
            "n_runs_selected": len(runs),
            "n_skipped": len(skip),
            "message": "All files already present and valid. No downloads needed.",
            "script_path": None,
        }

    has_aspera = any(dl.get("aspera") for dl in downloads)

    lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        f"# Download {len(downloads)} FASTQ file(s) from ENA",
        f"# Output: {dirpath}",
        "",
        f'OUTPUT_DIR="{dirpath}"',
        'mkdir -p "$OUTPUT_DIR"',
        'cd "$OUTPUT_DIR"',
        "",
        "# --- Detect fastest available download tool ---",
    ]

    if has_aspera:
        lines += [
            "ASPERA_KEY=${ASPERA_KEY:-${HOME}/.aspera/connect/etc/asperaweb_id_dsa.openssh}",
            'if command -v ascp &>/dev/null && [ -f "$ASPERA_KEY" ]; then',
            '  DL_METHOD="aspera"',
            '  echo "Using Aspera (fastest)"',
            'elif command -v aria2c &>/dev/null; then',
            '  DL_METHOD="aria2c"',
            '  echo "Using aria2c (multi-connection)"',
            "else",
            '  DL_METHOD="curl"',
            '  echo "Using curl (install aria2c or Aspera for faster downloads)"',
            "fi",
        ]
    else:
        lines += [
            'if command -v aria2c &>/dev/null; then',
            '  DL_METHOD="aria2c"',
            '  echo "Using aria2c (multi-connection)"',
            "else",
            '  DL_METHOD="curl"',
            '  echo "Using curl (install aria2c for faster downloads)"',
            "fi",
        ]

    lines.append("")

    for dl in downloads:
        lines.append(f'echo "Downloading {dl["filename"]} ..."')
        if has_aspera and dl.get("aspera"):
            lines += [
                'if [ "$DL_METHOD" = "aspera" ]; then',
                f'  ascp -QT -l 300m -P33001 -i "$ASPERA_KEY" \\\n'
                f'    {dl["aspera"]} ./',
                'elif [ "$DL_METHOD" = "aria2c" ]; then',
                f'  aria2c -x 4 -s 4 --retry-wait=5 -m 3 -q \\\n'
                f'    -o "{dl["filename"]}" "{dl["url"]}"',
                "else",
                f'  curl -fSL --retry 3 --retry-delay 5 \\\n'
                f'    -o "{dl["filename"]}" "{dl["url"]}"',
                "fi",
            ]
        else:
            lines += [
                'if [ "$DL_METHOD" = "aria2c" ]; then',
                f'  aria2c -x 4 -s 4 --retry-wait=5 -m 3 -q \\\n'
                f'    -o "{dl["filename"]}" "{dl["url"]}"',
                "else",
                f'  curl -fSL --retry 3 --retry-delay 5 \\\n'
                f'    -o "{dl["filename"]}" "{dl["url"]}"',
                "fi",
            ]
        lines.append("")

    md5_checks = [dl for dl in downloads if dl.get("md5")]
    if md5_checks:
        lines += [
            "# --- MD5 verification ---",
            'echo ""',
            'echo "Verifying checksums..."',
            "",
            "if command -v md5sum &>/dev/null; then",
            '  md5check() { md5sum "$1" | cut -d" " -f1; }',
            "elif command -v md5 &>/dev/null; then",
            '  md5check() { md5 -q "$1"; }',
            "else",
            '  echo "Warning: no md5 tool found. Skipping in-script verification."',
            '  echo "(Post-download validation will still run.)"',
            "  exit 0",
            "fi",
            "",
            "FAILED=0",
        ]
        for dl in md5_checks:
            lines += [
                f'ACTUAL=$(md5check "{dl["filename"]}")',
                f'if [ "$ACTUAL" != "{dl["md5"]}" ]; then',
                f'  echo "FAIL: {dl["filename"]}"',
                "  FAILED=$((FAILED + 1))",
                "else",
                f'  echo "OK: {dl["filename"]}"',
                "fi",
            ]
        lines += [
            "",
            'if [ "$FAILED" -gt 0 ]; then',
            '  echo "$FAILED file(s) failed checksum verification."',
            "  exit 1",
            "fi",
            'echo "All checksums verified."',
        ]

    script_content = "\n".join(lines) + "\n"
    paths.download_script.write_text(script_content)
    paths.download_script.chmod(0o755)

    total_bytes = sum(dl.get("bytes", 0) for dl in downloads)
    metadata["output_dir"] = str(dirpath)
    metadata["n_downloads"] = len(downloads)
    metadata["download_bytes"] = total_bytes
    paths.download_metadata.write_text(json.dumps(metadata, indent=2) + "\n")

    return {
        "n_files": len(downloads),
        "n_runs_selected": len(runs),
        "n_runs_total": len(metadata["runs"]),
        "n_skipped": len(skip),
        "total_gb": round(total_bytes / (1024**3), 2) if total_bytes else 0,
        "output_dir": str(dirpath),
        "script_path": str(paths.download_script),
    }


# --- Post-approval functions (called by run.py, not the agent loop) ---


def execute_download(script_path: str) -> Summary:
    """Execute the download script. Called by run.py after human approval."""
    path = Path(script_path)
    if not path.is_file():
        return {"error": "not_found", "message": f"Script not found: {script_path}"}

    result = subprocess.run(["bash", str(path)])

    return {"returncode": result.returncode, "success": result.returncode == 0}


def validate_downloads(runs: list[dict], download_dir: str) -> Summary:
    """Validate MD5 checksums of downloaded files. Strict pass/fail.

    Called by run.py after download execution. The checksums come from
    the metadata file (written by resolve_accession from ENA data),
    never from the LLM.
    """
    dirpath = Path(download_dir)
    results = []
    n_pass = 0
    n_fail = 0
    n_missing = 0

    for run in runs:
        for finfo in run.get("files", []):
            filename = finfo["filename"]
            expected_md5 = finfo.get("md5", "")
            filepath = dirpath / filename

            if not filepath.is_file():
                results.append({"filename": filename, "status": "missing"})
                n_missing += 1
            elif not expected_md5:
                results.append({"filename": filename, "status": "no_checksum"})
            else:
                actual_md5 = _compute_md5(filepath)
                if actual_md5 == expected_md5:
                    results.append({"filename": filename, "status": "pass"})
                    n_pass += 1
                else:
                    results.append(
                        {
                            "filename": filename,
                            "status": "fail",
                            "expected": expected_md5,
                            "actual": actual_md5,
                        }
                    )
                    n_fail += 1

    return {
        "download_dir": str(dirpath),
        "n_files": len(results),
        "n_pass": n_pass,
        "n_fail": n_fail,
        "n_missing": n_missing,
        "all_valid": n_fail == 0 and n_missing == 0,
        "files": results,
    }

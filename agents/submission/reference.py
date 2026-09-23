"""Shared reference cache: reuse the Salmon index nf-core built in an earlier run.

nf-core/rnaseq rebuilds the Salmon index (and re-stages the iGenomes GTF/FASTA) in every
run. With save_reference, a run publishes what it built to results/genome/; after a
successful run those files are copied here, and later Salmon-only runs point nf-core at
them instead of --genome. Only code reads and writes this layout — the LLM never sees
or relays these paths.

Layout (under data/reference/, gitignored):
  <genome>/nf-core-rnaseq-<revision>/
    genes.gtf         GTF nf-core used (the filtered one, if it filtered)
    transcripts.fa    transcript FASTA the index was built from
    salmon_index/     Salmon index
    reference.json    provenance: source run, versions, MD5s

Keyed by pipeline revision because each release pins its own Salmon version, and an
index must be read by the Salmon that built it.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any

from core import config

# Files nf-core writes into every Salmon index; a directory without them is incomplete.
_SALMON_INDEX_FILES = ("info.json", "versionInfo.json")


def reference_dir(genome: str, revision: str) -> Path:
    return config.REFERENCE_DIR / genome / f"nf-core-rnaseq-{revision}"


def salmon_only(extra_args: dict[str, Any]) -> bool:
    """The cache only covers Salmon-only runs (no genome FASTA or aligner index is cached)."""
    return bool(extra_args.get("skip_alignment")) and extra_args.get("pseudo_aligner") == "salmon"


def find_reference(genome: str, revision: str) -> Path | None:
    """The cached reference directory for this genome and revision, if complete."""
    d = reference_dir(genome, revision)
    complete = (
        (d / "reference.json").is_file()
        and (d / "genes.gtf").is_file()
        and (d / "transcripts.fa").is_file()
        and all((d / "salmon_index" / f).is_file() for f in _SALMON_INDEX_FILES)
    )
    return d if complete else None


def nf_reference_params(ref_dir: Path) -> dict[str, str]:
    """nf-core params that replace --genome when using a cached reference."""
    return {
        "gtf": str(ref_dir / "genes.gtf"),
        "transcript_fasta": str(ref_dir / "transcripts.fa"),
        "salmon_index": str(ref_dir / "salmon_index"),
    }


def cache_reference(genome: str, revision: str, results_dir: Path, run_name: str) -> dict[str, Any]:
    """Copy the reference a run published (save_reference) into the shared cache.

    Refuses anything ambiguous or incomplete. Copies into a temporary directory and
    renames it into place only when everything is there, so an interrupted copy never
    looks like a usable reference.
    """
    dest = reference_dir(genome, revision)
    if find_reference(genome, revision):
        return {"status": "already_cached", "path": str(dest)}

    genome_dir = results_dir / "genome"
    gtfs = sorted(genome_dir.glob("*.gtf"))
    filtered = [g for g in gtfs if g.name.endswith(".filtered.gtf")]
    gtf = filtered[0] if len(filtered) == 1 else gtfs[0] if len(gtfs) == 1 else None
    transcripts = sorted(genome_dir.glob("*transcripts.fa"))
    index = genome_dir / "index" / "salmon"

    problems = []
    if gtf is None:
        problems.append(f"expected one GTF in {genome_dir}, found {[g.name for g in gtfs]}")
    if len(transcripts) != 1:
        problems.append(f"expected one *transcripts.fa in {genome_dir}, found {[t.name for t in transcripts]}")
    missing = [f for f in _SALMON_INDEX_FILES if not (index / f).is_file()]
    if missing:
        problems.append(f"Salmon index {index} is missing {', '.join(missing)}")
    if dest.exists():
        problems.append(f"{dest} exists but is incomplete; remove it to re-cache")
    if problems:
        return {"status": "not_cached", "problems": problems}

    tmp = dest.parent / f".{dest.name}.tmp-{os.getpid()}"
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    try:
        shutil.copy2(gtf, tmp / "genes.gtf")
        shutil.copy2(transcripts[0], tmp / "transcripts.fa")
        shutil.copytree(index, tmp / "salmon_index")
        provenance = {
            "genome": genome,
            "source": f"iGenomes '{genome}' via nf-core --genome, processed by nf-core/rnaseq {revision}",
            "pipeline": config.NFCORE_PIPELINE,
            "revision": revision,
            "salmon_version": _salmon_version(results_dir),
            "built_by_run": run_name,
            "copied_from": {"gtf": gtf.name, "transcripts_fa": transcripts[0].name, "salmon_index": "index/salmon"},
            **_input_sources(results_dir),
            "md5": {"genes.gtf": _md5(tmp / "genes.gtf"), "transcripts.fa": _md5(tmp / "transcripts.fa")},
            "created": datetime.now().isoformat(timespec="seconds"),
        }
        (tmp / "reference.json").write_text(json.dumps(provenance, indent=2) + "\n")
        tmp.rename(dest)
    except BaseException:
        shutil.rmtree(tmp, ignore_errors=True)
        raise
    return {"status": "cached", "path": str(dest)}


def _input_sources(results_dir: Path) -> dict[str, str]:
    """The GTF and FASTA nf-core was given (e.g. iGenomes URLs), from its params file, so
    runs reusing this cache can still report which annotation it came from."""
    params_files = sorted((results_dir / "pipeline_info").glob("params_*.json"))
    if not params_files:
        return {}
    params = json.loads(params_files[-1].read_text())
    return {f"{key}_source": str(params[key]) for key in ("gtf", "fasta") if params.get(key)}


def _salmon_version(results_dir: Path) -> str | None:
    versions = results_dir / "pipeline_info" / "nf_core_rnaseq_software_mqc_versions.yml"
    if not versions.is_file():
        return None
    try:
        import yaml
        data = yaml.safe_load(versions.read_text()) or {}
    except Exception:
        return None
    for process, tools in data.items():
        if "SALMON_INDEX" in str(process) and isinstance(tools, dict) and "salmon" in tools:
            return str(tools["salmon"])
    return None


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

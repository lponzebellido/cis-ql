#!/usr/bin/env python3

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import platform
import re
import statistics
import subprocess
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "benchmarks" / "examples.json"
LOAD_PATTERN = re.compile(
    r'\bLOAD\s+(?:SEQUENCE|ANNOTATION|MATRIX|TRACK)\s+"([^"]+)"'
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_metadata(path: Path, display_path: str | None = None) -> dict:
    return {
        "path": display_path or str(path),
        "bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def validate_relative_path(value: str, label: str) -> Path:
    path = Path(value)
    if path.is_absolute() or ".." in path.parts:
        raise RuntimeError(f"{label} must stay inside its workspace: {value}")
    return path


def load_manifest(path: Path) -> list[dict]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schemaVersion") != 1:
        raise RuntimeError(f"Unsupported benchmark manifest: {path}")
    examples = payload.get("examples")
    if not isinstance(examples, list) or not examples:
        raise RuntimeError(f"Benchmark manifest has no examples: {path}")
    identifiers = [entry.get("id") for entry in examples]
    files = [entry.get("file") for entry in examples]
    if len(identifiers) != len(set(identifiers)):
        raise RuntimeError("Benchmark manifest contains duplicate IDs")
    if len(files) != len(set(files)):
        raise RuntimeError("Benchmark manifest contains duplicate files")
    for entry in examples:
        required = {"id", "file", "exports", "resultCounts"}
        if not required.issubset(entry):
            missing = sorted(required - set(entry))
            raise RuntimeError(
                f"Benchmark manifest entry is missing: {', '.join(missing)}"
            )
        if not isinstance(entry["exports"], list):
            raise RuntimeError(f"Invalid exports for example {entry['id']}")
        example_relative = validate_relative_path(
            entry["file"], "Example path"
        )
        for exported in entry["exports"]:
            validate_relative_path(exported, "Export path")
        example = ROOT / "cql_examples" / example_relative
        if not example.is_file():
            raise RuntimeError(f"Missing example declared in manifest: {example}")
    declared = set(files)
    discovered = {
        example.name for example in (ROOT / "cql_examples").glob("*.cql")
    }
    if declared != discovered:
        missing = sorted(discovered - declared)
        extra = sorted(declared - discovered)
        raise RuntimeError(
            f"Example manifest mismatch; unlisted={missing}, missing={extra}"
        )
    return examples


def select_examples(entries: list[dict], selectors: list[str] | None) -> list[dict]:
    if not selectors:
        return entries
    selected = []
    seen = set()
    for selector in selectors:
        matches = [
            entry for entry in entries
            if selector in {
                entry["id"], entry["file"], Path(entry["file"]).stem
            }
        ]
        if len(matches) != 1:
            raise RuntimeError(
                f"Expected one manifest entry matching '{selector}', "
                f"found {len(matches)}"
            )
        identifier = matches[0]["id"]
        if identifier not in seen:
            selected.append(matches[0])
            seen.add(identifier)
    return selected


def summarize_results(parsed: dict) -> dict:
    summary = {
        section: {
            name: len(values)
            for name, values in parsed.get(section, {}).items()
        }
        for section in ("resultSets", "gcProfiles")
        if parsed.get(section)
    }
    feature_groups = parsed.get("featureGroups", {})
    if feature_groups:
        summary["featureGroups"] = {
            name: {
                "groups": len(groups),
                "members": sum(group["memberCount"] for group in groups),
            }
            for name, groups in feature_groups.items()
        }
    return summary


def resolve_repository_path(relative: str) -> Path:
    path = (ROOT / validate_relative_path(relative, "Input path")).resolve()
    try:
        path.relative_to(ROOT.resolve())
    except ValueError as error:
        raise RuntimeError(f"Input leaves repository: {relative}") from error
    if not path.is_file():
        raise RuntimeError(f"Missing example input: {relative}")
    return path


def prepare_workspace(entry: dict, workspace: Path) -> list[dict]:
    example = ROOT / "cql_examples" / entry["file"]
    relative_inputs = list(dict.fromkeys(
        LOAD_PATTERN.findall(example.read_text(encoding="utf-8"))
    ))
    inputs = []
    for relative in relative_inputs:
        source = resolve_repository_path(relative)
        destination = workspace / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.symlink_to(source)
        inputs.append(file_metadata(source, relative))
    for relative in entry["exports"]:
        (workspace / relative).parent.mkdir(parents=True, exist_ok=True)
    return inputs


def execute_once(
    binary: Path, entry: dict, workspace: Path
) -> tuple[float, dict, dict, list[dict]]:
    result_path = workspace / ".cisql_results.json"
    paths = [result_path] + [workspace / path for path in entry["exports"]]
    for path in paths:
        if path.exists():
            path.unlink()
    example = ROOT / "cql_examples" / entry["file"]
    started = time.perf_counter()
    completed = subprocess.run(
        [str(binary), str(example)],
        cwd=workspace,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    elapsed = time.perf_counter() - started
    if completed.returncode != 0:
        raise RuntimeError(
            f"{entry['file']} failed: {completed.stderr.strip()}"
        )
    if not result_path.is_file():
        raise RuntimeError(f"{entry['file']} produced no results JSON")
    parsed = json.loads(result_path.read_text(encoding="utf-8"))
    summary = summarize_results(parsed)
    if summary != entry["resultCounts"]:
        raise RuntimeError(
            f"{entry['file']} result counts differ from the manifest: "
            f"expected {entry['resultCounts']}, observed {summary}"
        )
    exports = []
    for relative in entry["exports"]:
        path = workspace / relative
        if not path.is_file():
            raise RuntimeError(
                f"{entry['file']} did not create expected export: {relative}"
            )
        exports.append(file_metadata(path, relative))
    fingerprint = {
        "results": sha256_file(result_path),
        "exports": {item["path"]: item["sha256"] for item in exports},
    }
    return elapsed, summary, fingerprint, exports


def benchmark(
    binary: Path, entry: dict, repetitions: int, warmups: int
) -> dict:
    samples = []
    expected_fingerprint = None
    final_exports = []
    with tempfile.TemporaryDirectory(prefix="cisql-benchmark-") as temp:
        workspace = Path(temp)
        inputs = prepare_workspace(entry, workspace)
        for _ in range(warmups):
            execute_once(binary, entry, workspace)
        for _ in range(repetitions):
            elapsed, summary, fingerprint, exports = execute_once(
                binary, entry, workspace
            )
            if expected_fingerprint is None:
                expected_fingerprint = fingerprint
            elif fingerprint != expected_fingerprint:
                raise RuntimeError(
                    f"{entry['file']} produced nondeterministic artifacts"
                )
            samples.append(elapsed)
            final_exports = exports
    return {
        "id": entry["id"],
        "example": entry["file"],
        "inputs": inputs,
        "warmups": warmups,
        "repetitions": repetitions,
        "median_seconds": statistics.median(samples),
        "min_seconds": min(samples),
        "max_seconds": max(samples),
        "samples_seconds": samples,
        "result_counts": summary,
        "result_sha256": expected_fingerprint["results"],
        "exports": final_exports,
    }


def repository_state() -> dict:
    revision = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    status = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    return {
        "revision": revision.stdout.strip() if revision.returncode == 0 else None,
        "dirty": bool(status.stdout.strip()) if status.returncode == 0 else None,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, default=ROOT / "cisql")
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--examples", nargs="*")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    if args.warmups < 0:
        raise RuntimeError("Warm-up count cannot be negative")
    if args.repetitions < 1:
        raise RuntimeError("Repetition count must be positive")
    binary = args.binary.resolve()
    if not binary.is_file():
        raise RuntimeError(f"Missing Cis-QL binary: {binary}")
    manifest = args.manifest.resolve()
    entries = select_examples(load_manifest(manifest), args.examples)
    results = [
        benchmark(binary, entry, args.repetitions, args.warmups)
        for entry in entries
    ]

    payload = {
        "schema_version": 2,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "repository": repository_state(),
        "system": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "processor": platform.processor(),
            "python": platform.python_version(),
        },
        "binary": file_metadata(binary),
        "manifest": file_metadata(manifest),
        "results": results,
    }
    rendered = json.dumps(payload, indent=2)
    if args.output:
        args.output.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

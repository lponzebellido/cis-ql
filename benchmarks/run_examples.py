#!/usr/bin/env python3

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import platform
import statistics
import subprocess
import tempfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EXAMPLES = (
    "03", "04", "05", "06", "07", "09", "12", "13", "14", "15", "16"
)


def benchmark(binary: Path, example: Path, repetitions: int) -> dict:
    samples: list[float] = []
    expected_summary: dict[str, dict[str, int]] | None = None
    with tempfile.TemporaryDirectory(prefix="cisql-benchmark-") as temp:
        workspace = Path(temp)
        (workspace / "data_examples").symlink_to(ROOT / "data_examples",
                                                  target_is_directory=True)
        (workspace / "matrices").symlink_to(ROOT / "matrices",
                                            target_is_directory=True)
        for _ in range(repetitions):
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
                    f"{example.name} failed: {completed.stderr.strip()}"
                )
            results_path = workspace / ".cisql_results.json"
            if not results_path.exists():
                raise RuntimeError(f"{example.name} produced no results JSON")
            parsed = json.loads(results_path.read_text(encoding="utf-8"))
            summary = {
                section: {
                    name: len(values)
                    for name, values in parsed.get(section, {}).items()
                }
                for section in ("resultSets", "gcProfiles")
                if parsed.get(section)
            }
            if expected_summary is None:
                expected_summary = summary
            elif summary != expected_summary:
                raise RuntimeError(
                    f"{example.name} produced inconsistent results: "
                    f"{expected_summary} vs. {summary}"
                )
            samples.append(elapsed)
    return {
        "example": example.name,
        "repetitions": repetitions,
        "median_seconds": statistics.median(samples),
        "min_seconds": min(samples),
        "max_seconds": max(samples),
        "samples_seconds": samples,
        "result_counts": expected_summary,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, default=ROOT / "cisql")
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--examples", nargs="*", default=DEFAULT_EXAMPLES)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    results = []
    for prefix in args.examples:
        matches = sorted((ROOT / "cql_examples").glob(f"{prefix}_*.cql"))
        if len(matches) != 1:
            raise RuntimeError(f"Expected one example matching {prefix}")
        results.append(benchmark(args.binary.resolve(), matches[0],
                                 args.repetitions))

    payload = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "system": {
            "platform": platform.platform(),
            "processor": platform.processor(),
            "python": platform.python_version(),
        },
        "binary": str(args.binary.resolve()),
        "results": results,
    }
    rendered = json.dumps(payload, indent=2)
    if args.output:
        args.output.write_text(rendered + "\n", encoding="utf-8")
    print(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

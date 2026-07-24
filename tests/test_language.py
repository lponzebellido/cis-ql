#!/usr/bin/env python3
"""End-to-end correctness tests for the Cis-QL language."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BINARY = Path(os.environ.get("CISQL_BINARY", ROOT / "cisql"))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def run_query(workspace: Path, name: str, source: str) -> tuple[dict, str]:
    query = workspace / f"{name}.cql"
    query.write_text(source, encoding="utf-8")
    completed = subprocess.run(
        [str(BINARY), str(query)],
        cwd=workspace,
        text=True,
        capture_output=True,
        check=False,
    )
    require(
        completed.returncode == 0,
        f"{name}: return code {completed.returncode}\n"
        f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}",
    )
    result_path = workspace / ".cisql_results.json"
    require(result_path.exists(), f"{name}: missing JSON output")
    return json.loads(result_path.read_text(encoding="utf-8")), completed.stderr


def run_invalid_query(
    workspace: Path, name: str, source: str, expected_code: int
) -> str:
    query = workspace / f"{name}.cql"
    query.write_text(source, encoding="utf-8")
    completed = subprocess.run(
        [str(BINARY), str(query)],
        cwd=workspace,
        text=True,
        capture_output=True,
        check=False,
    )
    require(
        completed.returncode == expected_code,
        f"{name}: expected return code {expected_code}, got "
        f"{completed.returncode}\nstdout:\n{completed.stdout}\n"
        f"stderr:\n{completed.stderr}",
    )
    require(
        not (workspace / ".cisql_results.json").exists(),
        f"{name}: failed execution left a stale results file",
    )
    return completed.stderr


def main() -> int:
    require(BINARY.exists(), "Build cisql before running language tests")
    with tempfile.TemporaryDirectory(prefix="cisql-language-tests-") as temp:
        workspace = Path(temp)
        sequence = "A" * 1000 + "C" * 1000 + ("CG" * 1100)
        (workspace / "fixture.fasta").write_text(
            f">chr1\n{sequence}\n", encoding="utf-8"
        )
        (workspace / "fixture.gff3").write_text(
            "##gff-version 3\n"
            "chr1\ttest\tgene\t1\t1200\t.\t+\t.\tID=short\n"
            "chr1\ttest\tgene\t1301\t2900\t.\t+\t.\tID=long\n"
            "chr1\ttest\tgene\t3001\t3400\t.\t+\t.\tID=tiny\n",
            encoding="utf-8",
        )
        (workspace / "fixture.pwm").write_text(
            ">TEST test\n"
            "A [ 10 10 ]\n"
            "C [ 0 0 ]\n"
            "G [ 0 0 ]\n"
            "T [ 0 0 ]\n",
            encoding="utf-8",
        )

        data, _ = run_query(
            workspace,
            "decimal_units",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            "EXTRACT GENE WHERE LENGTH > 1.5 KB;\n",
        )
        extracted = next(
            values
            for key, values in data["resultSets"].items()
            if key.startswith("Extract_GENE")
        )
        require(len(extracted) == 1 and extracted[0]["name"] == "long",
                "decimal KB conversion")

        data, _ = run_query(
            workspace,
            "logical_filters",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            "EXTRACT GENE WHERE LENGTH < 500 BP OR NOT LENGTH < 1500 BP;\n",
        )
        extracted = next(
            values
            for key, values in data["resultSets"].items()
            if key.startswith("Extract_GENE")
        )
        require({item["name"] for item in extracted} == {"tiny", "long"},
                "OR and NOT semantics")

        data, _ = run_query(
            workspace,
            "similarity_reference",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            "EXTRACT GENE WHERE LENGTH > 1.5 KB "
            "AND SIMILARITY > 99 %;\n",
        )
        extracted = next(
            values
            for key, values in data["resultSets"].items()
            if key.startswith("Extract_GENE")
        )
        require(len(extracted) == 1 and extracted[0]["name"] == "long",
                "similarity reference follows non-similarity predicates")

        data, _ = run_query(
            workspace,
            "motif_length",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'FIND MOTIF "A{50}|C{150}" STRAND POSITIVE AS long_runs '
            'WHERE LENGTH > 100 BP;\n',
        )
        require(len(data["resultSets"]["long_runs"]) > 0 and
                all(item["end"] - item["start"] == 150
                    for item in data["resultSets"]["long_runs"]),
                "motif length filtering")

        data, _ = run_query(
            workspace,
            "analyze_alias",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            "ANALYZE CPG_ISLANDS AS islands;\n",
        )
        require(len(data["resultSets"]["islands"]) > 0,
                "ANALYZE alias propagation")

        data, _ = run_query(
            workspace,
            "gc_profile_filter",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            "ANALYZE GC_CONTENT WINDOW 1 KB AS high_gc "
            "WHERE GC_CONTENT > 50 %;\n",
        )
        require(set(data["gcProfiles"]) == {"high_gc"} and
                len(data["gcProfiles"]["high_gc"]) == 3 and
                all(window["gc"] > 50
                    for window in data["gcProfiles"]["high_gc"]),
                "GC profile filtering and alias propagation")

        data, _ = run_query(
            workspace,
            "decimal_window",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            "ANALYZE GC_CONTENT WINDOW 1.5 KB AS gc_windows;\n",
        )
        require(len(data["gcProfiles"]["gc_windows"]) == 2,
                "decimal ANALYZE window conversion")

        data, _ = run_query(
            workspace,
            "motif_gc_filter",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'FIND MOTIF "NN" STRAND POSITIVE AS gc_dinucleotides '
            "WHERE GC_CONTENT > 50 %;\n",
        )
        require(len(data["resultSets"]["gc_dinucleotides"]) > 0 and
                all(sum(base in "GC" for base in item["sequence"]) > 1
                    for item in data["resultSets"]["gc_dinucleotides"]),
                "motif GC_CONTENT filtering")

        data, _ = run_query(
            workspace,
            "foreach_aliases",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD MATRIX "fixture.pwm" AS first;\n'
            'LOAD MATRIX "fixture.pwm" AS second;\n'
            "FOREACH m IN [first, second] DO\n"
            "  SCAN m STRAND POSITIVE THRESHOLD 90 % AS sites;\n"
            "ENDFOR;\n",
        )
        require("sites_first" in data["resultSets"] and
                "sites_second" in data["resultSets"],
                "FOREACH preserves per-iteration result tracks")

        data, _ = run_query(
            workspace,
            "zero_pwm_threshold",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            "SCAN matrix STRAND POSITIVE THRESHOLD 0 % AS all_windows;\n",
        )
        require(len(data["resultSets"]["all_windows"]) == len(sequence) - 1,
                "explicit zero PWM threshold is not replaced by the default")

        data, _ = run_query(
            workspace,
            "complex_if",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            "IF GC_CONTENT > 50 % AND NOT GC_CONTENT < 50 % THEN\n"
            '  FIND MOTIF "CG" AS high_gc;\n'
            "ELSE\n"
            '  FIND MOTIF "TT" AS low_gc;\n'
            "ENDIF;\n",
        )
        require("high_gc" in data["resultSets"] and
                "low_gc" not in data["resultSets"],
                "compound IF condition")

        lexical_error = run_invalid_query(
            workspace,
            "lexical_error",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n@\n',
            2,
        )
        require("Lexical Error" in lexical_error,
                "lexical failures return a non-zero status")

        semantic_error = run_invalid_query(
            workspace,
            "wrong_alias_type",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            "SCAN genome THRESHOLD 80 %;\n",
            3,
        )
        require("expects a matrix alias" in semantic_error,
                "semantic alias type checking")

        (workspace / "mismatch.gff3").write_text(
            "##gff-version 3\n"
            "chrMissing\ttest\tgene\t1\t10\t.\t+\t.\tID=mismatch\n",
            encoding="utf-8",
        )
        runtime_error = run_invalid_query(
            workspace,
            "chromosome_mismatch",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD ANNOTATION "mismatch.gff3" AS annot;\n'
            'FIND MOTIF "AA" WITHIN 10 BP UPSTREAM FROM GENE;\n',
            4,
        )
        require("not present in the active FASTA" in runtime_error,
                "runtime chromosome compatibility checking")

        invalid_pattern = run_invalid_query(
            workspace,
            "invalid_pattern",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'FIND MOTIF "(";\n',
            4,
        )
        require("Invalid or empty motif pattern" in invalid_pattern,
                "invalid regular expressions are execution errors")

    print("All language tests passed.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AssertionError as error:
        print(f"FAILED: {error}", file=sys.stderr)
        raise SystemExit(1)

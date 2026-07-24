#!/usr/bin/env python3
"""Validate Cis-QL outputs against small independent reference implementations."""

from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BINARY = Path(os.environ.get("CISQL_BINARY", ROOT / "cisql"))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def run_query(workspace: Path, name: str, source: str) -> dict:
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
        f"{name} failed ({completed.returncode})\n"
        f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}",
    )
    return json.loads(
        (workspace / ".cisql_results.json").read_text(encoding="utf-8")
    )


def overlapping_exact_matches(sequence: str, motif: str) -> list[int]:
    return [
        start
        for start in range(len(sequence) - len(motif) + 1)
        if sequence[start : start + len(motif)].upper() == motif.upper()
    ]


def subtract_intervals(
    intervals: list[tuple[int, int]], masks: list[tuple[int, int]]
) -> list[tuple[int, int]]:
    result: list[tuple[int, int]] = []
    for start, end in intervals:
        fragments = [(start, end)]
        for mask_start, mask_end in masks:
            updated: list[tuple[int, int]] = []
            for frag_start, frag_end in fragments:
                if mask_end <= frag_start or mask_start >= frag_end:
                    updated.append((frag_start, frag_end))
                    continue
                if frag_start < mask_start:
                    updated.append((frag_start, mask_start))
                if mask_end < frag_end:
                    updated.append((mask_end, frag_end))
            fragments = updated
        result.extend(fragments)
    return result


def smith_waterman_similarity(first: str, second: str) -> float:
    previous = [0] * (len(second) + 1)
    best = 0
    for left in first.upper():
        current = [0]
        for index, right in enumerate(second.upper(), start=1):
            current.append(
                max(
                    0,
                    previous[index - 1] + (2 if left == right else -1),
                    previous[index] - 2,
                    current[index - 1] - 2,
                )
            )
            best = max(best, current[-1])
        previous = current
    normalizer = 2 * min(len(first), len(second))
    return 0.0 if normalizer == 0 else 100.0 * best / normalizer


def pssm_percent_scores(
    sequence: str, counts: list[list[float]]
) -> list[float]:
    length = len(counts[0])
    backgrounds = [0.25] * 4
    scores = [[0.0] * length for _ in range(4)]
    minimum = 0.0
    maximum = 0.0
    for position in range(length):
        total = sum(row[position] for row in counts)
        column = [
            math.log2(
                ((counts[base][position] + 0.1) / (total + 0.4))
                / backgrounds[base]
            )
            for base in range(4)
        ]
        minimum += min(column)
        maximum += max(column)
        for base in range(4):
            scores[base][position] = column[base]

    nucleotide = {"A": 0, "C": 1, "G": 2, "T": 3}
    percentages: list[float] = []
    for start in range(len(sequence) - length + 1):
        raw = 0.0
        valid = True
        for offset, base in enumerate(sequence[start : start + length]):
            if base.upper() not in nucleotide:
                valid = False
                break
            raw += scores[nucleotide[base.upper()]][offset]
        percentages.append(
            (raw - minimum) / (maximum - minimum) * 100.0 if valid else 0.0
        )
    return percentages


def validate_optional_bedtools(workspace: Path) -> str:
    bedtools = shutil.which("bedtools")
    if not bedtools:
        return "bedtools not installed (optional comparison skipped)"

    (workspace / "a.bed").write_text(
        "chr1\t0\t10\tgene\t0\t+\n", encoding="utf-8"
    )
    (workspace / "b.bed").write_text(
        "chr1\t3\t6\texon\t0\t+\n", encoding="utf-8"
    )
    completed = subprocess.run(
        [bedtools, "subtract", "-a", "a.bed", "-b", "b.bed"],
        cwd=workspace,
        text=True,
        capture_output=True,
        check=True,
    )
    bedtools_coordinates = [
        tuple(map(int, line.split("\t")[1:3]))
        for line in completed.stdout.splitlines()
    ]
    cisql_coordinates = [
        tuple(map(int, line.split("\t")[1:3]))
        for line in (workspace / "difference.bed")
        .read_text(encoding="utf-8")
        .splitlines()
    ]
    require(
        cisql_coordinates == bedtools_coordinates,
        "Cis-QL interval subtraction differs from bedtools subtract",
    )
    return "bedtools subtract coordinates agree"


def validate_optional_biopython() -> str:
    try:
        from Bio import Align  # type: ignore
    except ImportError:
        return "Biopython not installed (optional comparison skipped)"

    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -2
    first, second = "ACGT", "ACGA"
    expected = smith_waterman_similarity(first, second)
    observed = 100.0 * aligner.score(first, second) / (2 * min(len(first), len(second)))
    require(abs(expected - observed) < 1e-9,
            "independent Smith-Waterman score differs from Biopython")
    return "Biopython local-alignment score agrees"


def main() -> int:
    require(BINARY.exists(), "Build cisql before running reference validation")
    with tempfile.TemporaryDirectory(prefix="cisql-reference-") as temp:
        workspace = Path(temp)
        sequence = "ATATATAAACCGGTT"
        (workspace / "motifs.fasta").write_text(
            f">chr1\n{sequence}\n", encoding="utf-8"
        )
        motif_data = run_query(
            workspace,
            "motif_reference",
            'LOAD SEQUENCE "motifs.fasta" AS genome;\n'
            'FIND MOTIF "ATA" STRAND POSITIVE AS motif_hits;\n',
        )
        observed_motifs = [
            region["start"]
            for region in motif_data["resultSets"]["motif_hits"]
        ]
        require(
            observed_motifs == overlapping_exact_matches(sequence, "ATA"),
            "exact motif coordinates differ from independent reference",
        )
        print("[ok] overlapping exact motif coordinates")

        counts = [[10, 10], [0, 0], [0, 0], [0, 0]]
        (workspace / "aa.pwm").write_text(
            ">AA exact_AA\n"
            "A [ 10 10 ]\nC [ 0 0 ]\nG [ 0 0 ]\nT [ 0 0 ]\n",
            encoding="utf-8",
        )
        pwm_data = run_query(
            workspace,
            "pwm_reference",
            'LOAD SEQUENCE "motifs.fasta" AS genome;\n'
            'LOAD MATRIX "aa.pwm" AS aa;\n'
            "SCAN aa STRAND POSITIVE THRESHOLD 100 % AS pwm_hits;\n",
        )
        observed_pwm = [
            region["start"]
            for region in pwm_data["resultSets"]["pwm_hits"]
        ]
        expected_pwm = [
            index
            for index, score in enumerate(pssm_percent_scores(sequence, counts))
            if score >= 100.0 - 1e-9
        ]
        require(observed_pwm == expected_pwm,
                "PWM coordinates differ from independent log-odds reference")
        print("[ok] PWM/PSSM threshold coordinates")

        (workspace / "intervals.gff3").write_text(
            "##gff-version 3\n"
            "chr1\tref\tgene\t1\t10\t.\t+\t.\tID=gene\n"
            "chr1\tref\texon\t4\t6\t.\t+\t.\tID=exon\n",
            encoding="utf-8",
        )
        interval_data = run_query(
            workspace,
            "interval_reference",
            'LOAD ANNOTATION "intervals.gff3" AS annot;\n'
            "EXTRACT GENE AS genes;\n"
            "EXTRACT EXON AS exons;\n"
            "EXCEPT genes FROM exons AS difference;\n"
            'EXPORT difference TO "difference.bed" FORMAT BED;\n',
        )
        observed_intervals = [
            (region["start"], region["end"])
            for region in interval_data["resultSets"]["difference"]
        ]
        require(
            observed_intervals == subtract_intervals([(0, 10)], [(3, 6)]),
            "EXCEPT differs from independent geometric subtraction",
        )
        print("[ok] geometric interval subtraction")

        (workspace / "similarity.fasta").write_text(
            ">chr1\nACGTNNACGA\n", encoding="utf-8"
        )
        (workspace / "similarity.gff3").write_text(
            "##gff-version 3\n"
            "chr1\tref\tgene\t1\t4\t.\t+\t.\tID=reference\n"
            "chr1\tref\tgene\t7\t10\t.\t+\t.\tID=candidate\n",
            encoding="utf-8",
        )
        similarity = smith_waterman_similarity("ACGT", "ACGA")
        similarity_data = run_query(
            workspace,
            "similarity_reference",
            'LOAD SEQUENCE "similarity.fasta" AS genome;\n'
            'LOAD ANNOTATION "similarity.gff3" AS annot;\n'
            'EXTRACT GENE AS reference WHERE ID = "reference";\n'
            f"EXTRACT GENE AS similar WHERE SIMILARITY TO reference "
            f">= {similarity:.6f} %;\n",
        )
        require(
            {region["name"] for region
             in similarity_data["resultSets"]["similar"]}
            == {"reference", "candidate"},
            "explicit similarity filtering differs from independent reference",
        )
        print("[ok] normalized Smith-Waterman filter")

        print(f"[optional] {validate_optional_bedtools(workspace)}")
        print(f"[optional] {validate_optional_biopython()}")

    print("All independent reference validations passed.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AssertionError as error:
        print(f"FAILED: {error}", file=os.sys.stderr)
        raise SystemExit(1)

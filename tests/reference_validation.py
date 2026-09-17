#!/usr/bin/env python3

from __future__ import annotations

import json
import itertools
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


def select_overlapping_intervals(
    query: list[tuple[str, int, int]],
    reference: list[tuple[str, int, int]],
) -> list[tuple[str, int, int]]:
    return [
        interval
        for interval in query
        if any(
            interval[0] == candidate[0]
            and interval[1] < candidate[2]
            and interval[2] > candidate[1]
            for candidate in reference
        )
    ]


def interval_gap(
    query: tuple[str, int, int], reference: tuple[str, int, int]
) -> int:
    if query[1] < reference[2] and query[2] > reference[1]:
        return 0
    if query[2] <= reference[1]:
        return reference[1] - query[2]
    return query[1] - reference[2]


def select_nearest_intervals(
    query: list[tuple[str, int, int]],
    reference: list[tuple[str, int, int, str]],
    maximum_distance: int,
) -> list[tuple[tuple[str, int, int], tuple[str, int, int, str], int, bool]]:
    selected = []
    for interval in query:
        same_chromosome = [
            candidate for candidate in reference
            if candidate[0] == interval[0]
        ]
        if not same_chromosome:
            continue
        nearest = min(
            same_chromosome,
            key=lambda candidate: (
                interval_gap(interval, candidate[:3]), candidate[1],
                candidate[2], candidate[3],
            ),
        )
        distance = interval_gap(interval, nearest[:3])
        if distance <= maximum_distance:
            overlaps = interval[1] < nearest[2] and interval[2] > nearest[1]
            selected.append((interval, nearest, distance, overlaps))
    return selected


def count_overlaps_by_container(
    counted: list[tuple[str, int, int]],
    containers: list[tuple[str, int, int]],
) -> list[int]:
    return [
        sum(
            candidate[0] == container[0]
            and candidate[1] < container[2]
            and candidate[2] > container[1]
            for candidate in counted
        )
        for container in containers
    ]


def anchor_consensus(
    anchor: list[tuple[str, int, int]],
    supporting_sets: list[list[tuple[str, int, int]]],
    minimum_support: int,
) -> list[tuple[tuple[str, int, int], int]]:
    result = []
    for interval in anchor:
        observed_support = 1 + sum(
            any(
                interval[0] == candidate[0]
                and interval[1] < candidate[2]
                and interval[2] > candidate[1]
                for candidate in support_set
            )
            for support_set in supporting_sets
        )
        if observed_support >= minimum_support:
            result.append((interval, observed_support))
    return result


def define_homotypic_modules(
    sites: list[tuple[str, int, int, str]],
    minimum_spacing: int,
    maximum_spacing: int,
) -> list[tuple[str, int, int, int, str]]:
    modules = []
    for first, second in itertools.combinations(sites, 2):
        if first[0] != second[0]:
            continue
        spacing = interval_gap(first[:3], second[:3])
        if not minimum_spacing <= spacing <= maximum_spacing:
            continue
        orientation = "SAME" if first[3] == second[3] else "OPPOSITE"
        modules.append(
            (
                first[0], min(first[1], second[1]),
                max(first[2], second[2]), spacing, orientation,
            )
        )
    return modules


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
                ((counts[base][position] + 0.1 * backgrounds[base])
                 / (total + 0.1))
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


def exact_pssm_tail_probability(
    counts: list[list[float]], word: str,
    backgrounds: list[float] | None = None,
) -> float:
    background = backgrounds or [0.25] * 4
    length = len(counts[0])
    require(len(word) == length, "reference word width differs from PWM")
    scores = [[0.0] * length for _ in range(4)]
    for position in range(length):
        total = sum(row[position] for row in counts)
        for base in range(4):
            frequency = (
                counts[base][position] + 0.1 * background[base]
            ) / (total + 0.1)
            scores[base][position] = math.log2(
                frequency / background[base]
            )

    nucleotide = {"A": 0, "C": 1, "G": 2, "T": 3}
    target = sum(
        scores[nucleotide[base]][position]
        for position, base in enumerate(word.upper())
    )
    tail = 0.0
    for candidate in itertools.product(range(4), repeat=length):
        score = sum(scores[base][position]
                    for position, base in enumerate(candidate))
        if score >= target - 1e-12:
            probability = math.prod(background[base] for base in candidate)
            tail += probability
    return tail


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

    (workspace / "overlap_a.bed").write_text(
        "chr1\t0\t10\tgene\t0\t+\n"
        "chr1\t10\t15\ttouching_gene\t0\t+\n",
        encoding="utf-8",
    )
    (workspace / "overlap_b.bed").write_text(
        "chr1\t3\t6\tenhancer_a\t0\t+\n"
        "chr1\t8\t10\tenhancer_b\t0\t+\n",
        encoding="utf-8",
    )
    completed = subprocess.run(
        [bedtools, "intersect", "-a", "overlap_a.bed", "-b",
         "overlap_b.bed", "-u"],
        cwd=workspace,
        text=True,
        capture_output=True,
        check=True,
    )
    bedtools_supported = [
        tuple(map(int, line.split("\t")[1:3]))
        for line in completed.stdout.splitlines()
    ]
    cisql_supported = [
        tuple(map(int, line.split("\t")[1:3]))
        for line in (workspace / "supported.bed")
        .read_text(encoding="utf-8")
        .splitlines()
    ]
    require(cisql_supported == bedtools_supported,
            "Cis-QL OVERLAPS differs from bedtools intersect -u")
    return "bedtools subtract and overlap-selection coordinates agree"


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


def validate_optional_fimo(workspace: Path, cisql_pvalue: float) -> str:
    fimo = shutil.which("fimo")
    if not fimo:
        return "FIMO not installed (optional p-value comparison skipped)"

    meme_motif = workspace / "aa.meme"
    meme_motif.write_text(
        "MEME version 4\n\n"
        "ALPHABET= ACGT\n\n"
        "strands: +\n\n"
        "Background letter frequencies\n"
        "A 0.25 C 0.25 G 0.25 T 0.25\n\n"
        "MOTIF AA exact_AA\n"
        "letter-probability matrix: alength= 4 w= 2 nsites= 10 E= 0\n"
        "1 0 0 0\n"
        "1 0 0 0\n",
        encoding="utf-8",
    )
    completed = subprocess.run(
        [fimo, "--text", "--norc", "--bgfile", "--uniform--",
         "--thresh", "1", str(meme_motif),
         str(workspace / "motifs.fasta")],
        cwd=workspace,
        text=True,
        capture_output=True,
        check=False,
    )
    require(
        completed.returncode == 0,
        f"FIMO comparison failed\nstdout:\n{completed.stdout}\n"
        f"stderr:\n{completed.stderr}",
    )
    rows = [line.split("\t") for line in completed.stdout.splitlines()
            if line and not line.startswith("#")]
    require(len(rows) > 1, "FIMO comparison returned no motif sites")
    header = rows[0]
    pvalue_column = header.index("p-value")
    start_column = header.index("start")
    best = [float(row[pvalue_column]) for row in rows[1:]
            if int(row[start_column]) in (7, 8)]
    require(best, "FIMO comparison did not return the expected AA sites")
    require(all(abs(value - cisql_pvalue) < 1e-9 for value in best),
            "Cis-QL PWM p-value differs from FIMO")
    return "FIMO AA p-values agree"


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
            'FIND MOTIF "ATA" STRAND POSITIVE AS motif_hits;\n'
            "DEFINE MODULE FROM motif_hits WITH motif_hits "
            "SPACING 1 BP TO 1 BP ORDER ANY ORIENTATION SAME "
            "AS motif_modules;\n",
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
        expected_modules = define_homotypic_modules(
            [("chr1", start, start + 3, "+")
             for start in observed_motifs],
            1,
            1,
        )
        observed_modules = [
            (
                region["chr"], region["start"], region["end"],
                region["moduleEvidence"]["spacing"]["observed"],
                region["moduleEvidence"]["orientation"]["observed"],
            )
            for region in motif_data["resultSets"]["motif_modules"]
        ]
        require(observed_modules == expected_modules,
                "DEFINE MODULE differs from independent pair enumeration")
        print("[ok] constrained homotypic motif modules")

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
        cisql_pvalue = pwm_data["resultSets"]["pwm_hits"][0][
            "motifEvidence"
        ]["statistics"]["pValue"]
        exact_pvalue = exact_pssm_tail_probability(counts, "AA")
        require(abs(cisql_pvalue - exact_pvalue) < 1e-12,
                "PWM p-value differs from exhaustive null enumeration")
        print("[ok] PWM p-value against exhaustive null distribution")

        (workspace / "intervals.gff3").write_text(
            "##gff-version 3\n"
            "chr1\tref\tgene\t1\t10\t.\t+\t.\tID=gene\n"
            "chr1\tref\tgene\t11\t15\t.\t+\t.\tID=touching_gene\n"
            "chr1\tref\texon\t4\t6\t.\t+\t.\tID=exon\n"
            "chr1\tref\tenhancer\t4\t6\t.\t+\t.\tID=enhancer_a\n"
            "chr1\tref\tenhancer\t9\t10\t.\t+\t.\tID=enhancer_b\n",
            encoding="utf-8",
        )
        interval_data = run_query(
            workspace,
            "interval_reference",
            'LOAD ANNOTATION "intervals.gff3" AS annot;\n'
            'EXTRACT GENE AS genes WHERE ID = "gene";\n'
            "EXTRACT GENE AS all_genes;\n"
            "EXTRACT EXON AS exons;\n"
            "EXTRACT ENHANCER AS enhancers;\n"
            "EXCEPT genes FROM exons AS difference;\n"
            "OVERLAPS all_genes WITH enhancers AS supported;\n"
            "NEAR all_genes TO enhancers WITHIN 0 BP AS nearest;\n"
            "COUNT enhancers IN all_genes AS enhancer_counts;\n"
            'EXPORT difference TO "difference.bed" FORMAT BED;\n'
            'EXPORT supported TO "supported.bed" FORMAT BED;\n',
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
        observed_supported = [
            (region["chr"], region["start"], region["end"])
            for region in interval_data["resultSets"]["supported"]
        ]
        expected_supported = select_overlapping_intervals(
            [("chr1", 0, 10), ("chr1", 10, 15)],
            [("chr1", 3, 6), ("chr1", 8, 10)],
        )
        require(observed_supported == expected_supported,
                "OVERLAPS differs from independent interval semi-join")
        print("[ok] directional overlap selection")
        expected_nearest = select_nearest_intervals(
            [("chr1", 0, 10), ("chr1", 10, 15)],
            [("chr1", 3, 6, "enhancer_a"),
             ("chr1", 8, 10, "enhancer_b")],
            0,
        )
        observed_nearest = [
            (
                (region["chr"], region["start"], region["end"]),
                (
                    relation["reference"]["chr"],
                    relation["reference"]["start"],
                    relation["reference"]["end"],
                    relation["reference"]["name"],
                ),
                relation["distance"],
                relation["overlaps"],
            )
            for region in interval_data["resultSets"]["nearest"]
            for relation in [region["spatialRelation"]]
        ]
        require(observed_nearest == expected_nearest,
                "NEAR differs from independent nearest-interval selection")
        print("[ok] bounded nearest-reference selection")
        observed_counts = [
            region["countEvidence"]["count"]
            for region in interval_data["resultSets"]["enhancer_counts"]
        ]
        expected_counts = count_overlaps_by_container(
            [("chr1", 3, 6), ("chr1", 8, 10)],
            [("chr1", 0, 10), ("chr1", 10, 15)],
        )
        require(observed_counts == expected_counts == [2, 0],
                "COUNT differs from independent overlap aggregation")
        print("[ok] overlap counting by container")

        (workspace / "replicate_one.bed").write_text(
            "chr1\t0\t10\tanchor_one\t100\t.\n"
            "chr1\t20\t30\tanchor_two\t100\t.\n",
            encoding="utf-8",
        )
        (workspace / "replicate_two.bed").write_text(
            "chr1\t1\t4\tr2_first_a\t100\t.\n"
            "chr1\t5\t9\tr2_first_b\t100\t.\n"
            "chr1\t21\t25\tr2_second\t100\t.\n",
            encoding="utf-8",
        )
        (workspace / "replicate_three.bed").write_text(
            "chr1\t6\t12\tr3_first\t100\t.\n",
            encoding="utf-8",
        )
        consensus_data = run_query(
            workspace,
            "consensus_reference",
            'LOAD TRACK "replicate_one.bed" FORMAT BED EVIDENCE BINDING '
            'AS replicate_one;\n'
            'LOAD TRACK "replicate_two.bed" FORMAT BED EVIDENCE BINDING '
            'AS replicate_two;\n'
            'LOAD TRACK "replicate_three.bed" FORMAT BED EVIDENCE BINDING '
            'AS replicate_three;\n'
            'CONSENSUS FROM [replicate_one, replicate_two, replicate_three] '
            'ANCHOR replicate_one MIN_SUPPORT 3 AS strict_consensus;\n',
        )
        observed_consensus = consensus_data["resultSets"]["strict_consensus"]
        expected_consensus = anchor_consensus(
            [("chr1", 0, 10), ("chr1", 20, 30)],
            [
                [("chr1", 1, 4), ("chr1", 5, 9), ("chr1", 21, 25)],
                [("chr1", 6, 12)],
            ],
            3,
        )
        require(
            [
                ((region["chr"], region["start"], region["end"]),
                 region["consensusEvidence"]["observedSupport"])
                for region in observed_consensus
            ] == expected_consensus == [(('chr1', 0, 10), 3)] and
            [item["referenceSet"] for item in
             observed_consensus[0]["overlapEvidence"]]
            == ["replicate_two", "replicate_two", "replicate_three"],
            "CONSENSUS differs from independent distinct-set support count",
        )
        print("[ok] anchor-preserving distinct-set consensus")

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
        print(f"[optional] {validate_optional_fimo(workspace, cisql_pvalue)}")

    print("All independent reference validations passed.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AssertionError as error:
        print(f"FAILED: {error}", file=os.sys.stderr)
        raise SystemExit(1)

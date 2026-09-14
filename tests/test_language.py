#!/usr/bin/env python3

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
        (workspace / "data_examples").symlink_to(
            ROOT / "data_examples", target_is_directory=True
        )
        (workspace / "matrices").symlink_to(
            ROOT / "matrices", target_is_directory=True
        )
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
        (workspace / "alternate.fasta").write_text(
            ">chr1\n" + ("T" * len(sequence)) + "\n", encoding="utf-8"
        )
        (workspace / "statistics.fasta").write_text(
            ">chrStats\nAAAA\n", encoding="utf-8"
        )
        (workspace / "mixed_statistics.fasta").write_text(
            ">chrMixed\nAACA\n", encoding="utf-8"
        )
        (workspace / "alternate.gff3").write_text(
            "##gff-version 3\n"
            "chr1\ttest\tgene\t101\t200\t.\t-\t.\tID=alternate\n",
            encoding="utf-8",
        )

        data, _ = run_query(
            workspace,
            "explicit_promoters",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            "DEFINE PROMOTERS OF GENE FROM TSS "
            "UPSTREAM 100 BP DOWNSTREAM 20 BP AS promoters;\n",
        )
        promoters = data["resultSets"]["promoters"]
        require(
            [(item["name"], item["start"], item["end"])
             for item in promoters]
            == [
                ("short_promoter", 0, 20),
                ("long_promoter", 1200, 1320),
                ("tiny_promoter", 2900, 3020),
            ],
            "explicit positive-strand promoter coordinates",
        )
        require(all(len(item["sequence"]) == item["end"] - item["start"]
                    for item in promoters),
                "explicit promoters retain their active-genome sequence")

        data, _ = run_query(
            workspace,
            "negative_strand_promoter",
            'LOAD SEQUENCE "alternate.fasta" AS genome;\n'
            'LOAD ANNOTATION "alternate.gff3" AS annot;\n'
            "DEFINE PROMOTERS OF GENE FROM TSS "
            "UPSTREAM 50 BP DOWNSTREAM 10 BP AS promoters;\n",
        )
        negative = data["resultSets"]["promoters"]
        require(len(negative) == 1 and negative[0]["start"] == 190 and
                negative[0]["end"] == 250 and negative[0]["strand"] == "-",
                "explicit promoters are oriented on the negative strand")

        data, _ = run_query(
            workspace,
            "promoters_from_alias",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            'EXTRACT GENE AS selected WHERE ID = "long";\n'
            "DEFINE PROMOTERS OF selected FROM TSS "
            "UPSTREAM 100 BP DOWNSTREAM 20 BP AS promoters;\n",
        )
        require([item["name"] for item in data["resultSets"]["promoters"]]
                == ["long_promoter"],
                "promoters can be built from a filtered result set")

        data, _ = run_query(
            workspace,
            "no_implicit_promoters",
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            "EXTRACT PROMOTER AS annotated_promoters;\n"
            "EXTRACT TSS AS annotated_tss;\n",
        )
        require(data["resultSets"]["annotated_promoters"] == [] and
                data["resultSets"]["annotated_tss"] == [],
                "GFF import does not invent promoter or TSS features")

        semantic_error = run_invalid_query(
            workspace,
            "zero_promoter_window",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            "DEFINE PROMOTERS OF GENE FROM TSS "
            "UPSTREAM 0 BP DOWNSTREAM 0 BP AS promoters;\n",
            3,
        )
        require("cannot have both distances set to zero" in semantic_error,
                "zero-length promoter windows are rejected")

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
            "explicit_similarity_reference",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            'EXTRACT GENE AS reference WHERE ID = "short";\n'
            "EXTRACT GENE AS matches "
            "WHERE SIMILARITY TO reference >= 99 %;\n",
        )
        require(
            [item["name"] for item in data["resultSets"]["reference"]]
            == ["short"]
            and [item["name"] for item in data["resultSets"]["matches"]]
            == ["short"],
            "explicit similarity reference",
        )

        runtime_error = run_invalid_query(
            workspace,
            "ambiguous_similarity_reference",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            "EXTRACT GENE AS references WHERE LENGTH > 500 BP;\n"
            "EXTRACT GENE AS matches "
            "WHERE SIMILARITY TO references > 50 %;\n",
            4,
        )
        require("must contain exactly one region" in runtime_error,
                "ambiguous explicit similarity references are rejected")

        data, _ = run_query(
            workspace,
            "dataset_selection",
            'LOAD SEQUENCE "fixture.fasta" AS primary_genome;\n'
            'LOAD SEQUENCE "alternate.fasta" AS alternate_genome;\n'
            "USE SEQUENCE primary_genome;\n"
            'FIND MOTIF "AAAA" STRAND POSITIVE AS primary_hits;\n'
            "USE SEQUENCE alternate_genome;\n"
            'FIND MOTIF "AAAA" STRAND POSITIVE AS alternate_hits;\n'
            'LOAD ANNOTATION "fixture.gff3" AS primary_annot;\n'
            'LOAD ANNOTATION "alternate.gff3" AS alternate_annot;\n'
            "USE ANNOTATION primary_annot;\n"
            "EXTRACT GENE AS primary_genes;\n"
            "USE ANNOTATION alternate_annot;\n"
            "EXTRACT GENE AS alternate_genes;\n",
        )
        require(len(data["resultSets"]["primary_hits"]) > 0 and
                len(data["resultSets"]["alternate_hits"]) == 0,
                "explicit sequence dataset selection")
        require(len(data["resultSets"]["primary_genes"]) == 3 and
                [item["name"] for item
                 in data["resultSets"]["alternate_genes"]] == ["alternate"],
                "explicit annotation dataset selection")
        require(data["metadata"]["sequenceDataset"] == "alternate_genome" and
                data["metadata"]["annotationDataset"] == "alternate_annot" and
                data["metadata"]["coordinateSystem"]
                == "zero-based-half-open",
                "result provenance metadata")

        data, _ = run_query(
            workspace,
            "exports",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            'EXTRACT GENE AS genes WHERE ID = "short";\n'
            'EXPORT genes TO "genes.bed" FORMAT BED;\n'
            'EXPORT genes TO "genes.gff3" FORMAT GFF3;\n'
            'EXPORT genes TO "genes.tsv" FORMAT TSV;\n'
            "ANALYZE GC_CONTENT WINDOW 1 KB AS profile;\n"
            'EXPORT profile TO "profile.tsv" FORMAT TSV;\n',
        )
        require((workspace / "genes.bed").read_text(encoding="utf-8")
                == "chr1\t0\t1200\tshort\t0\t+\n",
                "BED export uses zero-based half-open coordinates")
        gff_lines = (workspace / "genes.gff3").read_text(
            encoding="utf-8"
        ).splitlines()
        require(gff_lines[0] == "##gff-version 3" and
                gff_lines[1].split("\t")[3:5] == ["1", "1200"],
                "GFF3 export converts starts to one-based coordinates")
        require((workspace / "genes.tsv").read_text(
                    encoding="utf-8"
                ).splitlines()[0]
                == "chromosome\tstart\tend\tstrand\ttype\tname\tlength"
                "\tmatrix_alias\tmatrix_id\tmatrix_name\tmatrix_source"
                "\traw_score\tscore_percent\tp_value\tq_value\ttested_positions"
                "\tp_value_method\tmultiple_testing_method\tscaled_score"
                "\tscore_range\tscore_scale\tscore_offset"
                "\tbackground_mode\tbackground_source"
                "\tbackground_a\tbackground_c\tbackground_g\tbackground_t"
                "\tbackground_estimation_pseudocount\tbackground_observed_bases"
                "\tbackground_strand_policy\tmotif_pseudocount"
                "\tsource_region\tsource_type\tsource_start\tsource_end"
                "\trelative_start",
                "region TSV export")
        profile_lines = (workspace / "profile.tsv").read_text(
            encoding="utf-8"
        ).splitlines()
        require(profile_lines[0] == "chromosome\tstart\tgc_percent" and
                profile_lines[1].startswith("chr1\t0\t"),
                "GC profile TSV export")

        runtime_error = run_invalid_query(
            workspace,
            "unsafe_export_path",
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            "EXTRACT GENE AS genes;\n"
            'EXPORT genes TO "../outside.bed" FORMAT BED;\n',
            4,
        )
        require("must be relative to the query workspace" in runtime_error,
                "export path traversal is rejected")

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
            "scoped_pwm_scan",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            'EXTRACT GENE AS selected WHERE ID = "short";\n'
            "DEFINE PROMOTERS OF selected FROM TSS "
            "UPSTREAM 0 BP DOWNSTREAM 10 BP AS promoters;\n"
            "SCAN matrix IN promoters STRAND POSITIVE "
            "THRESHOLD 100 % AS sites;\n"
            'EXPORT sites TO "sites.bed" FORMAT BED;\n'
            'EXPORT sites TO "sites.gff3" FORMAT GFF3;\n'
            'EXPORT sites TO "sites.tsv" FORMAT TSV;\n',
        )
        sites = data["resultSets"]["sites"]
        require(len(sites) == 9 and
                [site["start"] for site in sites] == list(range(9)),
                "SCAN IN searches only the selected region and remaps hits")
        require(all(site["type"] == "motif_hit" and
                    site["motifEvidence"]["matrixAlias"] == "matrix" and
                    site["motifEvidence"]["matrixId"] == "TEST" and
                    site["motifEvidence"]["matrixName"] == "test" and
                    site["motifEvidence"]["scorePercent"] == 100 and
                    abs(site["motifEvidence"]["statistics"]["pValue"]
                        - 0.0625) < 1e-12 and
                    abs(site["motifEvidence"]["statistics"]["qValue"]
                        - 0.0625) < 1e-12 and
                    site["motifEvidence"]["statistics"]["testedPositions"]
                    == 9 and
                    site["motifEvidence"]["statistics"]["pValueMethod"]
                    == "zero_order_dynamic_programming" and
                    site["motifEvidence"]["statistics"]
                    ["multipleTestingMethod"] == "benjamini_hochberg" and
                    site["motifEvidence"]["statistics"]["scoreRange"]
                    == 1000 and
                    site["motifEvidence"]["background"]["mode"]
                    == "uniform" and
                    site["motifEvidence"]["background"]["source"]
                    == "default" and
                    site["motifEvidence"]["motifPseudocount"] == 0.1 and
                    site["motifEvidence"]["sourceRegion"]["name"]
                    == "short_promoter" and
                    site["motifEvidence"]["sourceRegion"]["relativeStart"]
                    == site["start"]
                    for site in sites),
                "scoped PWM hits preserve scores and source-region evidence")
        bed_rows = (workspace / "sites.bed").read_text(
            encoding="utf-8"
        ).splitlines()
        require(len(bed_rows) == 9 and
                all(row.split("\t")[4] == "1000" for row in bed_rows),
                "motif BED export maps relative score to the BED scale")
        gff_rows = (workspace / "sites.gff3").read_text(
            encoding="utf-8"
        ).splitlines()
        require("MatrixID=TEST" in gff_rows[1] and
                "MatrixSource=fixture.pwm" in gff_rows[1] and
                "BackgroundMode=uniform" in gff_rows[1] and
                "BackgroundSource=default" in gff_rows[1] and
                "PValue=0.0625" in gff_rows[1] and
                "QValue=0.0625" in gff_rows[1] and
                "TestedPositions=9" in gff_rows[1] and
                "PValueMethod=zero_order_dynamic_programming" in gff_rows[1] and
                "MultipleTestingMethod=benjamini_hochberg" in gff_rows[1] and
                "MotifPseudocount=0.1" in gff_rows[1] and
                "SourceRegion=short_promoter" in gff_rows[1] and
                "SourceRegionType=promoter" in gff_rows[1] and
                "RelativeStart=0" in gff_rows[1],
                "motif GFF3 export retains evidence")
        tsv_rows = (workspace / "sites.tsv").read_text(
            encoding="utf-8"
        ).splitlines()
        first_tsv_site = tsv_rows[1].split("\t")
        require(len(tsv_rows) == 10 and len(first_tsv_site) == 37 and
                first_tsv_site[7:11]
                == ["matrix", "TEST", "test", "fixture.pwm"] and
                float(first_tsv_site[11]) > 0 and
                first_tsv_site[12] == "100" and
                abs(float(first_tsv_site[13]) - 0.0625) < 1e-12 and
                abs(float(first_tsv_site[14]) - 0.0625) < 1e-12 and
                first_tsv_site[15:18]
                == ["9", "zero_order_dynamic_programming",
                    "benjamini_hochberg"] and
                first_tsv_site[19] == "1000" and
                first_tsv_site[22:24] == ["uniform", "default"] and
                first_tsv_site[30] == "forward" and
                abs(float(first_tsv_site[31]) - 0.1) < 1e-12 and
                first_tsv_site[32:34] == ["short_promoter", "promoter"] and
                first_tsv_site[36] == "0",
                "motif TSV export retains evidence columns")

        data, _ = run_query(
            workspace,
            "directional_overlap_semijoin",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            'EXTRACT GENE AS selected WHERE ID = "short";\n'
            "DEFINE PROMOTERS OF selected FROM TSS "
            "UPSTREAM 0 BP DOWNSTREAM 10 BP AS promoters;\n"
            "SCAN matrix STRAND POSITIVE THRESHOLD 100 % AS sites;\n"
            "OVERLAPS sites WITH promoters AS promoter_sites;\n"
            "OVERLAPS promoters WITH sites AS supported_promoters;\n",
        )
        promoter_sites = data["resultSets"]["promoter_sites"]
        require(len(promoter_sites) == 10 and
                [site["start"] for site in promoter_sites]
                == list(range(10)) and
                promoter_sites[-1]["end"] == 11 and
                all(site["motifEvidence"]["matrixId"] == "TEST"
                    for site in promoter_sites),
                "OVERLAPS retains complete query intervals and PWM evidence")
        supported_promoters = data["resultSets"]["supported_promoters"]
        require(len(supported_promoters) == 1 and
                supported_promoters[0]["start"] == 0 and
                supported_promoters[0]["end"] == 10,
                "OVERLAPS emits a supported query once despite multiple hits")

        parser_error = run_invalid_query(
            workspace,
            "overlap_requires_with",
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            "OVERLAPS GENE AND CDS AS invalid;\n",
            2,
        )
        require("'WITH' (for OVERLAPS)" in parser_error,
                "OVERLAPS requires an explicit directional separator")

        showcase_source = (
            ROOT / "cql_examples" / "14_integrated_query.cql"
        ).read_text(encoding="utf-8")
        data, _ = run_query(
            workspace, "anthocyanin_regulatory_showcase", showcase_source
        )
        require(len(data["resultSets"]["supported_myb_sites"]) == 3 and
                [(site["start"], site["strand"])
                 for site in data["resultSets"]["promoter_myb_evidence"]]
                == [(150, "-"), (630, "+")] and
                [(site["start"], site["strand"])
                 for site in data["resultSets"]["enhancer_myb_evidence"]]
                == [(400, "+")] and
                all(site["motifEvidence"]["statistics"]["qValue"] <= 0.01
                    for site in
                    data["resultSets"]["promoter_myb_evidence"]),
                "integrated example separates promoter and enhancer MYB evidence")

        data, _ = run_query(
            workspace,
            "pvalue_scientific_threshold",
            'LOAD SEQUENCE "statistics.fasta" AS genome;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            "SCAN matrix STRAND POSITIVE PVALUE <= 1e-1 AS sites;\n",
        )
        pvalue_sites = data["resultSets"]["sites"]
        require(len(pvalue_sites) == 3 and
                all(abs(site["motifEvidence"]["statistics"]["pValue"]
                        - 0.0625) < 1e-12 and
                    abs(site["motifEvidence"]["statistics"]["qValue"]
                        - 0.0625) < 1e-12
                    for site in pvalue_sites),
                "PVALUE accepts scientific notation without implicit score filtering")

        data, _ = run_query(
            workspace,
            "statistical_threshold_replaces_default_score",
            'LOAD SEQUENCE "mixed_statistics.fasta" AS genome;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            "SCAN matrix STRAND POSITIVE PVALUE <= 1 AS sites;\n",
        )
        require(len(data["resultSets"]["sites"]) == 3 and
                any(site["motifEvidence"]["scorePercent"] < 75
                    for site in data["resultSets"]["sites"]),
                "statistical scans do not inherit the legacy 75% threshold")

        data, _ = run_query(
            workspace,
            "combined_score_and_pvalue_thresholds",
            'LOAD SEQUENCE "mixed_statistics.fasta" AS genome;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            "SCAN matrix STRAND POSITIVE THRESHOLD 100 % "
            "PVALUE <= 1 AS sites;\n",
        )
        require(len(data["resultSets"]["sites"]) == 1 and
                data["resultSets"]["sites"][0]["start"] == 0,
                "explicit score and statistical thresholds compose with AND")

        data, _ = run_query(
            workspace,
            "strict_pvalue_threshold",
            'LOAD SEQUENCE "statistics.fasta" AS genome;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            "SCAN matrix STRAND POSITIVE PVALUE < 6.25e-2 AS sites;\n",
        )
        require(data["resultSets"]["sites"] == [],
                "strict PVALUE thresholds exclude tied boundary values")

        data, _ = run_query(
            workspace,
            "qvalue_test_universe",
            'LOAD SEQUENCE "statistics.fasta" AS genome;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            "SCAN matrix QVALUE <= 1e-1 AS sites;\n",
        )
        require(data["resultSets"]["sites"] == [],
                "QVALUE filtering includes tests from both selected strands")

        data, _ = run_query(
            workspace,
            "inclusive_qvalue_threshold",
            'LOAD SEQUENCE "statistics.fasta" AS genome;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            "SCAN matrix QVALUE <= 0.125 AS sites;\n",
        )
        qvalue_sites = data["resultSets"]["sites"]
        require(len(qvalue_sites) == 3 and
                all(abs(site["motifEvidence"]["statistics"]["qValue"]
                        - 0.125) < 1e-12 and
                    site["motifEvidence"]["statistics"]["testedPositions"]
                    == 6 for site in qvalue_sites),
                "inclusive QVALUE thresholds retain boundary hits")

        semantic_error = run_invalid_query(
            workspace,
            "out_of_range_pvalue",
            'LOAD SEQUENCE "statistics.fasta" AS genome;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            "SCAN matrix PVALUE <= 1.1 AS sites;\n",
            3,
        )
        require("PVALUE threshold must be between 0 and 1"
                in semantic_error,
                "statistical thresholds are probability-bounded")

        parser_error = run_invalid_query(
            workspace,
            "duplicate_statistical_threshold",
            'LOAD SEQUENCE "statistics.fasta" AS genome;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            "SCAN matrix PVALUE <= 0.1 QVALUE <= 0.2 AS sites;\n",
            2,
        )
        require("only one PVALUE or QVALUE" in parser_error,
                "SCAN rejects ambiguous statistical threshold combinations")

        parser_error = run_invalid_query(
            workspace,
            "invalid_statistical_operator",
            'LOAD SEQUENCE "statistics.fasta" AS genome;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            "SCAN matrix PVALUE >= 0.1 AS sites;\n",
            2,
        )
        require("Expected '<' or '<=' after PVALUE" in parser_error,
                "statistical scan options reject reverse comparisons")

        lexical_error = run_invalid_query(
            workspace,
            "invalid_scientific_notation",
            'LOAD SEQUENCE "statistics.fasta" AS genome;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            "SCAN matrix PVALUE <= 1e AS sites;\n",
            2,
        )
        require("Invalid numeric exponent" in lexical_error,
                "malformed scientific notation is a lexical error")

        semantic_error = run_invalid_query(
            workspace,
            "overflowing_statistical_threshold",
            'LOAD SEQUENCE "statistics.fasta" AS genome;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            "SCAN matrix PVALUE <= 1e9999 AS sites;\n",
            3,
        )
        require("PVALUE threshold must be between 0 and 1"
                in semantic_error,
                "overflowing scientific notation is rejected")

        data, _ = run_query(
            workspace,
            "estimated_sequence_background",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            "SCAN matrix STRAND POSITIVE BACKGROUND FROM genome "
            "THRESHOLD 100 % AS sites;\n",
        )
        estimated = data["resultSets"]["sites"][0]["motifEvidence"]
        background = estimated["background"]
        require(background["mode"] == "estimated_zero_order" and
                background["source"] == "genome" and
                background["observedBases"] == len(sequence) and
                background["strandPolicy"] == "forward" and
                background["T"] > 0 and
                abs(sum(background[base] for base in "ACGT") - 1.0) < 1e-6,
                "BACKGROUND FROM sequence estimates and records composition")

        data, _ = run_query(
            workspace,
            "symmetric_sequence_background",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            "SCAN matrix BACKGROUND FROM genome "
            "THRESHOLD 100 % AS sites;\n",
        )
        background = data["resultSets"]["sites"][0]["motifEvidence"][
            "background"
        ]
        require(background["strandPolicy"] == "symmetric" and
                abs(background["A"] - background["T"]) < 1e-12 and
                abs(background["C"] - background["G"]) < 1e-12,
                "two-strand scans symmetrize background frequencies")

        data, _ = run_query(
            workspace,
            "explicit_uniform_background",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            "SCAN matrix STRAND POSITIVE BACKGROUND UNIFORM "
            "THRESHOLD 100 % AS sites;\n",
        )
        background = data["resultSets"]["sites"][0]["motifEvidence"][
            "background"
        ]
        require(background["mode"] == "uniform" and
                background["source"] == "explicit" and
                all(background[base] == 0.25 for base in "ACGT"),
                "explicit uniform backgrounds are distinguished from default")

        data, _ = run_query(
            workspace,
            "estimated_region_background",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            'EXTRACT GENE AS selected WHERE ID = "short";\n'
            "SCAN matrix STRAND POSITIVE BACKGROUND FROM selected "
            "THRESHOLD 100 % AS sites;\n",
        )
        background = data["resultSets"]["sites"][0]["motifEvidence"][
            "background"
        ]
        require(background["source"] == "selected" and
                background["observedBases"] == 1200 and
                background["A"] > 0.8,
                "BACKGROUND FROM region sets uses active-genome intervals")

        data, _ = run_query(
            workspace,
            "negative_source_relative_coordinates",
            'LOAD SEQUENCE "alternate.fasta" AS genome;\n'
            'LOAD ANNOTATION "alternate.gff3" AS annot;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            "DEFINE PROMOTERS OF GENE FROM TSS "
            "UPSTREAM 50 BP DOWNSTREAM 10 BP AS promoters;\n"
            "SCAN matrix IN promoters STRAND NEGATIVE "
            "THRESHOLD 100 % AS sites;\n",
        )
        negative_sites = data["resultSets"]["sites"]
        require(len(negative_sites) == 59 and
                negative_sites[0]["start"] == 190 and
                negative_sites[0]["motifEvidence"]["background"]
                ["strandPolicy"] == "reverse_complement" and
                negative_sites[0]["motifEvidence"]["sourceRegion"]
                ["relativeStart"] == 58 and
                negative_sites[-1]["start"] == 248 and
                negative_sites[-1]["motifEvidence"]["sourceRegion"]
                ["relativeStart"] == 0,
                "source-relative coordinates follow negative-strand orientation")

        semantic_error = run_invalid_query(
            workspace,
            "wrong_scan_target_type",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            "SCAN matrix IN genome THRESHOLD 80 % AS sites;\n",
            3,
        )
        require("SCAN IN expects a region or motif-hit set" in semantic_error,
                "SCAN IN rejects dataset aliases as region targets")

        semantic_error = run_invalid_query(
            workspace,
            "wrong_background_source_type",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            "SCAN matrix BACKGROUND FROM annot "
            "THRESHOLD 80 % AS sites;\n",
            3,
        )
        require("BACKGROUND FROM expects a sequence dataset or region set"
                in semantic_error,
                "BACKGROUND FROM rejects annotation dataset aliases")

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

        semantic_error = run_invalid_query(
            workspace,
            "wrong_dataset_type",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            "USE ANNOTATION genome;\n",
            3,
        )
        require("expects a ANNOTATION_DATA alias" in semantic_error,
                "dataset selector type checking")

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

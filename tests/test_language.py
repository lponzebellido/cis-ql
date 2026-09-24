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
        (workspace / "module.fasta").write_text(
            ">chrModule\nAATTAA\n", encoding="utf-8"
        )
        (workspace / "alternate.gff3").write_text(
            "##gff-version 3\n"
            "chr1\ttest\tgene\t101\t200\t.\t-\t.\tID=alternate\n",
            encoding="utf-8",
        )
        (workspace / "accessibility.narrowPeak").write_text(
            "chr1\t0\t50\topen_promoter\t500\t.\t12.5\t4.2\t3.8\t25\n"
            "chr1\t1250\t1350\topen_gene_edge\t200\t+\t7\t-1\t-1\t-1\n",
            encoding="utf-8",
        )
        (workspace / "binding.narrowPeak").write_text(
            "chr1\t10\t30\tbound_promoter\t800\t.\t18\t7\t5\t5\n"
            "chr1\t1260\t1280\tbound_gene_edge\t600\t.\t11\t5\t3\t10\n",
            encoding="utf-8",
        )
        (workspace / "candidate_regions.bed").write_text(
            "chr1\t10\t20\tcandidate_a\t100\t-\n",
            encoding="utf-8",
        )
        (workspace / "wrong_assembly.bed").write_text(
            "chrMissing\t0\t10\tmissing_chromosome\t100\t.\n",
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
            "load_regulatory_track",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            'LOAD TRACK "accessibility.narrowPeak" FORMAT NARROWPEAK '
            'EVIDENCE ACCESSIBILITY SAMPLE "leaf" CONTROL "input" '
            'ASSAY "ATAC-seq" REPLICATE "R1" CONDITION "pigmented" '
            'AS accessible;\n'
            'LOAD TRACK "binding.narrowPeak" FORMAT NARROWPEAK '
            'EVIDENCE BINDING ASSAY "DAP-seq" SAMPLE "leaf" '
            'CONDITION "pigmented" REPLICATE "R1" CONTROL "input" '
            'AS bound;\n'
            'SCAN matrix IN accessible STRAND POSITIVE THRESHOLD 100 % '
            'AS accessible_sites;\n'
            'EXTRACT accessible AS strong_accessibility '
            'WHERE TRACK_SCORE >= 500 AND SIGNAL_VALUE >= 10;\n'
            'EXTRACT accessible AS statistically_annotated_accessibility '
            'WHERE MINUS_LOG10_PVALUE >= 4 '
            'AND MINUS_LOG10_QVALUE >= 3;\n'
            'EXTRACT accessible AS leaf_accessibility '
            'WHERE EVIDENCE_CLASS = "ACCESSIBILITY" '
            'AND ASSAY = "ATAC-seq" AND SAMPLE = "leaf" '
            'AND CONDITION = "pigmented" AND REPLICATE = "R1" '
            'AND CONTROL = "input";\n'
            'EXTRACT GENE AS genes_with_track_score '
            'WHERE TRACK_SCORE >= 0;\n'
            'OVERLAPS accessible WITH bound AS accessible_and_bound;\n'
            'OVERLAPS accessible WITH GENE AS gene_accessible_peaks;\n'
            'INTERSECT accessible AND GENE AS clipped_accessible;\n'
            'EXPORT accessible TO "accessible.tsv" FORMAT TSV;\n'
            'EXPORT accessible TO "accessible.gff3" FORMAT GFF3;\n'
            'EXPORT accessible_and_bound TO "accessible_bound.tsv" '
            'FORMAT TSV;\n'
            'EXPORT accessible_and_bound TO "accessible_bound.gff3" '
            'FORMAT GFF3;\n',
        )
        accessible = data["resultSets"]["accessible"]
        accessible_sites = data["resultSets"]["accessible_sites"]
        accessible_and_bound = data["resultSets"]["accessible_and_bound"]
        clipped_accessible = data["resultSets"]["clipped_accessible"]
        require(len(accessible) == 2 and
                len(data["resultSets"]["gene_accessible_peaks"]) == 2 and
                len(accessible_sites) == 49 and
                all(site["trackEvidence"]["trackAlias"] == "accessible"
                    for site in accessible_sites) and
                all(site["trackEvidence"].get("peakPosition") == 25
                    for site in accessible_sites) and
                "trackEvidence" in clipped_accessible[0] and
                "trackEvidence" not in clipped_accessible[1] and
                accessible[0]["sequence"] == "A" * 50 and
                accessible[0]["trackEvidence"] == {
                    "trackAlias": "accessible",
                    "source": "accessibility.narrowPeak",
                    "format": "NARROWPEAK",
                    "evidenceClass": "ACCESSIBILITY",
                    "assay": "ATAC-seq",
                    "sample": "leaf",
                    "condition": "pigmented",
                    "replicate": "R1",
                    "control": "input",
                    "score": 500,
                    "signalValue": 12.5,
                    "minusLog10PValue": 4.2,
                    "minusLog10QValue": 3.8,
                    "peakOffset": 25,
                    "peakPosition": 25,
                } and
                accessible[1]["trackEvidence"] == {
                    "trackAlias": "accessible",
                    "source": "accessibility.narrowPeak",
                    "format": "NARROWPEAK",
                    "evidenceClass": "ACCESSIBILITY",
                    "assay": "ATAC-seq",
                    "sample": "leaf",
                    "condition": "pigmented",
                    "replicate": "R1",
                    "control": "input",
                    "score": 200,
                    "signalValue": 7,
                },
                "LOAD TRACK preserves narrowPeak evidence and is queryable")
        require(
            [region["name"] for region in
             data["resultSets"]["strong_accessibility"]]
            == ["open_promoter"] and
            [region["name"] for region in data["resultSets"]
             ["statistically_annotated_accessibility"]]
            == ["open_promoter"] and
            len(data["resultSets"]["leaf_accessibility"]) == 2 and
            data["resultSets"]["genes_with_track_score"] == [],
            "track evidence properties are executable WHERE filters",
        )
        require(
            len(accessible_and_bound) == 2 and
            all(region["trackEvidence"]["evidenceClass"]
                == "ACCESSIBILITY" for region in accessible_and_bound) and
            [item["reference"]["name"] for item in
             accessible_and_bound[0]["overlapEvidence"]]
            == ["bound_promoter"] and
            accessible_and_bound[0]["overlapEvidence"][0]
                ["referenceSet"] == "bound" and
            accessible_and_bound[0]["overlapEvidence"][0]
                ["trackEvidence"]["evidenceClass"] == "BINDING" and
            accessible_and_bound[0]["overlapEvidence"][0]
                ["trackEvidence"]["peakPosition"] == 15,
            "OVERLAPS retains query and reference regulatory evidence",
        )
        track_tsv = (workspace / "accessible.tsv").read_text(
            encoding="utf-8"
        ).splitlines()[1].split("\t")
        require(len(track_tsv) == 103 and track_tsv[81:96] == [
                    "accessible", "accessibility.narrowPeak", "NARROWPEAK",
                    "ACCESSIBILITY", "ATAC-seq", "leaf", "pigmented",
                    "R1", "input", "500", "12.5",
                    "4.2000000000000002",
                    "3.7999999999999998", "25", "25",
                ], "TSV export preserves typed narrowPeak evidence")
        combined_tsv = (workspace / "accessible_bound.tsv").read_text(
            encoding="utf-8"
        ).splitlines()[1].split("\t")
        combined_overlap = json.loads(combined_tsv[96])
        require(len(combined_tsv) == 103 and
                combined_overlap[0]["referenceSet"] == "bound" and
                combined_overlap[0]["trackEvidence"]["evidenceClass"]
                == "BINDING" and
                combined_overlap[0]["trackEvidence"]["peakPosition"] == 15,
                "TSV export preserves overlap evidence as structured JSON")
        track_gff = (workspace / "accessible.gff3").read_text(
            encoding="utf-8"
        ).splitlines()[1]
        require("TrackAlias=accessible" in track_gff and
                "TrackFormat=NARROWPEAK" in track_gff and
                "EvidenceClass=ACCESSIBILITY" in track_gff and
                "Assay=ATAC-seq" in track_gff and
                "Sample=leaf" in track_gff and
                "Condition=pigmented" in track_gff and
                "Replicate=R1" in track_gff and
                "Control=input" in track_gff and
                "SignalValue=12.5" in track_gff and
                "TrackMinusLog10PValue=4.2" in track_gff and
                "PeakPosition=25" in track_gff,
                "GFF3 export preserves narrowPeak evidence")
        combined_gff = (workspace / "accessible_bound.gff3").read_text(
            encoding="utf-8"
        ).splitlines()[1]
        require("OverlapEvidenceCount=1" in combined_gff and
                "OverlapEvidenceJSON=" in combined_gff and
                "BINDING" in combined_gff,
                "GFF3 export preserves overlapping track evidence")

        data, _ = run_query(
            workspace,
            "load_bed_track",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD TRACK "candidate_regions.bed" FORMAT BED '
            'EVIDENCE OTHER AS candidates;\n',
        )
        bed_region = data["resultSets"]["candidates"][0]
        require((bed_region["start"], bed_region["end"],
                 bed_region["strand"], bed_region["sequence"])
                == (10, 20, "-", "A" * 10) and
                bed_region["trackEvidence"] == {
                    "trackAlias": "candidates",
                    "source": "candidate_regions.bed",
                    "format": "BED",
                    "evidenceClass": "OTHER",
                    "score": 100,
                }, "LOAD TRACK supports BED6 and attaches genome sequence")

        parser_error = run_invalid_query(
            workspace,
            "track_requires_format",
            'LOAD TRACK "candidate_regions.bed" AS candidates;\n',
            2,
        )
        require("Expected 'FORMAT'" in parser_error,
                "LOAD TRACK requires an explicit input format")

        parser_error = run_invalid_query(
            workspace,
            "track_requires_evidence_class",
            'LOAD TRACK "candidate_regions.bed" FORMAT BED AS candidates;\n',
            2,
        )
        require("Expected 'EVIDENCE'" in parser_error,
                "LOAD TRACK requires an explicit evidence class")

        parser_error = run_invalid_query(
            workspace,
            "track_rejects_duplicate_metadata",
            'LOAD TRACK "candidate_regions.bed" FORMAT BED '
            'EVIDENCE OTHER REPLICATE "R1" REPLICATE "R2" '
            'AS candidates;\n',
            2,
        )
        require("Duplicate REPLICATE clause" in parser_error,
                "LOAD TRACK rejects contradictory duplicate metadata")

        semantic_error = run_invalid_query(
            workspace,
            "consensus_requires_listed_anchor",
            'LOAD TRACK "accessibility.narrowPeak" FORMAT NARROWPEAK '
            'EVIDENCE ACCESSIBILITY AS first;\n'
            'LOAD TRACK "binding.narrowPeak" FORMAT NARROWPEAK '
            'EVIDENCE BINDING AS second;\n'
            'CONSENSUS FROM [first, second] ANCHOR missing '
            'MIN_SUPPORT 2 AS invalid;\n',
            3,
        )
        require("must also appear in the input list" in semantic_error,
                "CONSENSUS requires an explicit listed anchor")

        semantic_error = run_invalid_query(
            workspace,
            "consensus_checks_support_threshold",
            'LOAD TRACK "accessibility.narrowPeak" FORMAT NARROWPEAK '
            'EVIDENCE ACCESSIBILITY AS first;\n'
            'LOAD TRACK "binding.narrowPeak" FORMAT NARROWPEAK '
            'EVIDENCE BINDING AS second;\n'
            'CONSENSUS FROM [first, second] ANCHOR first '
            'MIN_SUPPORT 3 AS invalid;\n',
            3,
        )
        require("MIN_SUPPORT must be between 2" in semantic_error,
                "CONSENSUS bounds support by its distinct input sets")

        semantic_error = run_invalid_query(
            workspace,
            "consensus_checks_reciprocal_overlap",
            'LOAD TRACK "accessibility.narrowPeak" FORMAT NARROWPEAK '
            'EVIDENCE ACCESSIBILITY AS first;\n'
            'LOAD TRACK "binding.narrowPeak" FORMAT NARROWPEAK '
            'EVIDENCE BINDING AS second;\n'
            'CONSENSUS FROM [first, second] ANCHOR first '
            'MIN_SUPPORT 2 MIN_RECIPROCAL_OVERLAP 101 % AS invalid;\n',
            3,
        )
        require("must be greater than 0% and at most 100%" in semantic_error,
                "CONSENSUS bounds reciprocal-overlap percentages")

        semantic_error = run_invalid_query(
            workspace,
            "consensus_checks_summit_distance",
            'LOAD TRACK "accessibility.narrowPeak" FORMAT NARROWPEAK '
            'EVIDENCE ACCESSIBILITY AS first;\n'
            'LOAD TRACK "binding.narrowPeak" FORMAT NARROWPEAK '
            'EVIDENCE BINDING AS second;\n'
            'CONSENSUS FROM [first, second] ANCHOR first '
            'MIN_SUPPORT 2 MAX_SUMMIT_DISTANCE 0.5 BP AS invalid;\n',
            3,
        )
        require("must resolve to a whole number of base pairs"
                in semantic_error,
                "CONSENSUS requires integral summit distances")

        runtime_error = run_invalid_query(
            workspace,
            "track_checks_active_assembly",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD TRACK "wrong_assembly.bed" FORMAT BED '
            'EVIDENCE OTHER AS candidates;\n',
            4,
        )
        require("absent from the active sequence dataset" in runtime_error,
                "LOAD TRACK rejects chromosome mismatches against active FASTA")

        semantic_error = run_invalid_query(
            workspace,
            "track_score_requires_number",
            'LOAD TRACK "candidate_regions.bed" FORMAT BED '
            'EVIDENCE OTHER AS candidates;\n'
            'EXTRACT candidates AS invalid WHERE TRACK_SCORE >= "high";\n',
            3,
        )
        require("TRACK_SCORE must be compared with a finite" in semantic_error,
                "track scores reject string thresholds")

        semantic_error = run_invalid_query(
            workspace,
            "track_signal_rejects_units",
            'LOAD TRACK "candidate_regions.bed" FORMAT BED '
            'EVIDENCE OTHER AS candidates;\n'
            'EXTRACT candidates AS invalid WHERE SIGNAL_VALUE >= 1 BP;\n',
            3,
        )
        require("SIGNAL_VALUE must be compared with a finite"
                in semantic_error,
                "track evidence thresholds reject genomic units")

        semantic_error = run_invalid_query(
            workspace,
            "track_metadata_requires_equality",
            'LOAD TRACK "candidate_regions.bed" FORMAT BED '
            'EVIDENCE OTHER AS candidates;\n'
            'EXTRACT candidates AS invalid '
            'WHERE EVIDENCE_CLASS >= "OTHER";\n',
            3,
        )
        require("EVIDENCE_CLASS supports only equality" in semantic_error,
                "track metadata uses exact string equality")

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
                "\trelative_start"
                "\tspatial_relation\treference_set\treference_chr"
                "\treference_start\treference_end\treference_strand"
                "\treference_type\treference_name\tdistance_bp"
                "\tmaximum_distance_bp\toverlaps"
                "\tcount_relation\tcounted_set\tcontainer_set\toverlap_count"
                "\tmodule_minimum_spacing_bp\tmodule_maximum_spacing_bp"
                "\tmodule_observed_spacing_bp\tmodule_order_policy"
                "\tmodule_observed_order\tmodule_orientation_policy"
                "\tmodule_observed_orientation"
                "\tfirst_set\tfirst_chr\tfirst_start\tfirst_end"
                "\tfirst_strand\tfirst_type\tfirst_name\tfirst_matrix_id"
                "\tfirst_raw_score\tfirst_p_value\tfirst_q_value"
                "\tsecond_set\tsecond_chr\tsecond_start\tsecond_end"
                "\tsecond_strand\tsecond_type\tsecond_name"
                "\tsecond_matrix_id\tsecond_raw_score\tsecond_p_value"
                "\tsecond_q_value"
                "\ttrack_alias\ttrack_source\ttrack_format\tevidence_class"
                "\tassay\tsample\tcondition\treplicate\tcontrol"
                "\ttrack_score"
                "\tsignal_value\ttrack_minus_log10_p_value"
                "\ttrack_minus_log10_q_value\tpeak_offset\tpeak_position"
                "\toverlap_evidence_json"
                "\tconsensus_anchor_set\tconsensus_minimum_support"
                "\tconsensus_observed_support"
                "\tconsensus_minimum_reciprocal_overlap_percent"
                "\tconsensus_maximum_summit_distance_bp"
                "\tconsensus_evidence_json",
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
        require(len(tsv_rows) == 10 and len(first_tsv_site) == 103 and
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
                first_tsv_site[36] == "0" and
                first_tsv_site[37:] == [""] * 66,
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
            "OVERLAPS promoters WITH sites AS supported_promoters;\n"
            'EXPORT promoter_sites TO "promoter_sites.tsv" FORMAT TSV;\n',
        )
        promoter_sites = data["resultSets"]["promoter_sites"]
        require(len(promoter_sites) == 10 and
                [site["start"] for site in promoter_sites]
                == list(range(10)) and
                promoter_sites[-1]["end"] == 11 and
                all(site["motifEvidence"]["matrixId"] == "TEST"
                    for site in promoter_sites) and
                all(site["overlapEvidence"][0]["referenceSet"]
                    == "promoters" for site in promoter_sites),
                "OVERLAPS retains complete query intervals and PWM evidence")
        supported_promoters = data["resultSets"]["supported_promoters"]
        require(len(supported_promoters) == 1 and
                supported_promoters[0]["start"] == 0 and
                supported_promoters[0]["end"] == 10 and
                len(supported_promoters[0]["overlapEvidence"]) == 10 and
                all(item["referenceSet"] == "sites" for item in
                    supported_promoters[0]["overlapEvidence"]),
                "OVERLAPS emits a supported query once despite multiple hits")
        promoter_overlap_tsv = (workspace / "promoter_sites.tsv").read_text(
            encoding="utf-8"
        ).splitlines()[1].split("\t")
        require(len(promoter_overlap_tsv) == 103 and
                json.loads(promoter_overlap_tsv[96])[0]["reference"]
                    ["name"] == "short_promoter",
                "OVERLAPS TSV evidence identifies the matching reference")

        parser_error = run_invalid_query(
            workspace,
            "overlap_requires_with",
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            "OVERLAPS GENE AND CDS AS invalid;\n",
            2,
        )
        require("'WITH' (for OVERLAPS)" in parser_error,
                "OVERLAPS requires an explicit directional separator")

        data, _ = run_query(
            workspace,
            "nearest_reference_evidence",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            "DEFINE PROMOTERS OF GENE FROM TSS "
            "UPSTREAM 100 BP DOWNSTREAM 20 BP AS promoters;\n"
            "NEAR promoters TO GENE WITHIN 0 BP AS linked_promoters;\n"
            'EXPORT linked_promoters TO "linked.gff3" FORMAT GFF3;\n'
            'EXPORT linked_promoters TO "linked.tsv" FORMAT TSV;\n',
        )
        linked_promoters = data["resultSets"]["linked_promoters"]
        first_relation = linked_promoters[0]["spatialRelation"]
        require(len(linked_promoters) == 3 and
                first_relation["relation"] == "NEAR" and
                first_relation["referenceSet"] == "GENE" and
                first_relation["reference"] == {
                    "chr": "chr1", "start": 0, "end": 1200,
                    "strand": "+", "type": "gene", "name": "short",
                } and
                first_relation["distance"] == 0 and
                first_relation["maximumDistance"] == 0 and
                first_relation["overlaps"] is True,
                "NEAR emits auditable nearest-reference evidence")
        linked_gff = (workspace / "linked.gff3").read_text(
            encoding="utf-8"
        ).splitlines()[1]
        require("SpatialRelation=NEAR" in linked_gff and
                "ReferenceSet=GENE" in linked_gff and
                "ReferenceName=short" in linked_gff and
                "Distance=0" in linked_gff and
                "MaximumDistance=0" in linked_gff and
                "Overlaps=true" in linked_gff,
                "GFF3 export retains nearest-reference evidence")
        linked_tsv = (workspace / "linked.tsv").read_text(
            encoding="utf-8"
        ).splitlines()[1].split("\t")
        require(len(linked_tsv) == 103 and linked_tsv[37:48] == [
                    "NEAR", "GENE", "chr1", "0", "1200", "+", "gene",
                    "short", "0", "0", "true",
                ] and linked_tsv[48:] == [""] * 55,
                "TSV export retains typed nearest-reference evidence")

        data, _ = run_query(
            workspace,
            "count_overlapping_motif_support",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            "DEFINE PROMOTERS OF GENE FROM TSS "
            "UPSTREAM 100 BP DOWNSTREAM 20 BP AS promoters;\n"
            "SCAN matrix STRAND POSITIVE THRESHOLD 100 % AS sites;\n"
            "COUNT sites IN promoters AS promoter_counts;\n"
            "COUNT sites IN promoters AS supported_promoters "
            "WHERE COUNT >= 1;\n"
            'EXPORT promoter_counts TO "counts.gff3" FORMAT GFF3;\n'
            'EXPORT promoter_counts TO "counts.tsv" FORMAT TSV;\n',
        )
        promoter_counts = data["resultSets"]["promoter_counts"]
        require([region["countEvidence"]["count"]
                 for region in promoter_counts] == [20, 0, 0] and
                promoter_counts[0]["countEvidence"] == {
                    "relation": "OVERLAPS",
                    "countedSet": "sites",
                    "containerSet": "promoters",
                    "count": 20,
                } and
                len(data["resultSets"]["supported_promoters"]) == 1 and
                data["resultSets"]["supported_promoters"][0]["name"]
                == "short_promoter" and
                data["resultSets"]["supported_promoters"][0]
                    ["countEvidence"]["count"] == 20,
                "COUNT preserves zeroes and filters promoters by motif support")
        count_gff = (workspace / "counts.gff3").read_text(
            encoding="utf-8"
        ).splitlines()[1]
        require("CountRelation=OVERLAPS" in count_gff and
                "CountedSet=sites" in count_gff and
                "ContainerSet=promoters" in count_gff and
                "OverlapCount=20" in count_gff,
                "GFF3 export retains count provenance")
        count_tsv = (workspace / "counts.tsv").read_text(
            encoding="utf-8"
        ).splitlines()[1].split("\t")
        require(len(count_tsv) == 103 and count_tsv[48:52] == [
                    "OVERLAPS", "sites", "promoters", "20",
                ] and count_tsv[52:] == [""] * 51,
                "TSV export retains count provenance")

        data, _ = run_query(
            workspace,
            "explicit_motif_module",
            'LOAD SEQUENCE "module.fasta" AS genome;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            "SCAN matrix STRAND POSITIVE THRESHOLD 100 % AS sites;\n"
            "DEFINE MODULE FROM sites WITH sites "
            "SPACING 2 BP TO 2 BP ORDER ANY ORIENTATION SAME "
            "AS homotypic_modules;\n"
            "DEFINE MODULE FROM sites WITH sites "
            "SPACING 2 BP TO 2 BP ORDER AS_WRITTEN ORIENTATION SAME "
            "AS ordered_modules;\n"
            "DEFINE MODULE FROM sites WITH sites "
            "SPACING 2 BP TO 2 BP ORDER ANY ORIENTATION OPPOSITE "
            "AS opposite_modules;\n"
            'EXPORT homotypic_modules TO "modules.gff3" FORMAT GFF3;\n'
            'EXPORT homotypic_modules TO "modules.tsv" FORMAT TSV;\n',
        )
        modules = data["resultSets"]["homotypic_modules"]
        module = modules[0]
        evidence = module["moduleEvidence"]
        require(len(modules) == 1 and
                (module["start"], module["end"], module["sequence"])
                == (0, 6, "AATTAA") and
                evidence["spacing"] == {
                    "minimum": 2, "maximum": 2, "observed": 2,
                } and
                evidence["order"] == {
                    "policy": "ANY", "observed": "FIRST_BEFORE_SECOND",
                } and
                evidence["orientation"] == {
                    "policy": "SAME", "observed": "SAME",
                } and
                [(member["start"], member["end"], member["strand"])
                 for member in evidence["members"]]
                == [(0, 2, "+"), (4, 6, "+")] and
                all(member["sourceSet"] == "sites" and
                    member["motifEvidence"]["matrixId"] == "TEST"
                    for member in evidence["members"]),
                "DEFINE MODULE emits one auditable homotypic pair")
        require(len(data["resultSets"]["ordered_modules"]) == 1 and
                data["resultSets"]["ordered_modules"][0]
                    ["moduleEvidence"]["order"]["policy"] == "AS_WRITTEN" and
                data["resultSets"]["opposite_modules"] == [],
                "module order and orientation policies are executable syntax")
        module_gff = (workspace / "modules.gff3").read_text(
            encoding="utf-8"
        ).splitlines()[1]
        require("ModuleObservedSpacing=2" in module_gff and
                "ModuleOrderPolicy=ANY" in module_gff and
                "ModuleOrientationPolicy=SAME" in module_gff and
                "FirstMatrixID=TEST" in module_gff and
                "SecondMatrixID=TEST" in module_gff,
                "GFF3 export retains module constraints and members")
        module_tsv = (workspace / "modules.tsv").read_text(
            encoding="utf-8"
        ).splitlines()[1].split("\t")
        require(len(module_tsv) == 103 and
                module_tsv[52:59] == [
                    "2", "2", "2", "ANY", "FIRST_BEFORE_SECOND",
                    "SAME", "SAME",
                ] and
                module_tsv[59:67] == [
                    "sites", "chrModule", "0", "2", "+", "motif_hit",
                    "sites_chrModule_0_plus", "TEST",
                ] and
                module_tsv[70:78] == [
                    "sites", "chrModule", "4", "6", "+", "motif_hit",
                    "sites_chrModule_4_plus", "TEST",
                ],
                "TSV export retains typed module evidence columns")

        semantic_error = run_invalid_query(
            workspace,
            "module_rejects_reversed_spacing",
            'LOAD SEQUENCE "module.fasta" AS genome;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            "SCAN matrix STRAND POSITIVE THRESHOLD 100 % AS sites;\n"
            "DEFINE MODULE FROM sites WITH sites "
            "SPACING 3 BP TO 2 BP ORDER ANY ORIENTATION SAME AS invalid;\n",
            3,
        )
        require("minimum spacing cannot exceed" in semantic_error,
                "module spacing ranges reject reversed bounds")

        semantic_error = run_invalid_query(
            workspace,
            "module_rejects_fractional_base_pairs",
            'LOAD SEQUENCE "module.fasta" AS genome;\n'
            'LOAD MATRIX "fixture.pwm" AS matrix;\n'
            "SCAN matrix STRAND POSITIVE THRESHOLD 100 % AS sites;\n"
            "DEFINE MODULE FROM sites WITH sites "
            "SPACING 0.5 BP TO 2 BP ORDER ANY ORIENTATION SAME AS invalid;\n",
            3,
        )
        require("whole numbers of base pairs" in semantic_error,
                "module spacing rejects sub-base-pair bounds")

        semantic_error = run_invalid_query(
            workspace,
            "module_rejects_sequence_alias",
            'LOAD SEQUENCE "module.fasta" AS genome;\n'
            "DEFINE MODULE FROM genome WITH genome "
            "SPACING 0 BP TO 2 BP ORDER ANY ORIENTATION ANY AS invalid;\n",
            3,
        )
        require("requires region or motif-hit aliases" in semantic_error,
                "DEFINE MODULE rejects non-region datasets")

        parser_error = run_invalid_query(
            workspace,
            "count_requires_in",
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            "COUNT GENE WITH PROMOTER AS invalid;\n",
            2,
        )
        require("Expected 'IN'" in parser_error,
                "COUNT requires an explicit container relation")

        semantic_error = run_invalid_query(
            workspace,
            "count_filter_requires_whole_number",
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            "COUNT GENE IN GENE AS invalid WHERE COUNT >= 0.5;\n",
            3,
        )
        require("finite, non-negative whole number" in semantic_error,
                "COUNT filters reject fractional thresholds")

        semantic_error = run_invalid_query(
            workspace,
            "count_rejects_dataset_alias",
            'LOAD SEQUENCE "fixture.fasta" AS genome;\n'
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            "COUNT genome IN GENE AS invalid;\n",
            3,
        )
        require("requires region or motif-hit sets" in semantic_error,
                "COUNT rejects non-region aliases")

        (workspace / "distance_units.gff3").write_text(
            "##gff-version 3\n"
            "chrU\ttest\tgene\t1\t10\t.\t+\t.\tID=query_gene\n"
            "chrU\ttest\tenhancer\t15\t20\t.\t.\t.\tID=reference_enhancer\n",
            encoding="utf-8",
        )
        data, _ = run_query(
            workspace,
            "near_converts_distance_units",
            'LOAD ANNOTATION "distance_units.gff3" AS annot;\n'
            "NEAR GENE TO ENHANCER WITHIN 0.004 KB AS linked;\n",
        )
        converted_relation = data["resultSets"]["linked"][0][
            "spatialRelation"
        ]
        require(converted_relation["distance"] == 4 and
                converted_relation["maximumDistance"] == 4,
                "NEAR converts explicit distance units to whole base pairs")

        parser_error = run_invalid_query(
            workspace,
            "near_requires_to",
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            "NEAR GENE WITH CDS WITHIN 5 BP AS invalid;\n",
            2,
        )
        require("'TO' (for NEAR)" in parser_error,
                "NEAR requires its directional separator")

        parser_error = run_invalid_query(
            workspace,
            "near_requires_unit",
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            "NEAR GENE TO CDS WITHIN 5 AS invalid;\n",
            2,
        )
        require("Expected BP, KB, or MB" in parser_error,
                "NEAR requires an explicit genomic distance unit")

        semantic_error = run_invalid_query(
            workspace,
            "near_requires_whole_base_pairs",
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            "NEAR GENE TO CDS WITHIN 0.0001 KB AS invalid;\n",
            3,
        )
        require("whole number of base pairs" in semantic_error,
                "NEAR rejects sub-base-pair distances")

        semantic_error = run_invalid_query(
            workspace,
            "near_rejects_non_finite_distance",
            'LOAD ANNOTATION "fixture.gff3" AS annot;\n'
            "NEAR GENE TO CDS WITHIN 1e999 MB AS invalid;\n",
            3,
        )
        require("finite, non-negative" in semantic_error,
                "NEAR rejects non-finite or overflowing distances")

        showcase_source = (
            ROOT / "cql_examples" / "08_integrated_anthocyanin_query.cql"
        ).read_text(encoding="utf-8")
        data, _ = run_query(
            workspace, "anthocyanin_regulatory_showcase", showcase_source
        )
        require(len(data["resultSets"]["supported_myb_sites"]) == 4 and
                [(site["start"], site["strand"])
                 for site in data["resultSets"]["promoter_myb_evidence"]]
                == [(150, "-"), (630, "+")] and
                [(site["start"], site["strand"])
                 for site in data["resultSets"]["enhancer_myb_evidence"]]
                == [(400, "+"), (420, "+")] and
                all(site["motifEvidence"]["statistics"]["qValue"] <= 0.01
                    for site in
                    data["resultSets"]["promoter_myb_evidence"]) and
                len(data["resultSets"]["enhancer_myb_modules"]) == 1 and
                data["resultSets"]["enhancer_myb_modules"][0]
                    ["moduleEvidence"]["spacing"]["observed"] == 11 and
                data["resultSets"]["enhancer_myb_modules"][0]
                    ["moduleEvidence"]["orientation"]
                == {"policy": "ANY", "observed": "SAME"} and
                [member["motifEvidence"]["matrixId"] for member in
                 data["resultSets"]["enhancer_myb_modules"][0]
                    ["moduleEvidence"]["members"]]
                == ["MA0054.1", "MA0054.1"] and
                len(data["resultSets"]
                    ["module_nearest_gene_hypotheses"]) == 1 and
                data["resultSets"]["module_nearest_gene_hypotheses"][0]
                    ["spatialRelation"]["reference"]["name"]
                    == "transporter_candidate" and
                data["resultSets"]["module_nearest_gene_hypotheses"][0]
                    ["spatialRelation"]["distance"] == 71 and
                data["resultSets"]["module_nearest_gene_hypotheses"][0]
                    ["spatialRelation"]["maximumDistance"] == 100 and
                data["resultSets"]["module_nearest_gene_hypotheses"][0]
                    ["spatialRelation"]["overlaps"] is False and
                [region["countEvidence"]["count"] for region in
                 data["resultSets"]["promoter_myb_counts"]] == [1, 1, 0] and
                [region["name"] for region in
                 data["resultSets"]["supported_promoter_candidates"]]
                == ["anthocyanin_enzyme_promoter",
                    "transporter_candidate_promoter"],
                "integrated example separates regulatory evidence, defines an auditable module, and records a provisional nearest-gene link")

        accessibility_example = (
            ROOT / "cql_examples" / "09_accessible_myb_evidence.cql"
        ).read_text(encoding="utf-8")
        data, _ = run_query(
            workspace, "anthocyanin_accessibility_evidence",
            accessibility_example
        )
        peak_counts = data["resultSets"]["accessible_peak_myb_counts"]
        candidate_links = data["resultSets"][
            "accessible_myb_candidate_links"
        ]
        require([peak["countEvidence"]["count"] for peak in peak_counts]
                == [1, 2, 1, 0] and
                all(peak.get("trackEvidence", {}).get("evidenceClass")
                    == "ACCESSIBILITY" and
                    peak["trackEvidence"].get("assay")
                    == "synthetic ATAC-seq-like fixture" and
                    peak["trackEvidence"].get("sample")
                    == "synthetic anthocyanin locus" and
                    peak["trackEvidence"].get("condition")
                    == "pigmented petal" and
                    peak["trackEvidence"].get("replicate") == "A1" and
                    peak["trackEvidence"].get("control") == "input"
                    for peak in peak_counts) and
                [peak["name"] for peak in candidate_links] == [
                    "anthocyanin_promoter_accessible",
                    "distal_enhancer_accessible",
                    "transporter_promoter_accessible",
                ] and
                all(link.get("trackEvidence", {}).get("evidenceClass")
                    == "ACCESSIBILITY" for link in candidate_links) and
                [peak["spatialRelation"]["distance"]
                 for peak in candidate_links] == [25, 60, 20],
                "accessibility example combines track, motif-count, and proximity evidence")

        multi_evidence_example = (
            ROOT / "cql_examples" /
            "10_accessible_bound_myb_candidates.cql"
        ).read_text(encoding="utf-8")
        data, _ = run_query(
            workspace, "anthocyanin_multi_track_evidence",
            multi_evidence_example
        )
        combined_peaks = data["resultSets"]["accessible_bound_peaks"]
        combined_counts = data["resultSets"][
            "accessible_bound_peak_myb_counts"
        ]
        combined_links = data["resultSets"][
            "multi_evidence_gene_hypotheses"
        ]
        require(
            len(combined_peaks) == 2 and
            [peak["countEvidence"]["count"] for peak in combined_counts]
            == [1, 2] and
            all(peak["trackEvidence"]["evidenceClass"] == "ACCESSIBILITY"
                and peak["trackEvidence"]["replicate"] == "A1"
                for peak in combined_counts) and
            all(len(peak["overlapEvidence"]) == 1 and
                peak["overlapEvidence"][0]["trackEvidence"]
                    ["evidenceClass"] == "BINDING" and
                peak["overlapEvidence"][0]["trackEvidence"]
                    ["replicate"] == "R1" and
                peak["overlapEvidence"][0]["trackEvidence"]
                    ["control"] == "mock"
                for peak in combined_counts) and
            combined_counts[0]["overlapEvidence"][0]["trackEvidence"]
                ["peakPosition"] == 160 and
            all(link["overlapEvidence"][0]["trackEvidence"]
                    ["evidenceClass"] == "BINDING"
                for link in combined_links) and
            [link["spatialRelation"]["distance"] for link in combined_links]
            == [25, 60],
            "multi-track example retains accessibility, binding, motif, and proximity evidence",
        )

        replicate_example = (
            ROOT / "cql_examples" /
            "11_replicate_supported_candidates.cql"
        ).read_text(encoding="utf-8")
        data, _ = run_query(
            workspace, "anthocyanin_replicate_supported_candidates",
            replicate_example
        )
        replicate_peaks = data["resultSets"][
            "replicate_supported_accessible_peaks"
        ]
        binding_consensus = data["resultSets"][
            "replicate_supported_binding_peaks"
        ]
        replicate_counts = data["resultSets"][
            "replicate_supported_peak_myb_counts"
        ]
        replicate_links = data["resultSets"][
            "replicate_supported_gene_hypotheses"
        ]
        direct_support = replicate_peaks[0]["overlapEvidence"][0]
        nested_support = direct_support["supportingEvidence"][0]
        require(
            len(replicate_peaks) == 2 and
            len(binding_consensus) == 2 and
            all(peak["consensusEvidence"] == {
                    "anchorSet": "strong_myb_binding_rep1",
                    "minimumSupport": 2,
                    "observedSupport": 2,
                    "minimumReciprocalOverlapPercent": 50,
                    "maximumSummitDistanceBp": 5,
                    "inputSets": ["strong_myb_binding_rep1",
                                  "strong_myb_binding_rep2"],
                } for peak in binding_consensus) and
            [peak["countEvidence"]["count"] for peak in replicate_counts]
            == [1, 2] and
            all(peak["trackEvidence"]["replicate"] == "A1"
                for peak in replicate_counts) and
            direct_support["referenceSet"]
            == "replicate_supported_binding_peaks" and
            direct_support["trackEvidence"]["replicate"] == "R1" and
            direct_support["consensusEvidence"]["minimumSupport"] == 2 and
            direct_support["consensusEvidence"]["observedSupport"] == 2 and
            len(direct_support["supportingEvidence"]) == 1 and
            nested_support["referenceSet"] == "strong_myb_binding_rep2" and
            nested_support["trackEvidence"]["replicate"] == "R2" and
            nested_support["trackEvidence"]["condition"]
            == "pigmented petal" and
            nested_support["trackEvidence"]["control"] == "mock" and
            [link["spatialRelation"]["distance"]
             for link in replicate_links] == [25, 60],
            "replicate example preserves direct and nested experimental support",
        )
        replicate_tsv = (
            workspace / "replicate_supported_peak_myb_counts.tsv"
        ).read_text(encoding="utf-8").splitlines()[1].split("\t")
        exported_support = json.loads(replicate_tsv[96])[0]
        require(
            exported_support["trackEvidence"]["replicate"] == "R1" and
            exported_support["consensusEvidence"]["observedSupport"] == 2 and
            exported_support["supportingEvidence"][0]["trackEvidence"]
                ["replicate"] == "R2",
            "TSV recursively exports replicate-support provenance",
        )
        consensus_tsv = (
            workspace / "replicate_supported_binding_peaks.tsv"
        ).read_text(encoding="utf-8").splitlines()[1].split("\t")
        require(
            len(consensus_tsv) == 103 and
            consensus_tsv[97:100] == [
                "strong_myb_binding_rep1", "2", "2",
            ] and
            consensus_tsv[100] == "50" and
            consensus_tsv[101] == "5" and
            json.loads(consensus_tsv[102]) == {
                "anchorSet": "strong_myb_binding_rep1",
                "minimumSupport": 2,
                "observedSupport": 2,
                "minimumReciprocalOverlapPercent": 50,
                "maximumSummitDistanceBp": 5,
                "inputSets": ["strong_myb_binding_rep1",
                              "strong_myb_binding_rep2"],
            },
            "TSV exports typed coordinate-consensus evidence",
        )
        consensus_gff = (
            workspace / "replicate_supported_binding_peaks.gff3"
        ).read_text(encoding="utf-8").splitlines()[1]
        require(
            "ConsensusAnchorSet=strong_myb_binding_rep1" in consensus_gff and
            "ConsensusMinimumSupport=2" in consensus_gff and
            "ConsensusObservedSupport=2" in consensus_gff and
            "ConsensusMinimumReciprocalOverlapPercent=50" in consensus_gff and
            "ConsensusMaximumSummitDistanceBp=5" in consensus_gff and
            "ConsensusEvidenceJSON=" in consensus_gff,
            "GFF3 exports typed coordinate-consensus evidence",
        )

        nearest_example = (
            ROOT / "cql_examples" / "06_nearest_gene_candidates.cql"
        ).read_text(encoding="utf-8")
        data, _ = run_query(
            workspace, "anthocyanin_nearest_gene_candidates", nearest_example
        )
        proximal = data["resultSets"]["gene_proximal_myb_candidates"]
        require([
                    (
                        site["start"],
                        site["spatialRelation"]["reference"]["name"],
                        site["spatialRelation"]["distance"],
                    )
                    for site in proximal
                ] == [
                    (150, "anthocyanin_enzyme", 41),
                    (630, "transporter_candidate", 30),
                ] and
                all(site["motifEvidence"]["matrixId"] == "MA0054.1"
                    for site in proximal),
                "nearest-gene example exposes two auditable MYB candidates")

        count_example = (
            ROOT / "cql_examples" / "05_count_promoter_support.cql"
        ).read_text(encoding="utf-8")
        data, _ = run_query(
            workspace, "anthocyanin_promoter_support_counts", count_example
        )
        require([region["countEvidence"]["count"] for region in
                 data["resultSets"]["promoter_myb_counts"]] == [1, 1, 0] and
                [region["name"] for region in
                 data["resultSets"]["promoters_with_myb_support"]]
                == ["anthocyanin_enzyme_promoter",
                    "transporter_candidate_promoter"] and
                all(region["countEvidence"]["count"] == 1 for region in
                    data["resultSets"]["promoters_with_myb_support"]),
                "promoter-count example retains zeroes and supported promoters")

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

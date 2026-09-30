LOAD TRACK "data_examples/binding_rep1_demo.narrowPeak"
    FORMAT NARROWPEAK
    EVIDENCE BINDING
    ASSAY "synthetic binding assay"
    SAMPLE "demo"
    CONDITION "condition_a"
    REPLICATE "R1"
    CONTROL "mock"
    AS replicate_1;

LOAD TRACK "data_examples/binding_rep2_demo.narrowPeak"
    FORMAT NARROWPEAK
    EVIDENCE BINDING
    ASSAY "synthetic binding assay"
    SAMPLE "demo"
    CONDITION "condition_a"
    REPLICATE "R2"
    CONTROL "mock"
    AS replicate_2;

EXTRACT replicate_1 AS filtered_1
    WHERE TRACK_SCORE >= 800 AND SIGNAL_VALUE >= 18;

EXTRACT replicate_2 AS filtered_2
    WHERE TRACK_SCORE >= 800 AND SIGNAL_VALUE >= 18;

CONSENSUS FROM [filtered_1, filtered_2]
    ANCHOR filtered_1
    MIN_SUPPORT 2
    MIN_RECIPROCAL_OVERLAP 50 %
    MAX_SUMMIT_DISTANCE 5 BP
    AS consensus
    WHERE SUPPORT_COUNT >= 2;

EXPORT consensus TO "consensus.tsv" FORMAT TSV;
EXPORT consensus TO "consensus.gff3" FORMAT GFF3;

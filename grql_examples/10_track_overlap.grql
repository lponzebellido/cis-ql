LOAD TRACK "data_examples/accessibility_demo.narrowPeak"
    FORMAT NARROWPEAK
    EVIDENCE ACCESSIBILITY
    ASSAY "synthetic accessibility assay"
    SAMPLE "demo"
    CONDITION "condition_a"
    REPLICATE "A1"
    CONTROL "input"
    AS accessibility;

LOAD TRACK "data_examples/binding_rep1_demo.narrowPeak"
    FORMAT NARROWPEAK
    EVIDENCE BINDING
    ASSAY "synthetic binding assay"
    SAMPLE "demo"
    CONDITION "condition_a"
    REPLICATE "R1"
    CONTROL "mock"
    AS binding;

EXTRACT accessibility AS open_regions
    WHERE TRACK_SCORE >= 600 AND SIGNAL_VALUE >= 10;

EXTRACT binding AS bound_regions
    WHERE TRACK_SCORE >= 800 AND SIGNAL_VALUE >= 18;

OVERLAPS open_regions WITH bound_regions AS overlaps;

EXPORT overlaps TO "track_overlap.tsv" FORMAT TSV;

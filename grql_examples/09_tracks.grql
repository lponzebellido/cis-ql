LOAD TRACK "data_examples/accessibility_demo.narrowPeak"
    FORMAT NARROWPEAK
    EVIDENCE ACCESSIBILITY
    ASSAY "synthetic accessibility assay"
    SAMPLE "demo"
    CONDITION "condition_a"
    REPLICATE "A1"
    CONTROL "input"
    AS peaks;

EXTRACT peaks AS filtered_peaks
    WHERE TRACK_SCORE >= 600
      AND SIGNAL_VALUE >= 10
      AND CONDITION = "condition_a"
      AND REPLICATE = "A1";

EXPORT filtered_peaks TO "tracks.tsv" FORMAT TSV;

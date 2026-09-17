LOAD SEQUENCE "data_examples/anthocyanin_regulatory_demo.fasta" AS genome;
LOAD ANNOTATION "data_examples/anthocyanin_regulatory_demo.gff3" AS annotation;
LOAD MATRIX "matrices/MA0054.1_myb.Ph3.pwm" AS myb_matrix;

LOAD TRACK "data_examples/anthocyanin_accessibility_demo.narrowPeak"
    FORMAT NARROWPEAK
    EVIDENCE ACCESSIBILITY
    ASSAY "synthetic ATAC-seq-like fixture"
    SAMPLE "synthetic anthocyanin locus"
    CONDITION "pigmented petal"
    REPLICATE "A1"
    CONTROL "input"
    AS accessibility_peaks;

LOAD TRACK "data_examples/anthocyanin_myb_binding_demo.narrowPeak"
    FORMAT NARROWPEAK
    EVIDENCE BINDING
    ASSAY "synthetic DAP-seq-like fixture"
    SAMPLE "synthetic anthocyanin locus"
    CONDITION "pigmented petal"
    REPLICATE "R1"
    CONTROL "mock"
    AS myb_binding_rep1;

LOAD TRACK "data_examples/anthocyanin_myb_binding_rep2_demo.narrowPeak"
    FORMAT NARROWPEAK
    EVIDENCE BINDING
    ASSAY "synthetic DAP-seq-like fixture"
    SAMPLE "synthetic anthocyanin locus"
    CONDITION "pigmented petal"
    REPLICATE "R2"
    CONTROL "mock"
    AS myb_binding_rep2;

EXTRACT accessibility_peaks
    AS strong_accessibility_peaks
    WHERE TRACK_SCORE >= 600 AND SIGNAL_VALUE >= 10
      AND CONDITION = "pigmented petal" AND REPLICATE = "A1";

EXTRACT myb_binding_rep1
    AS strong_myb_binding_rep1
    WHERE TRACK_SCORE >= 800 AND SIGNAL_VALUE >= 18
      AND CONDITION = "pigmented petal" AND REPLICATE = "R1"
      AND CONTROL = "mock";

EXTRACT myb_binding_rep2
    AS strong_myb_binding_rep2
    WHERE TRACK_SCORE >= 800 AND SIGNAL_VALUE >= 18
      AND CONDITION = "pigmented petal" AND REPLICATE = "R2"
      AND CONTROL = "mock";

CONSENSUS FROM [strong_myb_binding_rep1, strong_myb_binding_rep2]
    ANCHOR strong_myb_binding_rep1
    MIN_SUPPORT 2
    AS replicate_supported_binding_peaks
    WHERE SUPPORT_COUNT >= 2;

OVERLAPS strong_accessibility_peaks WITH replicate_supported_binding_peaks
    AS replicate_supported_accessible_peaks;

SCAN myb_matrix
    BACKGROUND FROM genome
    QVALUE <= 0.01
    AS supported_myb_sites;

COUNT supported_myb_sites IN replicate_supported_accessible_peaks
    AS replicate_supported_peak_myb_counts;

EXTRACT replicate_supported_peak_myb_counts
    AS replicate_supported_regulatory_candidates
    WHERE COUNT >= 1;

NEAR replicate_supported_regulatory_candidates TO GENE WITHIN 100 BP
    AS replicate_supported_gene_hypotheses;

EXPORT replicate_supported_peak_myb_counts
    TO "replicate_supported_peak_myb_counts.tsv" FORMAT TSV;
EXPORT replicate_supported_binding_peaks
    TO "replicate_supported_binding_peaks.tsv" FORMAT TSV;
EXPORT replicate_supported_binding_peaks
    TO "replicate_supported_binding_peaks.gff3" FORMAT GFF3;
EXPORT replicate_supported_gene_hypotheses
    TO "replicate_supported_gene_hypotheses.gff3" FORMAT GFF3;

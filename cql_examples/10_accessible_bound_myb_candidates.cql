LOAD SEQUENCE "data_examples/anthocyanin_regulatory_demo.fasta" AS genome;
LOAD ANNOTATION "data_examples/anthocyanin_regulatory_demo.gff3" AS annotation;
LOAD MATRIX "matrices/MA0054.1_myb.Ph3.pwm" AS myb_matrix;

LOAD TRACK "data_examples/anthocyanin_accessibility_demo.narrowPeak"
    FORMAT NARROWPEAK
    EVIDENCE ACCESSIBILITY
    ASSAY "synthetic ATAC-seq-like fixture"
    SAMPLE "synthetic anthocyanin locus"
    AS accessibility_peaks;

LOAD TRACK "data_examples/anthocyanin_myb_binding_demo.narrowPeak"
    FORMAT NARROWPEAK
    EVIDENCE BINDING
    ASSAY "synthetic DAP-seq-like fixture"
    SAMPLE "synthetic anthocyanin locus"
    AS myb_binding_peaks;

EXTRACT accessibility_peaks
    AS strong_accessibility_peaks
    WHERE TRACK_SCORE >= 600 AND SIGNAL_VALUE >= 10;

EXTRACT myb_binding_peaks
    AS strong_myb_binding_peaks
    WHERE TRACK_SCORE >= 800 AND SIGNAL_VALUE >= 18;

OVERLAPS strong_accessibility_peaks WITH strong_myb_binding_peaks
    AS accessible_bound_peaks;

SCAN myb_matrix
    BACKGROUND FROM genome
    QVALUE <= 0.01
    AS supported_myb_sites;

COUNT supported_myb_sites IN accessible_bound_peaks
    AS accessible_bound_peak_myb_counts;

EXTRACT accessible_bound_peak_myb_counts
    AS multi_evidence_regulatory_candidates
    WHERE COUNT >= 1;

NEAR multi_evidence_regulatory_candidates TO GENE WITHIN 100 BP
    AS multi_evidence_gene_hypotheses;

EXPORT accessible_bound_peak_myb_counts
    TO "accessible_bound_peak_myb_counts.tsv" FORMAT TSV;
EXPORT multi_evidence_gene_hypotheses
    TO "multi_evidence_gene_hypotheses.gff3" FORMAT GFF3;

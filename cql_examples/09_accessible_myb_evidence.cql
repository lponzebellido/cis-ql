LOAD SEQUENCE "data_examples/anthocyanin_regulatory_demo.fasta" AS genome;
LOAD ANNOTATION "data_examples/anthocyanin_regulatory_demo.gff3" AS annotation;
LOAD MATRIX "matrices/MA0054.1_myb.Ph3.pwm" AS myb_matrix;
LOAD TRACK "data_examples/anthocyanin_accessibility_demo.narrowPeak"
    FORMAT NARROWPEAK
    EVIDENCE ACCESSIBILITY
    ASSAY "synthetic ATAC-seq-like fixture"
    SAMPLE "synthetic anthocyanin locus"
    AS accessibility_peaks;

SCAN myb_matrix
    BACKGROUND FROM genome
    QVALUE <= 0.01
    AS supported_myb_sites;

COUNT supported_myb_sites IN accessibility_peaks
    AS accessible_peak_myb_counts;

EXTRACT accessible_peak_myb_counts
    AS accessible_peaks_with_myb_support
    WHERE COUNT >= 1;

NEAR accessible_peaks_with_myb_support TO GENE WITHIN 100 BP
    AS accessible_myb_candidate_links;

EXPORT accessible_peak_myb_counts
    TO "accessible_peak_myb_counts.tsv" FORMAT TSV;
EXPORT accessible_myb_candidate_links
    TO "accessible_myb_candidate_links.gff3" FORMAT GFF3;

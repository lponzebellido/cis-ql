LOAD SEQUENCE "data_examples/anthocyanin_regulatory_demo.fasta" AS genome;
LOAD MATRIX "matrices/MA0054.1_myb.Ph3.pwm" AS myb_matrix;

SCAN myb_matrix
    BACKGROUND FROM genome
    QVALUE <= 0.01
    AS statistically_supported_myb_sites;

EXPORT statistically_supported_myb_sites
    TO "myb_sites.tsv" FORMAT TSV;

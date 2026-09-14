LOAD SEQUENCE "data_examples/anthocyanin_regulatory_demo.fasta" AS genome;
LOAD MATRIX "matrices/MA0054.1_myb.Ph3.pwm" AS myb_matrix;

SCAN myb_matrix
    THRESHOLD 90 %
    AS high_scoring_myb_sites;

EXPORT high_scoring_myb_sites TO "high_scoring_myb_sites.tsv" FORMAT TSV;

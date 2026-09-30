LOAD SEQUENCE "data_examples/regulatory_demo.fasta" AS genome;
LOAD MATRIX "matrices/MA0054.1_myb.Ph3.pwm" AS matrix;

SCAN matrix
    THRESHOLD 90 %
    AS hits;

EXPORT hits TO "pwm_hits.tsv" FORMAT TSV;

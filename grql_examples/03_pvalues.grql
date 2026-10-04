LOAD SEQUENCE "data_examples/regulatory_demo.fasta" AS genome;
LOAD MATRIX "matrices/MA0054.1_myb.Ph3.pwm" AS matrix;

SCAN matrix
    BACKGROUND FROM genome
    QVALUE <= 0.01
    AS hits;

EXPORT hits TO "pvalue_hits.tsv" FORMAT TSV;

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tbp_matrix;

SCAN tbp_matrix STRAND POSITIVE
    BACKGROUND FROM genome
    PVALUE <= 1e-4
    AS candidate_tbp_sites;

EXTRACT candidate_tbp_sites;

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tbp_matrix;

SCAN tbp_matrix STRAND POSITIVE THRESHOLD 80 % AS tbp_sites;

EXTRACT tbp_sites;

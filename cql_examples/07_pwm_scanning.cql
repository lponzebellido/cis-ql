/*
 * 07_pwm_scanning.cql
 *
 * This script performs probabilistic motif scanning using a Position Weight Matrix (PWM).
 * A JASPAR format matrix for TBP (MA0108.1) is loaded using LOAD MATRIX. The SCAN statement 
 * evaluates log-odds scores along the sequence and retains hits exceeding 80% of the matrix's 
 * maximum theoretical score.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tbp_matrix;

SCAN tbp_matrix STRAND POSITIVE THRESHOLD 80 % AS tbp_sites;

EXTRACT tbp_sites;

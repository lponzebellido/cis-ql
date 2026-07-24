/*
 * 16_foreach_batch_scan.cql
 *
 * This script demonstrates batch iteration using the FOREACH statement.
 * The loop iterates through a collection of loaded PWM matrices (tbp, sp1, ctcf), 
 * executing a scan at 80% threshold for each matrix in sequence without repeating code blocks.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tbp;
LOAD MATRIX "matrices/MA0079.5_SP1.pwm" AS sp1;
LOAD MATRIX "matrices/MA0139.1_CTCF.pwm" AS ctcf;

FOREACH m IN [tbp, sp1, ctcf] DO
    SCAN m THRESHOLD 80 % AS tf_sites;
ENDFOR;

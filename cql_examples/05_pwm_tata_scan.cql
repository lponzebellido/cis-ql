/*
 * 05_pwm_tata_scan.cql
 * 
 * Purpose: Scans the genome using a Position Weight Matrix (PWM) for the TATA-binding protein (TBP).
 * Expectation: Binding sites that match the TBP probability model with at least 85% confidence.
 * Why: Transcription factor binding sites degenerate over evolutionary time and rarely match an exact string. PWMs are the standard way to score and identify these fuzzy binding sites.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tbp_matrix;

SCAN tbp_matrix THRESHOLD 85 % AS tbp_sites;

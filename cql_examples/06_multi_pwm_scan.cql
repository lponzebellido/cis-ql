/*
 * 06_multi_pwm_scan.cql
 * 
 * Purpose: Scans for multiple transcription factors simultaneously using different PWMs.
 * Expectation: Distinct tracks for CTCF (insulator) and SP1 (activator) binding sites across the genome.
 * Why: Complex gene expression is often controlled by multiple transcription factors. Scanning for several factors at once is the first step in identifying co-regulatory modules.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD MATRIX "matrices/MA0139.1_CTCF.pwm" AS ctcf_matrix;
LOAD MATRIX "matrices/MA0079.5_SP1.pwm" AS sp1_matrix;

SCAN ctcf_matrix THRESHOLD 85 % AS ctcf_sites;
SCAN sp1_matrix THRESHOLD 80 % AS sp1_sites;

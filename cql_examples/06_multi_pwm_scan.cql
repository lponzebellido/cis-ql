// 06_multi_pwm_scan.cql
// Scans the genome with multiple TF matrices to profile the cis-regulatory landscape.
// Expect: binding site predictions for TBP, SP1, and CTCF across both strands.

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

// TBP: TATA-box binding protein — marks core promoter elements
LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tata_pwm;

// SP1: GC-box binding factor — common in housekeeping gene promoters
LOAD MATRIX "matrices/MA0079.5_SP1.pwm" AS sp1_pwm;

// CTCF: chromatin boundary / insulator element
LOAD MATRIX "matrices/MA0139.1_CTCF.pwm" AS ctcf_pwm;

// Scan each matrix on both strands; adjust threshold per factor confidence
SCAN tata_pwm THRESHOLD 80% AS tata_sites;
SCAN sp1_pwm  THRESHOLD 75% AS gc_boxes;
SCAN ctcf_pwm THRESHOLD 85% AS ctcf_sites;

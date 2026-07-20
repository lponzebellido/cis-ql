/*
 * 18_chromatin_insulator_analysis.cql
 * 
 * Purpose: Identifies potential chromatin boundary elements by combining PWM scans for CTCF with CpG island detection.
 * Expectation: Subsets of CTCF binding sites that overlap with CpG islands, and weak/cooperative CTCF sites that are non-canonical.
 * Why: CTCF sites inside CpG islands often have different methylation dynamics than those outside. Subtracting high-confidence sites from relaxed ones helps researchers find weak but cooperative binding events.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;
LOAD MATRIX "matrices/MA0139.1_CTCF.pwm" AS ctcf_pwm;

ANALYZE CPG_ISLANDS AS cpg_islands;

SCAN ctcf_pwm THRESHOLD 90 % AS ctcf_high_confidence;
SCAN ctcf_pwm THRESHOLD 75 % AS ctcf_relaxed;

INTERSECT ctcf_relaxed AND cpg_islands;

EXCEPT ctcf_high_confidence FROM ctcf_relaxed;

/*
 * 10_pwm_plus_set_ops.cql
 * 
 * Purpose: Combines probabilistic PWM scanning with set subtraction to filter results.
 * Expectation: Strong TBP binding sites that are located outside of any annotated CDS regions.
 * Why: This models a typical bioinformatics workflow: finding putative binding sites mathematically, then using structural annotations to discard those located in unlikely regions (like inside protein-coding genes).
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;
LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tbp_matrix;

SCAN tbp_matrix THRESHOLD 85 % AS all_tbp_sites;

EXCEPT all_tbp_sites FROM CDS;

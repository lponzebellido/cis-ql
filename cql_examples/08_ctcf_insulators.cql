/*
 * 08_ctcf_insulators.cql
 *
 * This script combines PWM scanning with set subtraction.
 * Candidate CTCF binding sites are identified genome-wide using a JASPAR matrix. The EXCEPT statement 
 * subtracts any sites that overlap annotated Coding Sequences (CDS), isolating candidates in intergenic regions.
 */

LOAD SEQUENCE "data_examples/ecoli2.fna" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;
LOAD MATRIX "matrices/MA0139.1_CTCF.pwm" AS ctcf_matrix;

SCAN ctcf_matrix THRESHOLD 75 % AS ctcf_sites;

EXCEPT ctcf_sites FROM CDS;

/*
 * 08_set_union.cql
 * 
 * Purpose: Merges two separate result sets into a single track using the UNION operator.
 * Expectation: A unified set of all TATAAT and TTGACA exact motif hits across the genome.
 * Why: Union allows you to group different regulatory elements or motifs into a single logical "regulatory regions" track for simplified downstream analysis.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

FIND MOTIF "TATAAT" STRAND POSITIVE AS minus10;
FIND MOTIF "TTGACA" STRAND POSITIVE AS minus35;

UNION minus10 AND minus35;

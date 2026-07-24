/*
 * 13_complex_where.cql
 *
 * This script demonstrates combining multiple filtering properties in a single WHERE clause.
 * The query extracts genes that satisfy two simultaneous constraints: a physical length greater than 800 BP 
 * and a sequence similarity score greater than 65% when compared against reference alignments.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

EXTRACT GENE WHERE LENGTH > 800 BP AND SIMILARITY > 65 %;

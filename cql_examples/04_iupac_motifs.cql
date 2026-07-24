/*
 * 04_iupac_motifs.cql
 *
 * This script searches for nucleotide motifs containing IUPAC degenerate base codes.
 * The pattern "TATAWAW" uses the IUPAC code 'W' (which represents A or T) to match variable 
 * Pribnow box (-10 promoter) sequences. Cis-QL translates IUPAC ambiguity characters into 
 * standard regular expression patterns during execution.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

FIND MOTIF "TATAWAW" STRAND POSITIVE AS tata_boxes;

EXTRACT tata_boxes;

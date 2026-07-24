/*
 * 11_strand_search.cql
 *
 * This script demonstrates strand-specific motif searching.
 * Specifying STRAND POSITIVE or STRAND NEGATIVE restricts motif matching to a single strand of 
 * the double-stranded DNA sequence. Here, Pribnow box motifs (-10 elements) are queried independently 
 * on each strand.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

FIND MOTIF "TATAAT" STRAND POSITIVE AS sense_promoters;
FIND MOTIF "TATAAT" STRAND NEGATIVE AS antisense_promoters;

EXTRACT sense_promoters;

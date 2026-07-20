/*
 * 11_strand_specific_search.cql
 * 
 * Purpose: Demonstrates how to restrict searches to specific DNA strands.
 * Expectation: Two distinct tracks for the TATAAT motif, one for the positive (forward) strand and one for the negative (reverse) strand.
 * Why: DNA is double-stranded and transcription is directional. Knowing which strand a motif is on is critical for predicting the direction of transcription for a gene.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

FIND MOTIF "TATAAT" STRAND POSITIVE AS pos_strand_promoters;
FIND MOTIF "TATAAT" STRAND NEGATIVE AS neg_strand_promoters;

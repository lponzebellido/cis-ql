/*
 * 02_literal_motif_search.cql
 * 
 * Purpose: Searches for exact literal motifs (the classic sigma-70 -10 and -35 promoter boxes) in the genome.
 * Expectation: Exact-match genomic coordinates for 'TATAAT' and 'TTGACA' on the positive strand.
 * Why: Literal motif searching is the simplest way to find known, highly conserved regulatory sequences without using probabilistic models.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

FIND MOTIF "TATAAT" STRAND POSITIVE AS minus10_box;
FIND MOTIF "TTGACA" STRAND POSITIVE AS minus35_box;

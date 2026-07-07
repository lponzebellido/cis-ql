// 11_strand_specific_search.cql
// Demonstrates strand-specific motif searches across both strands independently.
// Expect: separate tracks for the motif on the plus and minus strands.

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

// SOS response LexA binding box (SOS box): CTGT(N8)ACAG
FIND MOTIF "CTGTN{8}ACAG" STRAND POSITIVE AS lexA_fwd;
FIND MOTIF "CTGTN{8}ACAG" STRAND NEGATIVE AS lexA_rev;

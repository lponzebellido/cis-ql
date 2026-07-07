// 02_literal_motif_search.cql
// Searches for the classic sigma-70 promoter elements (-10 and -35 boxes) in E. coli.
// Expect: exact-match hits for TATAAT and TTGACA across the genome.

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

// The -10 box (Pribnow box) on the positive strand
FIND MOTIF "TATAAT" STRAND POSITIVE AS minus10_box;

// The -35 element on the positive strand
FIND MOTIF "TTGACA" STRAND POSITIVE AS minus35_box;

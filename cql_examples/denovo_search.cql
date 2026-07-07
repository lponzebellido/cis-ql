// denovo_search.cql
// De novo discovery of sigma-70 -10 box motifs with no annotation file required.
// Works directly on the raw FASTA sequence. Expect: hits for TATAAT on the plus strand.

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

// Exact TATAAT search on positive strand only
FIND MOTIF "TATAAT" STRAND POSITIVE AS sigma70_minus10;

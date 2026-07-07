// 08_set_union.cql
// Merges two motif search results into a single combined feature track using UNION.
// Expect: a unified set of all TATAAT and TTGACA hits across the genome.

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

// Scan for both canonical sigma-70 promoter elements
FIND MOTIF "TATAAT" STRAND POSITIVE AS minus10;
FIND MOTIF "TTGACA" STRAND POSITIVE AS minus35;

// Combine them into one track for joint visualization
UNION minus10 AND minus35;

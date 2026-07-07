// 09_set_except.cql
// Filters a result set by subtracting regions found in another set using EXCEPT.
// Expect: promoter motif hits that do NOT fall within annotated exon regions.

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

// Find all Pribnow box occurrences on the positive strand
FIND MOTIF "TATAAT" STRAND POSITIVE AS pribnow_hits;

// Remove any hits that overlap with exon regions — likely false positives
EXCEPT pribnow_hits FROM EXON;

/*
 * 09_set_except.cql
 * 
 * Purpose: Filters a set of regions by subtracting regions found in another set using EXCEPT.
 * Expectation: Pribnow box motif hits that do not fall within any annotated exon regions.
 * Why: Subtracting sets is useful for filtering out false positives. A promoter motif found inside an exon is less likely to be a functional promoter than one found in an intergenic region.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

FIND MOTIF "TATAAT" STRAND POSITIVE AS pribnow_hits;

EXCEPT pribnow_hits FROM EXON;

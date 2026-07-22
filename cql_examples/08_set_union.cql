/*
 * 08_set_union.cql
 *
 * What it does: Merges the -10 box (TATAAT) and -35 box (TTGACA) hits found upstream of
 *               all CDS regions into a single unified track of sigma-70 promoter elements.
 *
 * Expected output: A combined track containing all -10 and -35 hits upstream of CDS,
 *                  with each region labeled by its source (minus35_all or pribnow_all).
 *
 * Why it matters: The two elements of a sigma-70 promoter are separated by ~17 bp and
 *                 are always searched independently. Merging them into one track lets
 *                 you visualize the complete promoter picture in one view without
 *                 switching between result sets.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

FIND MOTIF "TTGACA" WITHIN 200 BP UPSTREAM FROM CDS AS minus35_all;

FIND MOTIF "TATAAT" WITHIN 50 BP UPSTREAM FROM CDS AS pribnow_all;

UNION minus35_all AND pribnow_all;

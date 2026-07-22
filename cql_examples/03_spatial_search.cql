/*
 * 03_spatial_search.cql
 *
 * What it does: Searches for the sigma-70 -10 box (TATAAT) in the region immediately
 *               upstream of two specific genes: thrA and thrL, which are the first genes
 *               of the threonine biosynthesis operon in E. coli.
 *
 * Expected output: A small number of TATAAT hits within 150 bp upstream of each gene.
 *                  Results are strand-specific to the gene's orientation.
 *
 * Why it matters: Regulatory elements are spatial — they only function when positioned
 *                 near the right gene. Restricting a search to the upstream window of a
 *                 specific gene eliminates genome-wide noise and focuses on functional candidates.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

FIND MOTIF "TATAAT" WITHIN 150 BP UPSTREAM FROM GENE "thrA" AS thrA_minus10;

FIND MOTIF "TATAAT" WITHIN 150 BP UPSTREAM FROM GENE "thrL" AS thrL_minus10;

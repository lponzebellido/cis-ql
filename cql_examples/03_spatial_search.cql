/*
 * 03_spatial_search.cql
 * 
 * Purpose: Searches for a specific motif restricted to the upstream region of specific genes.
 * Expectation: Motif hits located strictly within 500 bp upstream of the 'thrA' and 'thrL' genes.
 * Why: Regulatory elements typically reside upstream of the gene they regulate. Spatial constraints reduce false positives by only searching biologically relevant regions.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

FIND MOTIF "TATAAT" WITHIN 500 BP UPSTREAM FROM GENE "thrA" AS thrA_promoter;
FIND MOTIF "TATAAT" WITHIN 500 BP UPSTREAM FROM GENE "thrL" AS thrL_promoter;

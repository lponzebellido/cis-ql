// 03_spatial_search.cql
// Finds the sigma-70 -10 box (TATAAT) in the upstream region of specific genes.
// Expect: motif hits located within 500 bp upstream of thrA and thrL.

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

// Look for the Pribnow box upstream of the thrA operon start
FIND MOTIF "TATAAT" WITHIN 500 BP UPSTREAM FROM GENE "thrA" AS thrA_promoter;

// Same search for thrL (the leader peptide of the thr operon)
FIND MOTIF "TATAAT" WITHIN 500 BP UPSTREAM FROM GENE "thrL" AS thrL_promoter;

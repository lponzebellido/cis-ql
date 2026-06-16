LOAD SEQUENCE "examples/test.fasta" AS ambiente_1;

FIND MOTIF "TATA[AT]A[AT]" STRAND POSITIVE AS putative_promoters;

FIND MOTIF "ATG(...)*?(TAA|TAG|TGA)" 
WITHIN 200 BP DOWNSTREAM FROM putative_promoters 
AS putative_genes 
WHERE LENGTH >= 15 BP;

EXTRACT putative_genes;


LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

FIND MOTIF "TATA[AT]A[AT]" STRAND POSITIVE AS putative_promoters;

FIND MOTIF "ATG(...)*?(TAA|TAG|TGA)"
  WITHIN 200 BP DOWNSTREAM FROM putative_promoters
  AS putative_orfs
  WHERE LENGTH >= 15 BP;

EXTRACT putative_orfs;



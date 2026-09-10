LOAD SEQUENCE "data_examples/ecoli.fasta" AS sample;

FIND MOTIF "TATA[AT]A[AT]" STRAND POSITIVE AS putative_promoters;

FIND MOTIF "ATG(...)*?(TAA|TAG|TGA)"
    WITHIN 2500 BP DOWNSTREAM FROM putative_promoters
    AS candidate_orfs
    WHERE LENGTH > 600 BP;

EXTRACT candidate_orfs;

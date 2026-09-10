LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

FIND MOTIF "TATAAT" STRAND POSITIVE AS sense_promoters;
FIND MOTIF "TATAAT" STRAND NEGATIVE AS antisense_promoters;

EXTRACT sense_promoters;

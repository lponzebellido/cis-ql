LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

FIND MOTIF "TATAWAW" STRAND POSITIVE AS tata_boxes;

EXTRACT tata_boxes;

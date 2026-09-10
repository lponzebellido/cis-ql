LOAD SEQUENCE "data_examples/ecoli2.fna" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

FIND MOTIF "TTGACA" WITHIN 200 BP UPSTREAM FROM CDS AS minus35_box;
FIND MOTIF "TATAAT" WITHIN 50 BP UPSTREAM FROM CDS AS minus10_box;

UNION minus35_box AND minus10_box;

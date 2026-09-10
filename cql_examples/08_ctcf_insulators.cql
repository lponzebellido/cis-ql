LOAD SEQUENCE "data_examples/ecoli2.fna" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;
LOAD MATRIX "matrices/MA0139.1_CTCF.pwm" AS ctcf_matrix;

SCAN ctcf_matrix THRESHOLD 75 % AS ctcf_sites;

EXCEPT ctcf_sites FROM CDS;

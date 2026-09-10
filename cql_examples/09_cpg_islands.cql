LOAD SEQUENCE "data_examples/ecoli2.fna" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;
LOAD MATRIX "matrices/MA0079.5_SP1.pwm" AS sp1;

ANALYZE CPG_ISLANDS AS cpg_islands;
SCAN sp1 THRESHOLD 80 % AS sp1_sites;

INTERSECT sp1_sites AND cpg_islands;

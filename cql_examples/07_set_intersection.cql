// 07_set_intersection.cql
// Identifies genomic regions shared between two feature sets using INTERSECT.
// Expect: genes that overlap with annotated CDS regions (high similarity >= 90%).

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

// Find all regions where a gene and its CDS overlap with at least 90% similarity
INTERSECT GENE AND CDS WHERE SIMILARITY >= 90.0%;

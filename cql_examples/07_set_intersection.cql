/*
 * 07_set_intersection.cql
 * 
 * Purpose: Finds regions where two features overlap using the INTERSECT operator.
 * Expectation: A set of regions representing where annotated genes overlap with annotated CDS (coding sequences).
 * Why: Intersection is crucial for identifying functionally critical overlapping regions, such as genes that encode proteins versus non-coding genes.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

INTERSECT GENE AND CDS;

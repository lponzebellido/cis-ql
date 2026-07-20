/*
 * 13_similarity_filter.cql
 * 
 * Purpose: Filters extracted features based on their sequence similarity to the first feature in the set.
 * Expectation: A list of exon sequences that have at least 80% similarity to the first annotated exon.
 * Why: Similarity filtering (using Smith-Waterman alignment) is essential for finding families of related genes or elements that may share a common evolutionary origin but have diverged in sequence.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

EXTRACT EXON WHERE SIMILARITY >= 80 %;

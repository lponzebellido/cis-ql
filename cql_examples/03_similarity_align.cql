/*
 * 03_similarity_align.cql
 *
 * This script applies a sequence similarity filter using local pairwise alignment.
 * The WHERE SIMILARITY clause evaluates candidate features against a reference sequence using 
 * the Smith-Waterman local alignment algorithm in C++. Here, extracted genes are retained 
 * only if their sequence similarity score exceeds 70%.
 */

LOAD SEQUENCE "data_examples/ecoli2.fna" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

EXTRACT GENE WHERE LENGTH > 1 KB AND SIMILARITY > 70 %;

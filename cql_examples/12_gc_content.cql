/*
12_gc_content.cql

This script profiles GC content variation across the genome using a sliding window.
 * The ANALYZE GC_CONTENT WINDOW 1 KB statement calculates the percentage of G and C bases 
 * within non-overlapping 1000 BP windows. The resulting dataset can be rendered as a continuous track.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

ANALYZE GC_CONTENT WINDOW 1 KB AS gc_profile;

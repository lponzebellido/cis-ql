/*
 * 14_gc_analysis.cql
 * 
 * Purpose: Performs a genome-wide analysis of GC content and detects CpG islands.
 * Expectation: A continuous GC content profile using 100 bp sliding windows, and a set of regions classified as CpG islands.
 * Why: GC content variations are heavily linked to functional regions in a genome, and CpG islands are critical markers for gene regulation and promoters in many species.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

ANALYZE GC_CONTENT WINDOW 100 BP AS ecoli_gc_profile;
ANALYZE CPG_ISLANDS AS cpg_islands;

EXTRACT PROMOTER WHERE LENGTH > 100 BP;

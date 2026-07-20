/*
 * 17_coregulated_operon_discovery.cql
 * 
 * Purpose: Simulates a workflow to identify genes co-regulated by a specific mechanism (e.g., heat shock) by combining core promoter elements and regulatory motifs.
 * Expectation: A highly refined set of putative repressor binding sites located in active promoter regions.
 * Why: Genes are often regulated not by a single motif, but by a specific arrangement of motifs. This script demonstrates how to chain multiple searches and set operations to find complex regulatory modules.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

ANALYZE GC_CONTENT WINDOW 500 BP AS gc_profile;

FIND MOTIF "AGGAGG" WITHIN 50 BP UPSTREAM FROM CDS AS shine_dalgarno;
FIND MOTIF "TATAAT" WITHIN 100 BP UPSTREAM FROM CDS AS minus10_boxes;

UNION shine_dalgarno AND minus10_boxes;

FIND MOTIF "CTGT[AT]{4,8}ACAG" WITHIN 150 BP UPSTREAM FROM CDS AS repressor_sites;

EXTRACT repressor_sites WHERE LENGTH > 12 BP;

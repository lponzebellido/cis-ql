/*
 * 15_cpg_promoter_intersection.cql
 * 
 * Purpose: Identifies which annotated promoters overlap directly with high-confidence CpG islands.
 * Expectation: A filtered list of strong candidate promoter regions that coincide with CpG islands.
 * Why: In many genomes (especially vertebrates), a large percentage of active promoters are associated with CpG islands. Intersecting these two sets yields highly confident regulatory regions.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

ANALYZE CPG_ISLANDS AS all_cpg_islands;

EXTRACT PROMOTER WHERE LENGTH > 100 BP;

INTERSECT all_cpg_islands AND PROMOTER;

/*
 * 09_cpg_islands.cql
 *
 * This script identifies CpG islands and computes their overlap with transcription factor binding sites.
 * The ANALYZE CPG_ISLANDS statement calculates GC content and observed-to-expected CpG ratios. 
 * The INTERSECT statement then identifies SP1 binding sites that co-occur within these CpG islands.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;
LOAD MATRIX "matrices/MA0079.5_SP1.pwm" AS sp1;

ANALYZE CPG_ISLANDS AS cpg_islands;
SCAN sp1 THRESHOLD 80 % AS sp1_sites;

INTERSECT sp1_sites AND cpg_islands;

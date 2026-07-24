/*
 * 14_integrated_query.cql
 *
 * This script integrates multiple analysis types into a single query pipeline.
 * It computes genome-wide GC content and CpG islands, searches for upstream promoter elements 
 * (-35, -10) and ribosome binding sites (Shine-Dalgarno: AGGAGG), scans for TBP and SP1 binding sites, 
 * excludes coding region overlaps using EXCEPT, and extracts structural genes over 2 KB.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;
LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tbp;
LOAD MATRIX "matrices/MA0079.5_SP1.pwm" AS sp1;

ANALYZE GC_CONTENT WINDOW 1 KB AS gc_landscape;
ANALYZE CPG_ISLANDS AS cpg_islands;

FIND MOTIF "TTGACA" WITHIN 200 BP UPSTREAM FROM CDS AS minus35;
FIND MOTIF "TATAAT" WITHIN 50 BP UPSTREAM FROM CDS AS minus10;
FIND MOTIF "AGGAGG" WITHIN 30 BP UPSTREAM FROM CDS AS shine_dalgarno;

UNION minus35 AND minus10;

SCAN tbp THRESHOLD 85 % AS tbp_sites;
SCAN sp1 THRESHOLD 80 % AS sp1_sites;

EXCEPT tbp_sites FROM CDS;
EXCEPT sp1_sites FROM CDS;

EXTRACT GENE WHERE LENGTH > 2 KB;

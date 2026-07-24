/*
 * 10_promoter_union.cql
 *
 * This script merges distinct feature tracks using the UNION set operator.
 * Two motif searches identify bacterial promoter elements: -35 boxes (TTGACA) within 200 BP upstream 
 * of CDS, and -10 boxes (TATAAT) within 50 BP upstream of CDS. The UNION statement consolidates 
 * both feature sets into a single result track.
 */

LOAD SEQUENCE "data_examples/ecoli2.fna" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

FIND MOTIF "TTGACA" WITHIN 200 BP UPSTREAM FROM CDS AS minus35_box;
FIND MOTIF "TATAAT" WITHIN 50 BP UPSTREAM FROM CDS AS minus10_box;

UNION minus35_box AND minus10_box;

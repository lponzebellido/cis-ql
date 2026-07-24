/*
 * 05_spatial_promoters.cql
 *
 * This script restricts motif searches to a specific genomic region relative to annotated features.
 * The WITHIN 200 BP UPSTREAM FROM CDS clause constrains the search for the bacterial -35 consensus 
 * sequence (TTGACA) to the 200 base pair window preceding each Coding Sequence (CDS).
 */

LOAD SEQUENCE "data_examples/ecoli2.fna" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

FIND MOTIF "TTGACA" WITHIN 200 BP UPSTREAM FROM CDS AS minus35_promoters;

EXTRACT minus35_promoters;

/*
 * 01_extract_genes.cql
 *
 * This script extracts annotated gene features based on a minimum length threshold.
 * The EXTRACT statement filters genomic entities loaded from a GFF3 file. Here, the WHERE LENGTH > 1.5 KB
 * clause retains only genes longer than 1500 base pairs. Supported units include BP, KB, and MB.
 */

LOAD SEQUENCE "data_examples/ecoli2.fna" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

EXTRACT GENE WHERE LENGTH > 1.5 KB;

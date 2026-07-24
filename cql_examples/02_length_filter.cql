/*
 * 02_length_filter.cql
 *
 * This script demonstrates combining range conditions in a WHERE clause.
 * Using relational operators (>=, <=, >, <, =) and logical operators (AND, OR), queries can specify 
 * bounded intervals. In this query, genes are filtered to include only those with lengths between 
 * 500 BP and 3 KB.
 */

LOAD SEQUENCE "data_examples/ecoli2.fna" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

EXTRACT GENE WHERE LENGTH >= 500 BP AND LENGTH <= 3 KB;

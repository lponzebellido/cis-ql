/*
 * 01_load_and_extract.cql
 * 
 * Purpose: Demonstrates how to load a genome and extract specific features by size.
 * Expectation: Lists of genes, coding sequences (CDS), and promoters that meet the length criteria.
 * Why: This is a fundamental operation in bioinformatics to filter out small, likely false-positive features or to isolate features of a specific length for downstream analysis.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

EXTRACT GENE WHERE LENGTH > 1000 BP;
EXTRACT CDS WHERE LENGTH > 300 BP AND LENGTH < 1000 BP;
EXTRACT PROMOTER WHERE LENGTH > 100 BP;

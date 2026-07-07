// 01_load_and_extract.cql
// Loads an annotated E. coli genome and extracts genomic features by type and length.
// Expect: lists of genes, CDSs, and promoters that meet the size filters.

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

// All genes longer than 1 kb
EXTRACT GENE WHERE LENGTH > 1000 BP;

// Coding sequences between 300 bp and 1 kb
EXTRACT CDS WHERE LENGTH > 300 BP AND LENGTH < 1000 BP;

// Promoter regions longer than 100 bp
EXTRACT PROMOTER WHERE LENGTH > 100 BP;

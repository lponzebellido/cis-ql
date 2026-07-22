/*
 * 01_load_and_extract.cql
 *
 * What it does: Loads the E. coli genome and annotation, then filters annotated features
 *               by size. Each EXTRACT statement targets a different feature type.
 *
 * Expected output: ~586 genes > 1500 bp, ~4000 CDS between 150 and 600 bp,
 *                  and all annotated promoter regions longer than 200 bp.
 *
 * Why it matters: A length filter is one of the most basic operations in genomics.
 *                 Short ORFs are often noise or artifacts; applying a minimum length
 *                 threshold quickly narrows a search to biologically relevant features.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

EXTRACT GENE WHERE LENGTH > 1500 BP;

EXTRACT CDS WHERE LENGTH > 150 BP AND LENGTH < 600 BP;

EXTRACT PROMOTER WHERE LENGTH > 200 BP;

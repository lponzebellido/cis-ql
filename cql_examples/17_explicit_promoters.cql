/*
 * Explicit promoter construction for the v3 regulatory foundation.
 * No promoter length is inferred while loading the annotation.
 */

LOAD SEQUENCE "data_examples/ecoli2.fna" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annotation;

DEFINE PROMOTERS OF GENE
  FROM TSS
  UPSTREAM 200 BP
  DOWNSTREAM 25 BP
  AS explicit_promoters;

EXTRACT explicit_promoters;

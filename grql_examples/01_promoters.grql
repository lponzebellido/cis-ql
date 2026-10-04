LOAD SEQUENCE "data_examples/ecoli2.fna" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annotation;

EXTRACT GENE AS thra
  WHERE ID = "gene-b0002";

DEFINE PROMOTERS OF thra
    FROM TSS
    UPSTREAM 60 BP
    DOWNSTREAM 20 BP
    AS promoters;

EXPORT promoters TO "promoters.bed" FORMAT BED;

LOAD SEQUENCE "data_examples/anthocyanin_regulatory_demo.fasta" AS genome;
LOAD ANNOTATION "data_examples/anthocyanin_regulatory_demo.gff3" AS annotation;

DEFINE PROMOTERS OF GENE
  FROM TSS
    UPSTREAM 100 BP
    DOWNSTREAM 20 BP
    AS explicit_promoters;

EXPORT explicit_promoters TO "explicit_promoters.bed" FORMAT BED;

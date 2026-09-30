LOAD SEQUENCE "data_examples/ecoli2.fna" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annotation;

DEFINE PROMOTERS OF GENE
    FROM TSS
    UPSTREAM 60 BP
    DOWNSTREAM 20 BP
    AS promoters;

FIND MOTIF "TATAAT" AS sites;

OVERLAPS sites WITH promoters AS promoter_sites;

EXPORT promoter_sites TO "overlap.tsv" FORMAT TSV;

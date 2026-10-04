LOAD SEQUENCE "data_examples/ecoli2.fna" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annotation;

DEFINE PROMOTERS OF GENE
    FROM TSS
    UPSTREAM 60 BP
    DOWNSTREAM 20 BP
    AS promoters;

FIND MOTIF "TATAAT" AS sites;

COUNT sites IN promoters AS motif_counts;

EXTRACT motif_counts AS nonzero_counts
    WHERE COUNT >= 1;

EXPORT motif_counts TO "counts.tsv" FORMAT TSV;
EXPORT nonzero_counts TO "nonzero_counts.gff3" FORMAT GFF3;

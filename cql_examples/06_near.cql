LOAD SEQUENCE "data_examples/ecoli2.fna" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annotation;

FIND MOTIF "TATAAT" AS sites;

NEAR sites TO GENE WITHIN 30 BP AS nearby_sites;

EXPORT nearby_sites TO "near.tsv" FORMAT TSV;

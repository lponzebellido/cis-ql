LOAD SEQUENCE "data_examples/ecoli2.fna" AS genome;

FIND MOTIF "TATAAT" STRAND POSITIVE AS first_sites;
FIND MOTIF "AGGAGG" STRAND POSITIVE AS second_sites;

DEFINE MODULE FROM first_sites WITH second_sites
    SPACING 5 BP TO 30 BP
    ORDER ANY
    ORIENTATION ANY
    AS modules;

EXPORT modules TO "modules.tsv" FORMAT TSV;
EXPORT modules TO "modules.gff3" FORMAT GFF3;

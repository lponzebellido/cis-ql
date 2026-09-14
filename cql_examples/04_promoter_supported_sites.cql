LOAD SEQUENCE "data_examples/anthocyanin_regulatory_demo.fasta" AS genome;
LOAD ANNOTATION "data_examples/anthocyanin_regulatory_demo.gff3" AS annotation;
LOAD MATRIX "matrices/MA0054.1_myb.Ph3.pwm" AS myb_matrix;

DEFINE PROMOTERS OF GENE
    FROM TSS
    UPSTREAM 100 BP
    DOWNSTREAM 20 BP
    AS candidate_promoters;

SCAN myb_matrix
    BACKGROUND FROM genome
    QVALUE <= 0.01
    AS supported_myb_sites;

OVERLAPS supported_myb_sites WITH candidate_promoters
    AS promoter_supported_myb_sites;

EXPORT promoter_supported_myb_sites
    TO "promoter_supported_myb_sites.tsv" FORMAT TSV;

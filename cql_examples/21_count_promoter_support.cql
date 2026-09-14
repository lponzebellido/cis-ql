LOAD SEQUENCE "data_examples/anthocyanin_regulatory_demo.fasta" AS genome;
LOAD ANNOTATION "data_examples/anthocyanin_regulatory_demo.gff3" AS annotation;
LOAD MATRIX "matrices/MA0054.1_myb.Ph3.pwm" AS myb_matrix;

DEFINE PROMOTERS OF GENE
    FROM TSS
    UPSTREAM 100 BP
    DOWNSTREAM 20 BP
    AS candidate_promoters;

SCAN myb_matrix
    IN candidate_promoters
    BACKGROUND FROM genome
    QVALUE <= 0.01
    AS supported_myb_sites;

COUNT supported_myb_sites IN candidate_promoters
    AS promoter_myb_counts;

EXTRACT promoter_myb_counts AS promoters_with_myb_support
    WHERE COUNT >= 1;

EXPORT promoter_myb_counts TO "promoter_myb_counts.tsv" FORMAT TSV;
EXPORT promoters_with_myb_support
    TO "promoters_with_myb_support.gff3" FORMAT GFF3;

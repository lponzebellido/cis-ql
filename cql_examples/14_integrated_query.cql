LOAD SEQUENCE "data_examples/anthocyanin_regulatory_demo.fasta" AS genome;
LOAD ANNOTATION "data_examples/anthocyanin_regulatory_demo.gff3" AS annotation;
LOAD MATRIX "matrices/MA0054.1_myb.Ph3.pwm" AS myb_matrix;

DEFINE PROMOTERS OF GENE
    FROM TSS
    UPSTREAM 100 BP
    DOWNSTREAM 20 BP
    AS candidate_promoters;

EXTRACT ENHANCER AS candidate_enhancers;

SCAN myb_matrix
    BACKGROUND FROM genome
    QVALUE <= 0.01
    AS supported_myb_sites;

OVERLAPS supported_myb_sites WITH candidate_promoters
    AS promoter_myb_evidence;

OVERLAPS supported_myb_sites WITH candidate_enhancers
    AS enhancer_myb_evidence;

EXPORT promoter_myb_evidence TO "promoter_myb_evidence.tsv" FORMAT TSV;
EXPORT enhancer_myb_evidence TO "enhancer_myb_evidence.tsv" FORMAT TSV;

LOAD SEQUENCE "data_examples/anthocyanin_regulatory_demo.fasta" AS genome;
LOAD ANNOTATION "data_examples/anthocyanin_regulatory_demo.gff3" AS annotation;
LOAD MATRIX "matrices/MA0054.1_myb.Ph3.pwm" AS myb_matrix;

SCAN myb_matrix
    BACKGROUND FROM genome
    QVALUE <= 0.01
    AS supported_myb_sites;

NEAR supported_myb_sites TO GENE WITHIN 50 BP
    AS gene_proximal_myb_candidates;

EXPORT gene_proximal_myb_candidates
    TO "gene_proximal_myb_candidates.tsv" FORMAT TSV;

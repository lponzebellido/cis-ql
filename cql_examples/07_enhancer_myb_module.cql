LOAD SEQUENCE "data_examples/anthocyanin_regulatory_demo.fasta" AS genome;
LOAD ANNOTATION "data_examples/anthocyanin_regulatory_demo.gff3" AS annotation;
LOAD MATRIX "matrices/MA0054.1_myb.Ph3.pwm" AS myb_matrix;

EXTRACT ENHANCER AS candidate_enhancers;

SCAN myb_matrix
    BACKGROUND FROM genome
    QVALUE <= 0.01
    AS supported_myb_sites;

OVERLAPS supported_myb_sites WITH candidate_enhancers
    AS enhancer_myb_sites;

DEFINE MODULE
    FROM enhancer_myb_sites WITH enhancer_myb_sites
    SPACING 5 BP TO 30 BP
    ORDER ANY
    ORIENTATION ANY
    AS enhancer_myb_modules;

EXPORT enhancer_myb_modules TO "enhancer_myb_modules.tsv" FORMAT TSV;
EXPORT enhancer_myb_modules TO "enhancer_myb_modules.gff3" FORMAT GFF3;

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

COUNT supported_myb_sites IN candidate_promoters
    AS promoter_myb_counts;

EXTRACT promoter_myb_counts AS supported_promoter_candidates
    WHERE COUNT >= 1;

DEFINE MODULE
    FROM enhancer_myb_evidence WITH enhancer_myb_evidence
    SPACING 5 BP TO 30 BP
    ORDER ANY
    ORIENTATION ANY
    AS enhancer_myb_modules;

NEAR enhancer_myb_modules TO GENE WITHIN 100 BP
    AS module_nearest_gene_hypotheses;

EXPORT promoter_myb_evidence TO "promoter_myb_evidence.tsv" FORMAT TSV;
EXPORT enhancer_myb_modules TO "enhancer_myb_modules.tsv" FORMAT TSV;
EXPORT promoter_myb_counts TO "promoter_myb_counts.tsv" FORMAT TSV;
EXPORT supported_promoter_candidates
    TO "supported_promoter_candidates.gff3" FORMAT GFF3;
EXPORT module_nearest_gene_hypotheses
    TO "module_nearest_gene_hypotheses.tsv" FORMAT TSV;

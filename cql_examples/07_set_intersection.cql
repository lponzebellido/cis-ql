/*
 * 07_set_intersection.cql
 *
 * What it does: Finds the genomic overlap between two sets: TBP binding sites (from a
 *               PWM scan) and all annotated CDS regions. The result is binding sites that
 *               fall inside protein-coding sequences — likely spurious or internal starts.
 *
 * Expected output: A set of TBP-like sites located within CDS annotations.
 *
 * Why it matters: Intersection is the primary way to ask "which elements from set A
 *                 are also in set B?" In this case, a TBP site inside a CDS is a red flag
 *                 — it suggests the match may be a coding artifact rather than a true
 *                 promoter element.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;
LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tbp;

SCAN tbp THRESHOLD 80 % AS tbp_sites;

INTERSECT tbp_sites AND CDS;

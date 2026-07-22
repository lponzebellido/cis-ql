/*
 * 19_duplicate_gene_families.cql
 *
 * What it does: Extracts all annotated exons from the genome and applies Smith-Waterman
 *               similarity filtering at two thresholds: 70% and 90%.
 *               The 70% threshold identifies broadly related exons from the same gene
 *               family; the 90% threshold captures near-identical copies, indicative
 *               of very recent duplication events.
 *               A third extraction grabs all large genes (> 3000 bp) for comparison.
 *
 * Expected output: A filtered set of exons at each similarity threshold, showing the
 *                  graduated relationship between duplicated sequences. Large genes are
 *                  provided as a reference track.
 *
 * Why it matters: Paralogous gene families are common in bacteria — especially for
 *                 transposases, toxin-antitoxin systems, and metabolic enzymes.
 *                 Graded similarity filtering lets researchers distinguish between
 *                 ancient duplications (low similarity) and recent duplications (high
 *                 similarity), providing insight into genome evolution history.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

EXTRACT EXON WHERE SIMILARITY >= 70 %;

EXTRACT EXON WHERE SIMILARITY >= 90 %;

EXTRACT GENE WHERE LENGTH > 3000 BP;

/*
 * 13_similarity_filter.cql
 *
 * What it does: Extracts all annotated exon regions from the genome, then applies a
 *               sequence similarity filter using Smith-Waterman local alignment.
 *               Only exons with at least 80% similarity to the first exon in the set
 *               are kept.
 *
 * Expected output: A small set of exons that are highly similar in sequence to the
 *                  first annotated exon — these are likely exons from duplicated genes
 *                  or paralogs.
 *
 * Why it matters: Gene duplication is a key mechanism of genome evolution. Filtering
 *                 features by sequence similarity identifies families of related elements
 *                 that share a common ancestor but may have accumulated mutations since.
 *                 Smith-Waterman alignment is the gold standard for this comparison.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

EXTRACT EXON WHERE SIMILARITY >= 80 %;

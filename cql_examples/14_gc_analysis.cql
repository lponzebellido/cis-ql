/*
 * 14_gc_analysis.cql
 *
 * What it does: Performs a genome-wide structural analysis of the E. coli genome.
 *               A sliding 100 bp window computes the GC percentage at every position,
 *               producing a continuous coverage track. A second analysis detects CpG
 *               islands using the Gardiner-Garden & Frommer algorithm.
 *
 * Expected output: ~46,000 windows in the GC content profile. Many CpG island candidates
 *                  (E. coli has uniformly high GC content ~51%, so the algorithm detects
 *                  many regions). A set of all promoters longer than 100 bp is also extracted
 *                  for visual comparison against the GC landscape.
 *
 * Why it matters: GC content is one of the most informative structural features of a genome.
 *                 High-GC windows correlate with coding density and gene expression level.
 *                 CpG islands, in bacteria, mark horizontally transferred regions or
 *                 unusually structured regulatory zones.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

ANALYZE GC_CONTENT WINDOW 100 BP AS gc_100bp;

ANALYZE CPG_ISLANDS AS cpg_islands;

EXTRACT PROMOTER WHERE LENGTH > 100 BP;

/*
 * 18_genome_gc_landscape.cql
 *
 * What it does: Profiles the GC content of the entire E. coli genome at three different
 *               resolutions: 100 bp (fine), 1 KB (medium), and 5 KB (coarse).
 *               Alongside this, extracts annotated genes and CDS as reference tracks
 *               for overlay comparison in the track viewer.
 *
 * Expected output: Three GC content profiles at different window sizes, plus gene and
 *                  CDS annotation tracks for visual context. Large genes (> 2000 bp) are
 *                  also highlighted as a separate track.
 *
 * Why it matters: Genome architecture is not uniform. GC content varies regionally due
 *                 to evolutionary pressures, horizontal gene transfer, and structural
 *                 constraints. Comparing GC profiles at multiple resolutions reveals both
 *                 fine-grained codon biases and large-scale genomic islands or AT-rich
 *                 regions that indicate non-native DNA.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

ANALYZE GC_CONTENT WINDOW 100 BP AS gc_fine;

ANALYZE GC_CONTENT WINDOW 1 KB AS gc_medium;

ANALYZE GC_CONTENT WINDOW 5 KB AS gc_coarse;

EXTRACT GENE WHERE LENGTH > 2000 BP;

EXTRACT CDS WHERE LENGTH > 1000 BP;

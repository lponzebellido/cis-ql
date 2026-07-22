/*
 * 09_set_except.cql
 *
 * What it does: Finds all TBP-like binding sites across the genome (relaxed threshold),
 *               then removes those that fall inside coding sequences (CDS).
 *               What remains are TBP sites in intergenic or non-coding regions —
 *               the most likely candidates for real promoter elements.
 *
 * Expected output: Thousands of TBP-like sites with all hits inside CDS removed,
 *                  leaving only intergenic / regulatory-region candidates.
 *
 * Why it matters: Most of the genome is coding sequence. A broad scan will produce
 *                 many hits inside genes, which are almost certainly not functional
 *                 promoters. Subtracting CDS hits is a simple, powerful filter that
 *                 dramatically improves signal-to-noise.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;
LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tbp;

SCAN tbp THRESHOLD 80 % AS tbp_all;

EXCEPT tbp_all FROM CDS;

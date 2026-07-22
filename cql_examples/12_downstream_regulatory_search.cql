/*
 * 12_downstream_regulatory_search.cql
 *
 * What it does: Searches for two motifs characteristic of rho-independent transcription
 *               terminators in the 200 bp immediately after every CDS.
 *               First: a GC-rich stem pattern (G{2,5}C{2,5}) that forms a hairpin.
 *               Second: the Shine-Dalgarno sequence (AGGAG) in the 20 bp before CDS,
 *               which is shown as a reference track for comparison.
 *
 * Expected output: Dozens of GC-stem candidates downstream of CDS, and ~100 Shine-
 *                  Dalgarno hits upstream of CDS.
 *
 * Why it matters: Rho-independent terminators work by forming an RNA hairpin followed
 *                 by a poly-U run. Scanning downstream of every CDS maps where
 *                 transcription is expected to stop, which helps reconstruct operon
 *                 boundaries without any additional experimental data.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

FIND MOTIF "G{2,5}C{2,5}" WITHIN 200 BP DOWNSTREAM FROM CDS AS gc_stems;

FIND MOTIF "AGGAG" WITHIN 20 BP UPSTREAM FROM CDS AS shine_dalgarno;

EXTRACT shine_dalgarno WHERE LENGTH >= 5 BP;

/*
 * 04_regex_denovo.cql
 *
 * What it does: Uses a regular expression to find degenerate variants of the -10 box
 *               across the entire genome on the positive strand, then looks downstream
 *               of those hits for potential short open reading frames (ORFs) using a
 *               regex that captures ATG...stop codon patterns.
 *
 * Expected output: A set of IUPAC-degenerate -10 box hits, followed by a set of
 *                  putative ORFs downstream of them, filtered to at least 15 bp.
 *
 * Why it matters: Real promoter -10 boxes are not always the perfect TATAAT consensus.
 *                 Using a degenerate regex like TATA[AT]A[AT] recovers additional
 *                 candidates that exact-match search would miss. This is de novo discovery
 *                 — no annotation file required.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

FIND MOTIF "TATA[AT]A[AT]" STRAND POSITIVE AS putative_minus10;

FIND MOTIF "ATG(...)*?(TAA|TAG|TGA)"
    WITHIN 200 BP DOWNSTREAM FROM putative_minus10
    AS candidate_orfs
    WHERE LENGTH >= 15 BP;

EXTRACT candidate_orfs;

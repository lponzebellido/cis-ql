/*
 * 02_literal_motif_search.cql
 *
 * What it does: Searches for the two core elements of a bacterial sigma-70 promoter:
 *               the -10 box (Pribnow box, consensus TATAAT) and the -35 box (consensus TTGACA),
 *               both on the positive strand.
 *
 * Expected output: Hundreds of exact TATAAT hits and TTGACA hits across the whole genome,
 *                  returned as separate named tracks.
 *
 * Why it matters: The sigma-70 factor controls transcription of most housekeeping genes in E. coli.
 *                 Scanning the genome for its binding motifs is the first step in building a
 *                 map of all active promoters.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

FIND MOTIF "TATAAT" STRAND POSITIVE AS minus10_box;

FIND MOTIF "TTGACA" STRAND POSITIVE AS minus35_box;

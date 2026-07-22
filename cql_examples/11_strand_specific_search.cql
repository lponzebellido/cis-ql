/*
 * 11_strand_specific_search.cql
 *
 * What it does: Searches for the -10 box consensus (TATAAT) separately on the positive
 *               and negative strands, producing two independent result tracks.
 *
 * Expected output: One track of TATAAT hits on the forward strand and another on the
 *                  reverse strand. The reverse-strand hits mark promoters for genes
 *                  transcribed in the opposite direction.
 *
 * Why it matters: E. coli genes are encoded on both strands. A promoter for a forward-
 *                 strand gene points right; one for a reverse-strand gene points left.
 *                 Separating them by strand is essential for correctly assigning a
 *                 promoter motif to the gene it controls.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

FIND MOTIF "TATAAT" STRAND POSITIVE AS fwd_minus10;

FIND MOTIF "TATAAT" STRAND NEGATIVE AS rev_minus10;

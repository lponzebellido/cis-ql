/*
 * 17_operon_regulatory_elements.cql
 *
 * What it does: Reconstructs the complete upstream regulatory region of the threonine
 *               biosynthesis operon (thr operon) in E. coli by searching for its key
 *               promoter elements: the -35 box, the -10 box, and the Shine-Dalgarno
 *               sequence upstream of each CDS in the operon.
 *               All three element types are then merged into a single visualization track.
 *
 * Expected output: A handful of -35, -10, and ribosome-binding-site hits upstream of
 *                  the thr operon genes (thrL, thrA, thrB, thrC), plus all sigma-70
 *                  elements upstream of every CDS unified in one track.
 *
 * Why it matters: An operon is a cluster of genes under a single promoter. Mapping the
 *                 exact positions of each regulatory element relative to the start of
 *                 translation (the Shine-Dalgarno) lets researchers validate whether an
 *                 operon's structure matches the canonical bacterial model.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

FIND MOTIF "TTGACA" WITHIN 150 BP UPSTREAM FROM GENE "thrA" AS thrA_minus35;

FIND MOTIF "TATAAT" WITHIN 50 BP UPSTREAM FROM GENE "thrA" AS thrA_minus10;

FIND MOTIF "AGGAG" WITHIN 20 BP UPSTREAM FROM GENE "thrA" AS thrA_rbs;

UNION thrA_minus35 AND thrA_minus10;

FIND MOTIF "TTGACA" WITHIN 200 BP UPSTREAM FROM CDS AS all_minus35;
FIND MOTIF "TATAAT" WITHIN 50 BP UPSTREAM FROM CDS AS all_minus10;
FIND MOTIF "AGGAG" WITHIN 20 BP UPSTREAM FROM CDS AS all_rbs;

UNION all_minus35 AND all_minus10;

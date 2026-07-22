/*
 * 15_promoter_architecture.cql
 *
 * What it does: Reconstructs the two-element architecture of sigma-70 promoters.
 *               Searches for -35 boxes (TTGACA) in the 200 bp window upstream of CDS,
 *               and -10 boxes (TATAAT) in the 50 bp window upstream of CDS.
 *               Then merges both into a single "sigma-70 elements" track, and separately
 *               computes which -35 boxes also overlap with an SP1-like GC-box motif
 *               (a sign of an unusually structured promoter).
 *
 * Expected output: ~196 -35 hits and ~104 Shine-Dalgarno-level -10 hits.
 *                  A merged view of all sigma-70 elements, plus any SP1-like sites that
 *                  co-occur upstream of CDS.
 *
 * Why it matters: Bacterial promoters have a conserved two-element structure spaced ~17 bp
 *                 apart. Being able to query both elements and visualize them together is
 *                 essential for reconstructing the regulatory architecture of any operon.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;
LOAD MATRIX "matrices/MA0079.5_SP1.pwm" AS sp1;

FIND MOTIF "TTGACA" WITHIN 200 BP UPSTREAM FROM CDS AS minus35;

FIND MOTIF "TATAAT" WITHIN 50 BP UPSTREAM FROM CDS AS minus10;

UNION minus35 AND minus10;

SCAN sp1 THRESHOLD 85 % AS sp1_sites;

INTERSECT sp1_sites AND PROMOTER;

/*
 * 10_pwm_plus_set_ops.cql
 *
 * What it does: Scans for SP1 and TBP binding sites at 85% confidence, then uses
 *               set subtraction to isolate only the sites that fall outside of any
 *               coding sequence. Both filtered results are produced as separate tracks.
 *
 * Expected output: Two tracks of intergenic binding sites — one for SP1, one for TBP —
 *                  with all in-CDS hits removed.
 *
 * Why it matters: Combining PWM scanning with set subtraction is a core bioinformatics
 *                 workflow. It models the real question: "Where could this transcription
 *                 factor bind that is not inside a gene?" The result is a clean, filtered
 *                 set of regulatory candidates ready for experimental validation.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;
LOAD MATRIX "matrices/MA0079.5_SP1.pwm" AS sp1;
LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tbp;

SCAN sp1 THRESHOLD 85 % AS sp1_sites;
SCAN tbp THRESHOLD 85 % AS tbp_sites;

EXCEPT sp1_sites FROM CDS;

EXCEPT tbp_sites FROM CDS;

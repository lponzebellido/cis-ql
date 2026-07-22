/*
 * 05_pwm_tata_scan.cql
 *
 * What it does: Scans the E. coli genome using the JASPAR TBP (TATA-binding protein)
 *               Position Weight Matrix at 80% score threshold.
 *               A second scan at 90% threshold retrieves only the strongest hits.
 *
 * Expected output: Tens of thousands of 12 bp sites at 80%, and a much smaller subset
 *                  at 90% — high-confidence TATA-like elements.
 *
 * Why it matters: Unlike literal string search, a PWM scores every possible position
 *                 against a probability model derived from known binding sites.
 *                 This recovers degenerate sites that exact matching would miss, and
 *                 the threshold controls the trade-off between sensitivity and specificity.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tbp;

SCAN tbp THRESHOLD 80 % AS tbp_relaxed;

SCAN tbp THRESHOLD 90 % AS tbp_strict;

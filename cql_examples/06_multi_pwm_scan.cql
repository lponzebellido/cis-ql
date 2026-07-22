/*
 * 06_multi_pwm_scan.cql
 *
 * What it does: Scans the genome simultaneously for two transcription factors with
 *               different roles: SP1 (GC-box binding activator) and CTCF (genome
 *               organizer / insulator). Both are scanned at 85% confidence.
 *
 * Expected output: Two separate tracks — one for SP1 candidate sites, one for CTCF
 *                  candidate sites — covering thousands of positions across the genome.
 *
 * Why it matters: Profiling multiple transcription factors in parallel is the foundation
 *                 of ChIP-seq analysis and regulatory element annotation. Seeing both
 *                 tracks at once reveals whether they share genomic locations or occupy
 *                 distinct regions.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD MATRIX "matrices/MA0079.5_SP1.pwm" AS sp1;
LOAD MATRIX "matrices/MA0139.1_CTCF.pwm" AS ctcf;

SCAN sp1 THRESHOLD 85 % AS sp1_sites;

SCAN ctcf THRESHOLD 85 % AS ctcf_sites;

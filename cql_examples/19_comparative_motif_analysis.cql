/*
 * 19_comparative_motif_analysis.cql
 * 
 * Purpose: Identifies potential cis-regulatory modules (CRMs) by looking for co-occurring motifs of different transcription factors (TBP and SP1).
 * Expectation: A unified track showing clustered SP1 and TBP binding sites across the genome, alongside GC-rich windows.
 * Why: Transcription factors rarely act alone in complex organisms. By scanning for multiple factors and uniting their tracks, researchers can identify CRMs where multiple regulatory proteins cluster together to drive gene expression.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tbp_matrix;
LOAD MATRIX "matrices/MA0079.5_SP1.pwm" AS sp1_matrix;

SCAN tbp_matrix THRESHOLD 85 % AS tbp_sites;
SCAN sp1_matrix THRESHOLD 85 % AS sp1_sites;

ANALYZE GC_CONTENT WINDOW 200 BP AS gc_windows;

UNION tbp_sites AND sp1_sites;

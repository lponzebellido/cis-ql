/*
 * 16_gc_sp1_binding.cql
 * 
 * Purpose: Combines GC content profiling with the discovery of SP1 transcription factor binding sites.
 * Expectation: A macroscopic GC content landscape (1 KB windows) alongside putative SP1 binding sites longer than 10 bp.
 * Why: SP1 is a transcription factor that typically binds to GC-rich regions (GC-boxes). Profiling both simultaneously allows researchers to visually correlate GC spikes with TF binding events.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD MATRIX "matrices/MA0079.5_SP1.pwm" AS sp1_matrix;

ANALYZE GC_CONTENT WINDOW 1 KB AS large_scale_gc;

SCAN sp1_matrix THRESHOLD 80 % AS putative_sp1_sites;

EXTRACT putative_sp1_sites WHERE LENGTH >= 10 BP;

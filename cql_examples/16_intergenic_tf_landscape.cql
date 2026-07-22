/*
 * 16_intergenic_tf_landscape.cql
 *
 * What it does: Builds a clean map of transcription factor binding candidates in
 *               intergenic space. Scans for SP1 and TBP sites, then removes all hits
 *               that land inside coding sequences (CDS). Adds a GC content profile
 *               to visualize the sequence composition context of the surviving sites.
 *
 * Expected output: Thousands of TBP and SP1 sites with CDS hits subtracted, plus a
 *                  500 bp sliding-window GC profile to contextualize them.
 *
 * Why it matters: In regulatory genomics, the most important binding sites are those
 *                 in intergenic regions — not inside genes. This workflow is a standard
 *                 first-pass filter to focus on functionally relevant regulatory candidates
 *                 while also providing the GC landscape to visually confirm that the
 *                 surviving sites cluster in GC-rich regions, as expected for SP1.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;
LOAD MATRIX "matrices/MA0079.5_SP1.pwm" AS sp1;
LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tbp;

ANALYZE GC_CONTENT WINDOW 500 BP AS gc_landscape;

SCAN sp1 THRESHOLD 85 % AS sp1_all;
SCAN tbp THRESHOLD 85 % AS tbp_all;

EXCEPT sp1_all FROM CDS;

EXCEPT tbp_all FROM CDS;

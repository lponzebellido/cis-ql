/*
 * 20_mobile_genetic_elements.cql
 *
 * What it does: Identifies signatures of mobile genetic elements (transposons, insertion
 *               sequences) in the E. coli genome. Searches for terminal inverted repeats
 *               (TIR) — short palindromic sequences at the boundaries of mobile elements —
 *               and compares them with annotated regions to filter for those outside genes.
 *
 *               Additionally, computes a GC content profile: mobile elements are often
 *               AT-rich compared to the host genome, creating dips in the GC landscape.
 *
 * Expected output: Hits for the IS1 TIR sequence (TGGTCGGTGCCTTAAGC variant), a GC
 *                  content profile at 200 bp resolution, and all annotated CDS subtracted
 *                  to leave only intergenic hits — the expected location of TIRs.
 *
 * Why it matters: Insertion sequences are the most abundant mobile elements in E. coli
 *                 and have a major impact on genome plasticity and gene disruption.
 *                 Mapping their terminal repeats allows researchers to identify all
 *                 copies of a given IS element in a genome, even if they were not
 *                 annotated in the GFF file.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

ANALYZE GC_CONTENT WINDOW 200 BP AS gc_is_context;

FIND MOTIF "TGGTCG" STRAND POSITIVE AS tir_fwd;
FIND MOTIF "TGGTCG" STRAND NEGATIVE AS tir_rev;

UNION tir_fwd AND tir_rev;

EXCEPT tir_fwd FROM CDS;

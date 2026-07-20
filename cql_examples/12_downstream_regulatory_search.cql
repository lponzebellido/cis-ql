/*
 * 12_downstream_regulatory_search.cql
 * 
 * Purpose: Finds a regulatory motif downstream of a gene to identify transcription terminators.
 * Expectation: Rho-independent terminator candidates (GC-rich hairpins followed by poly-T runs) within 200 bp after each CDS.
 * Why: While promoters are upstream, terminators are downstream. Spatial constraints allow us to specifically target the regions where terminators are biologically functional.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

FIND MOTIF "G{3,6}[AT]{3,6}T{4,8}" WITHIN 200 BP DOWNSTREAM FROM CDS AS terminators;

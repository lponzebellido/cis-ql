/*
 * 21_riboswitch_discovery.cql
 * 
 * Purpose: Discovers potential riboswitches or attenuators by looking for highly conserved RNA structural motifs in the 5' UTR of specific metabolic genes.
 * Expectation: A set of conserved structural motifs located strictly in the untranslated regions preceding genes.
 * Why: Riboswitches are non-coding RNA elements that regulate gene expression. Searching within the 5' UTR using pattern matching allows researchers to discover elements that regulate translation based on metabolite binding.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

// Riboswitches often feature specific hairpin loops. We use a regex to approximate this.
FIND MOTIF "G{3,6}.{4,8}C{3,6}" WITHIN 150 BP UPSTREAM FROM GENE AS putative_riboswitches;

EXTRACT putative_riboswitches WHERE LENGTH >= 20 BP;

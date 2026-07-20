/*
 * 20_crispr_array_discovery.cql
 * 
 * Purpose: Identifies putative CRISPR arrays by searching for clustered, short repeating sequences interspersed with unique spacers.
 * Expectation: Regions showing dense clustering of a specific repeat sequence (like the 29bp CRISPR repeat in E. coli).
 * Why: CRISPR-Cas systems are adaptive immune systems in bacteria. Finding repeating patterns near genes like 'cas9' allows researchers to identify novel CRISPR arrays and their associated spacers de novo.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

FIND MOTIF "CGGTTTATCCCCGCTGGCGCGGGGAACTC" STRAND POSITIVE AS crispr_repeats;

EXTRACT crispr_repeats;

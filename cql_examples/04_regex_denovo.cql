/*
 * 04_regex_denovo.cql
 * 
 * Purpose: Discovers Open Reading Frames (ORFs) de novo using regular expressions, without relying on annotations.
 * Expectation: Putative start-to-stop ORFs on the positive strand, filtered to be at least 15 bp long.
 * Why: Regular expressions allow for flexible pattern matching. This example simulates finding novel genes by looking for an ATG start codon followed by a valid stop codon downstream of a putative promoter.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

FIND MOTIF "TATA[AT]A[AT]" STRAND POSITIVE AS putative_promoters;

FIND MOTIF "ATG(...)*?(TAA|TAG|TGA)"
  WITHIN 200 BP DOWNSTREAM FROM putative_promoters
  AS putative_orfs
  WHERE LENGTH >= 15 BP;

EXTRACT putative_orfs;

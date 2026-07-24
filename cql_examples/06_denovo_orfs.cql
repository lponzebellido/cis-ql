/*
 * 06_denovo_orfs.cql
 *
 * This script performs de novo discovery of open reading frames (ORFs) on raw FASTA data.
 * First, putative promoter elements are identified using a regular expression pattern. Next, 
 * a downstream search locates potential coding regions bounded by start (ATG) and stop codons 
 * (TAA, TAG, TGA) within 2500 BP downstream of those promoters, filtered to a minimum length of 600 BP.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS sample;

FIND MOTIF "TATA[AT]A[AT]" STRAND POSITIVE AS putative_promoters;

FIND MOTIF "ATG(...)*?(TAA|TAG|TGA)"
    WITHIN 2500 BP DOWNSTREAM FROM putative_promoters
    AS candidate_orfs
    WHERE LENGTH > 600 BP;

EXTRACT candidate_orfs;

// 04_regex_denovo.cql
// De novo ORF discovery using regular expressions on a raw sequence with no annotation.
// Expect: putative start-to-stop ORFs on the positive strand, filtered to >= 15 bp.

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

// Find consensus promoter-like patterns using IUPAC ambiguity codes
FIND MOTIF "TATA[AT]A[AT]" STRAND POSITIVE AS putative_promoters;

// Find short ORFs (ATG...stop) in the 200 bp window downstream of each candidate promoter
FIND MOTIF "ATG(...)*?(TAA|TAG|TGA)"
  WITHIN 200 BP DOWNSTREAM FROM putative_promoters
  AS putative_orfs
  WHERE LENGTH >= 15 BP;

EXTRACT putative_orfs;

// regex_test.cql
// De novo ORF prediction chained with spatial search using regular expressions.
// Expect: short putative ORFs (>= 15 bp) located within 200 bp of candidate promoter patterns.

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

// Detect TATA-like promoter consensus sequences using an IUPAC regex
FIND MOTIF "TATA[AT]A[AT]" STRAND POSITIVE AS putative_promoters;

// Look for start-to-stop ORFs immediately downstream of each promoter hit
FIND MOTIF "ATG(...)*?(TAA|TAG|TGA)"
  WITHIN 200 BP DOWNSTREAM FROM putative_promoters
  AS putative_genes
  WHERE LENGTH >= 15 BP;

EXTRACT putative_genes;

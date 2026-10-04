LOAD SEQUENCE "data_examples/regulatory_demo.fasta" AS genome;
LOAD ANNOTATION "data_examples/transcripts_demo.gff3" AS annotation;

EXTRACT GENE AS gene
  WHERE ID = "gene_1";

EXTRACT TRANSCRIPTS OF gene
  TYPES ["mRNA", "transcript"]
  SELECT ALL
  AS transcripts;

EXTRACT TRANSCRIPTS OF gene
  TYPES ["mRNA", "transcript"]
  SELECT ID "tx2"
  AS tx2;

EXTRACT TRANSCRIPTS OF gene
  TYPES ["mRNA", "transcript"]
  SELECT ATTRIBUTE "tag" = "canonical"
  AS canonical_transcript;

DEFINE PROMOTERS OF transcripts
  FROM TSS
  UPSTREAM 50 BP
  DOWNSTREAM 10 BP
  AS transcript_promoters;

EXPORT transcripts TO "transcripts.tsv" FORMAT TSV;
EXPORT canonical_transcript TO "canonical_transcript.gff3" FORMAT GFF3;
EXPORT transcript_promoters TO "transcript_promoters.bed" FORMAT BED;

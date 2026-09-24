LOAD SEQUENCE "data_examples/ecoli2.fna" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annotation;

EXTRACT FEATURE AS thra_gene
  WHERE TYPE = "gene"
    AND ID = "gene-b0002"
    AND NAME = "thrA"
    AND ATTRIBUTE "locus_tag" = "b0002";

EXTRACT FEATURE AS thra_coding_sequence
  WHERE TYPE = "CDS"
    AND PARENT = "gene-b0002"
    AND PHASE = "0"
    AND ATTRIBUTE "protein_id" = "AAC73113.1";

EXPORT thra_gene TO "thra_gene.tsv" FORMAT TSV;
EXPORT thra_coding_sequence TO "thra_coding_sequence.gff3" FORMAT GFF3;

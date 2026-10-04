LOAD ANNOTATION "data_examples/genomic.gff" AS annotation;

EXTRACT FEATURE AS coding_records
  WHERE TYPE = "CDS"
    AND STRAND = "+"
    AND PHASE = "0";

EXTRACT coding_records AS thra_cds
  WHERE ATTRIBUTE "protein_id" = "AAC73113.1";

EXPORT thra_cds TO "filters.tsv" FORMAT TSV;

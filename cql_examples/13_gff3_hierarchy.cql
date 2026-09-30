LOAD SEQUENCE "data_examples/ecoli2.fna" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annotation;

VALIDATE ANNOTATION annotation AS hierarchy_report;

EXTRACT FEATURE AS thra_gene
  WHERE TYPE = "gene"
    AND ID = "gene-b0002"
    AND NAME = "thrA"
    AND ATTRIBUTE "locus_tag" = "b0002";

EXTRACT CHILDREN OF thra_gene AS thra_direct_features;

EXTRACT DESCENDANTS OF thra_gene AS thra_coding_sequence
  WHERE TYPE = "CDS"
    AND PHASE = "0"
    AND ATTRIBUTE "protein_id" = "AAC73113.1";

EXTRACT PARENTS OF thra_coding_sequence AS thra_direct_parents;

EXTRACT ANCESTORS OF thra_coding_sequence AS thra_ancestor_genes
  WHERE TYPE = "gene";

EXTRACT FEATURE AS multipart_cds_parts
  WHERE ID = "cds-gnl|b4623|CDS4514";

GROUP multipart_cds_parts BY ID AS multipart_cds_features;
EXTRACT MEMBERS OF multipart_cds_features AS recovered_cds_parts;

EXPORT thra_gene TO "thra_gene.tsv" FORMAT TSV;
EXPORT thra_direct_features TO "thra_direct_features.tsv" FORMAT TSV;
EXPORT thra_coding_sequence TO "thra_coding_sequence.gff3" FORMAT GFF3;
EXPORT thra_ancestor_genes TO "thra_ancestor_genes.tsv" FORMAT TSV;
EXPORT hierarchy_report TO "hierarchy_report.tsv" FORMAT TSV;
EXPORT multipart_cds_features TO "multipart_cds_features.tsv" FORMAT TSV;
EXPORT recovered_cds_parts TO "recovered_cds_parts.gff3" FORMAT GFF3;

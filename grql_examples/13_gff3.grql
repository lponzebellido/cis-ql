LOAD SEQUENCE "data_examples/ecoli2.fna" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annotation;

VALIDATE ANNOTATION annotation AS report;

EXTRACT FEATURE AS thra
  WHERE TYPE = "gene"
    AND ID = "gene-b0002"
    AND NAME = "thrA"
    AND ATTRIBUTE "locus_tag" = "b0002";

EXTRACT CHILDREN OF thra AS children;

EXTRACT DESCENDANTS OF thra AS cds
  WHERE TYPE = "CDS"
    AND PHASE = "0"
    AND ATTRIBUTE "protein_id" = "AAC73113.1";

EXTRACT PARENTS OF cds AS parents;

EXTRACT ANCESTORS OF cds AS genes
  WHERE TYPE = "gene";

EXTRACT FEATURE AS cds_parts
  WHERE ID = "cds-gnl|b4623|CDS4514";

GROUP cds_parts BY ID AS cds_groups;
EXTRACT MEMBERS OF cds_groups AS restored_parts;

EXPORT thra TO "thra.tsv" FORMAT TSV;
EXPORT children TO "children.tsv" FORMAT TSV;
EXPORT cds TO "cds.gff3" FORMAT GFF3;
EXPORT genes TO "genes.tsv" FORMAT TSV;
EXPORT report TO "gff3_report.tsv" FORMAT TSV;
EXPORT cds_groups TO "cds_groups.tsv" FORMAT TSV;
EXPORT restored_parts TO "cds_parts.gff3" FORMAT GFF3;

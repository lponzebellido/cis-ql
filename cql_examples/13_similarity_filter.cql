// 13_length_filters.cql
// Filters annotated genomic features by length to isolate regions in specific size ranges.
// Expect: long promoters, short CDSs, and annotated exons.

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

// Long promoter regions (> 150 bp)
EXTRACT PROMOTER WHERE LENGTH > 150 BP;

// Short coding sequences — candidates for small regulatory proteins
EXTRACT CDS WHERE LENGTH < 300 BP;

// All annotated exons
EXTRACT EXON;

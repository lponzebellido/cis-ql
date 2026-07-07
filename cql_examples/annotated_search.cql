// annotated_search.cql
// Classic annotation-driven workflow: spatial motif search relative to a known gene,
// followed by feature extraction and set intersection using the GFF3 annotation.
// Expect: TATA-like hits near thrA, filtered promoters, and gene-CDS overlaps.

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

// Find the Pribnow box within 500 bp upstream of gene thrA
FIND MOTIF "TATAAT" WITHIN 500 BP UPSTREAM FROM GENE "thrA" AS thrA_tata;

// Extract all promoters between 100 and 500 bp
EXTRACT PROMOTER WHERE LENGTH > 100 BP AND LENGTH < 500 BP;

// Find where annotated genes and CDS overlap at high similarity
INTERSECT GENE AND CDS WHERE SIMILARITY >= 90.5%;

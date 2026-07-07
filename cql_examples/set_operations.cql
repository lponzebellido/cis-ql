// set_operations.cql
// Demonstrates INTERSECT, UNION, and EXCEPT on motif search results.
// Expect: overlapping Pribnow/conserved hits, a merged track, and a filtered set minus exons.

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

// Sigma-70 -10 box candidates
FIND MOTIF "TATAAT" STRAND POSITIVE AS pribnow;

// Conserved bacterial regulatory element (riboswitch-like)
FIND MOTIF "GATAAG" STRAND POSITIVE AS conserved;

// INTERSECT: positions where both motifs co-occur (potential dual-regulated regions)
INTERSECT pribnow AND conserved WHERE LENGTH > 5 BP;

// UNION: all hits from either motif in a single combined track
UNION pribnow AND conserved;

// EXCEPT: remove Pribnow hits that fall inside annotated exons (likely false positives)
EXCEPT pribnow FROM EXON;

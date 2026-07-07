// 12_downstream_regulatory_search.cql
// Finds a regulatory motif downstream of a gene — useful for finding transcription terminators.
// Expect: rho-independent terminator candidates (GC-rich stems) within 200 bp after each CDS.

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;

// Rho-independent terminators have a GC-rich hairpin followed by a poly-T run
FIND MOTIF "G{3,6}[AT]{3,6}T{4,8}" WITHIN 200 BP DOWNSTREAM FROM CDS AS terminators;

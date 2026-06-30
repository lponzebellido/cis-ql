LOAD SEQUENCE "data_examples/ecoli.fasta" AS ecoli_genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS ecoli_annot;

// Encuentra motivos reguladores
FIND MOTIF "TATAAT" WITHIN 500 BP UPSTREAM FROM GENE AS promoters;

// Encuentra regiones conservadas
FIND MOTIF "GATAAG" STRAND POSITIVE AS conserved_regions;

// Intersección: promotores que además están en una región conservada
INTERSECT promoters AND conserved_regions WHERE LENGTH > 10 BP;

// Diferencia: promotores que NO están en un exón (para filtrar falsos positivos)
EXCEPT promoters FROM EXON WHERE SIMILARITY >= 90.0%;

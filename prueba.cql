// Archivo de prueba para Cis-QL
LOAD SEQUENCE "ecoli.fasta" AS ecoli_genome;
FIND MOTIF "MYB" WITHIN 2000 BP UPSTREAM FROM GENE "X";
EXTRACT PROMOTER WHERE LENGTH > 500 BP;

/* Comentariooooos */
INTERSECT EXON AND INTRON WHERE SIMILARITY >= 90.5%;
FIND MOTIF "X" STRAND POSITIVE CHR "1";

@ // Esto deberia generar un TKN_ERROR pero continuar

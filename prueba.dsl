// Archivo de prueba para Cis-QL
LOAD SEQUENCE "genome.fasta" AS genome_data;
FIND MOTIF "MYB" WITHIN 2000 BP UPSTREAM FROM GENE "X";
EXTRACT PROMOTER WHERE LENGTH > 500 BP;

/* Probando las nuevas entidades y operaciones */
INTERSECT EXON AND INTRON WHERE SIMILARITY >= 90.5%;
FIND MOTIF "X" STRAND POSITIVE CHR "1";

@ // Esto deberia generar un TKN_ERROR pero continuar

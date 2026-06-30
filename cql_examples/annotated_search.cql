// Carga del genoma y su archivo de anotaciones GFF3
LOAD SEQUENCE "data_examples/ecoli.fasta" AS ecoli_genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS ecoli_annot;

// Búsqueda espacial relativa: Buscar la caja TATA exactamente antes del gen "thrA"
FIND MOTIF "TATAAT" WITHIN 20000 BP UPSTREAM FROM GENE "thrA";

// Extracción y filtrado utilizando metadatos del GFF3
EXTRACT PROMOTER WHERE LENGTH > 100 BP AND LENGTH < 500 BP;

// Operación lógica apoyada en las coordenadas geográficas del genoma
INTERSECT GENE AND CDS WHERE SIMILARITY >= 90.5%;

// Carga de la secuencia cruda, sin GFF3
LOAD SEQUENCE "data_examples/ecoli.fasta" AS raw_seq;

// Búsqueda directa del elemento cis-regulador en todo el cromosoma
FIND MOTIF "GATAAG" STRAND POSITIVE AS gata_factors;

// Filtrar
EXTRACT gata_factors WHERE SIMILARITY >= 85.0%;

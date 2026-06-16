// De Novo Discovery Test
LOAD SEQUENCE "examples/ecoli.fasta" AS ecoli;
FIND MOTIF "ATG" STRAND NEGATIVE AS Met;
FIND MOTIF "ATGCGA" AS binding_sites;

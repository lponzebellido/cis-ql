LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

ANALYZE GC_CONTENT WINDOW 1 KB AS gc_profile;

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tbp;
LOAD MATRIX "matrices/MA0079.5_SP1.pwm" AS sp1;

IF GC_CONTENT > 50 % THEN
    SCAN sp1 THRESHOLD 80 % AS gc_sites;
ELSE
    SCAN tbp THRESHOLD 80 % AS at_sites;
ENDIF;

// pwm_scan.cql
// Probabilistic TATA-box discovery using the TBP PWM (JASPAR MA0108.1), then
// extracted and filtered to only high-confidence sites.
// Expect: TBP binding site predictions above 80% of the maximum PWM score.

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

// Load TBP matrix (TATA-Binding Protein — core promoter element)
LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tata_pwm;

// Scan positive strand only — high threshold for specificity
SCAN tata_pwm STRAND POSITIVE THRESHOLD 80% AS tata_sites;

// Pull out the result for visualization
EXTRACT tata_sites;

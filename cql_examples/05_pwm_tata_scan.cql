// 05_pwm_tata_scan.cql
// Probabilistic TATA-box discovery using the TBP Position Weight Matrix (JASPAR MA0108.1).
// Expect: high-confidence TATA-box binding sites scoring >= 80% of the max PWM score.

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

// Load the TBP (TATA-Binding Protein) matrix from the JASPAR database
LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tata_pwm;

// Scan positive strand only; keep only sites above 80% of theoretical max score
SCAN tata_pwm STRAND POSITIVE THRESHOLD 80% AS tata_sites;

// Extract the results for visualization
EXTRACT tata_sites;

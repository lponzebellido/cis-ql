// 10_pwm_plus_set_ops.cql
// Combines PWM scanning with set intersection to find regions where two TFs co-occur.
// Expect: positions where a TATA box and a GC box are predicted within the same area.

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tata_pwm;
LOAD MATRIX "matrices/MA0079.5_SP1.pwm" AS sp1_pwm;

// High-confidence scans on both strands
SCAN tata_pwm THRESHOLD 80% AS tata_sites;
SCAN sp1_pwm  THRESHOLD 75% AS gc_boxes;

// Intersect: regions where both TF binding sites co-occur — characteristic of core promoters
INTERSECT tata_sites AND gc_boxes;

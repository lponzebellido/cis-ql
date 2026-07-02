/* Ejemplo: Búsqueda de elementos cis-reguladores con PWM (Position Weight Matrix) */

// 1. Cargar la secuencia genómica
LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

// 2. Cargar matrices PWM de factores de transcripción (formato JASPAR)
LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tata_pwm;

// 3. Escanear: identificar sitios de unión a TBP (TATA box) con score >= 80%
SCAN tata_pwm STRAND POSITIVE THRESHOLD 80% AS tata_sites;

// 4. Extraer los resultados finales
EXTRACT tata_sites;

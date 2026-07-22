/*
 * 21_comprehensive_regulatory_map.cql
 *
 * What it does: Builds a comprehensive view of the regulatory landscape of E. coli K-12
 *               in a single script. It combines:
 *                 (1) Genome-wide GC content at 1 KB resolution.
 *                 (2) CpG island detection for global structural context.
 *                 (3) Sigma-70 -35 and -10 box mapping upstream of all CDS.
 *                 (4) Ribosome binding site (Shine-Dalgarno) mapping.
 *                 (5) TBP and SP1 PWM scans filtered to intergenic regions.
 *                 (6) All large genes (> 2 KB) as a structural reference track.
 *
 * Expected output: Six tracks covering all major regulatory signal types.
 *                  This is the closest Cis-QL can get to a full regulatory map
 *                  from a raw genome and annotation file alone.
 *
 * Why it matters: In practice, building a complete regulatory map of a bacterium requires
 *                 integrating multiple evidence types simultaneously. This script shows how
 *                 Cis-QL can replace a multi-step pipeline of shell scripts and Python code
 *                 with a single, readable query that produces all tracks at once.
 */

LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;
LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tbp;
LOAD MATRIX "matrices/MA0079.5_SP1.pwm" AS sp1;

ANALYZE GC_CONTENT WINDOW 1 KB AS gc_1kb;

ANALYZE CPG_ISLANDS AS cpg_islands;

FIND MOTIF "TTGACA" WITHIN 200 BP UPSTREAM FROM CDS AS sigma70_minus35;
FIND MOTIF "TATAAT" WITHIN 50 BP UPSTREAM FROM CDS AS sigma70_minus10;
FIND MOTIF "AGGAG" WITHIN 20 BP UPSTREAM FROM CDS AS ribosome_binding_sites;

UNION sigma70_minus35 AND sigma70_minus10;

SCAN tbp THRESHOLD 85 % AS tbp_sites;
SCAN sp1 THRESHOLD 85 % AS sp1_sites;

EXCEPT tbp_sites FROM CDS;
EXCEPT sp1_sites FROM CDS;

EXTRACT GENE WHERE LENGTH > 2000 BP;

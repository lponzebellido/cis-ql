# Cis-QL: Cis-Regulatory Query Language & Studio

**Cis-QL** is a specialized Domain-Specific Language (DSL) and desktop environment engineered for bioinformatics, computational genomics, and cis-regulatory element discovery. Built on a native C++11 core, it enables researchers to search, filter, extract, and manipulate genomic regions, regulatory motifs, and sequence features using a human-readable, declarative syntax.

---

## Overview & Purpose

Identifying cis-regulatory elements (such as promoters, enhancers, and transcription factor binding sites) within large genomic datasets traditionally requires stringing together multiple command-line utilities (e.g., `bedtools`, `grep`, `awk`) or writing custom scripts in Python or Perl. These imperative approaches often result in complex, hard-to-maintain pipelines that are difficult to reproduce.

**Cis-QL provides a unified, declarative interface.** Instead of writing imperative code to parse files, perform string matching, compute spatial distances, and calculate set intersections, researchers declare the desired biological criteria. The underlying C++ execution engine handles parsing, spatial indexing, multi-threaded alignment, and interval arithmetic automatically.

---

## What Problems Does Cis-QL Solve?

- **Transcription Factor Binding Site Discovery:** Integrates Position Weight Matrices (PWMs) from standard databases (such as JASPAR) to identify TF binding sites using probabilistic log-odds scoring rather than rigid exact matching.
- **Metagenomic Bioprospecting:** Enables *de novo* motif and open reading frame (ORF) discovery on raw, unannotated FASTA contigs from environmental samples lacking curated GFF3 annotations.
- **Spatial Relative Queries:** Simplifies relative proximity searches (e.g., locating specific consensus motifs within designated base-pair windows upstream or downstream of coding sequences).
- **Pipeline Unification & Reproducibility:** Condenses multi-step bioinformatic shell workflows into concise, shareable, and self-documenting query scripts.

---

## Primary Operational Modes

Cis-QL operates across three primary modes:

1. **Annotation-Driven Analysis:** Queries established `GFF3` annotation files to analyze known genes, exons, CDS, and regulatory features.
2. **De Novo Discovery Mode:** Constructs a virtual annotation layer in memory directly from raw `.fasta` sequence data using literal strings, IUPAC ambiguity codes, or Regular Expressions.
3. **Probabilistic PWM Scanning:** Scans sequences using Position Weight Matrices with log-odds PSSM scoring, evaluating hit thresholds as percentages of maximum theoretical scores.

---

## Key Technical Features

- **Native IUPAC Degeneration Engine:** Translates IUPAC nucleotide ambiguity codes (`R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`, `H`, `V`, `N`) automatically into regular expression search patterns (e.g., `TATAWAW` translates to `TATA[AT]A[AT]`).
- **High-Performance Multi-Chromosome Parallelism:** Leverages `std::async` to execute sequence scanning across multiple contigs and chromosomes in parallel.
- **Sweep-Line Interval Algebra:** Evaluates set operations (`INTERSECT`, `UNION`, `EXCEPT`) on genomic intervals using an $O(N \log N)$ two-pointer sweep-line algorithm.
- **Pairwise Smith-Waterman Local Alignment:** Performs dynamic programming local sequence alignment in native C++ to filter candidate features by sequence similarity (`WHERE SIMILARITY > 70 %`).
- **Control Flow & Scripting (v2.0):** Supports conditional execution (`IF / ELSE`) based on sequence metrics and batch iteration (`FOREACH`) over matrix collections.
- **Cis-QL Studio (GUI):** Desktop workspace built with Electron and React 18, featuring a visual editor, interactive track canvas and live JSON data export.

---

## Compilation and Execution

### Building the C++ Binary

Cis-QL requires a standard C++11 compiler and `make`:

```bash
make clean
make
```

### Running Queries via Command Line

Execute a `.cql` script using the `cisql` binary:

```bash
./cisql cql_examples/01_extract_genes.cql
```

Use the `--debug` flag to inspect compilation phases, including token stream, Abstract Syntax Tree (AST), Symbol Table, Intermediate Representation (IR), and execution steps:

```bash
./cisql cql_examples/14_integrated_query.cql --debug
```

### Launching Cis-QL Studio

To run the interactive desktop graphical environment:

```bash
cd cisql-studio
npm start
```

---

## Language Reference

### 1. Loading Data (`LOAD`)

Load sequence files (FASTA), annotation files (GFF3), and matrix files (JASPAR format):

```sql
LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;
LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tbp_matrix;
```

### 2. Motif Searching & Spatial Conditions (`FIND MOTIF`)

Locate exact motifs, regular expressions, or IUPAC degenerate strings, with optional spatial constraints relative to other features:

```sql
-- IUPAC degenerate motif search on the positive strand
FIND MOTIF "TATAWAW" STRAND POSITIVE AS tata_boxes;

-- Spatial constraint: motif within 200 BP upstream of coding sequences
FIND MOTIF "TTGACA" WITHIN 200 BP UPSTREAM FROM CDS AS minus35_promoters;

-- De novo ORF search downstream of putative promoters
FIND MOTIF "ATG(...)*?(TAA|TAG|TGA)"
    WITHIN 300 BP DOWNSTREAM FROM tata_boxes
    AS candidate_orfs
    WHERE LENGTH > 600 BP;
```

### 3. Position Weight Matrix Scanning (`SCAN`)

Scan loaded sequences using Position Weight Matrices with log-odds scoring:

```sql
SCAN tbp_matrix STRAND POSITIVE THRESHOLD 80 % AS tbp_sites;
```

### 4. Biological & Structural Analysis (`ANALYZE`)

Calculate sliding-window GC content profiles or identify CpG islands:

```sql
ANALYZE GC_CONTENT WINDOW 1 KB AS gc_profile;
ANALYZE CPG_ISLANDS AS cpg_islands;
```

### 5. Set Operations (`INTERSECT`, `UNION`, `EXCEPT`)

Combine or filter interval sets using high-speed interval algebra:

```sql
INTERSECT sp1_sites AND cpg_islands;
UNION minus35_box AND minus10_box;
EXCEPT ctcf_sites FROM CDS;
```

### 6. Feature Extraction & Filtering (`EXTRACT`, `WHERE`)

Filter genomic entities by physical length or alignment similarity:

```sql
EXTRACT GENE WHERE LENGTH >= 500 BP AND LENGTH <= 3 KB;
EXTRACT GENE WHERE LENGTH > 1 KB AND SIMILARITY > 70 %;
```

### 7. Control Flow (`IF`, `FOREACH`)

Control query execution paths and iterate over collections:

```sql
-- Conditional execution based on sequence properties
IF GC_CONTENT > 50 % THEN
    SCAN sp1 THRESHOLD 80 % AS gc_sites;
ELSE
    SCAN tbp THRESHOLD 80 % AS at_sites;
ENDIF;

-- Batch processing over matrix lists
FOREACH m IN [tbp, sp1, ctcf] DO
    SCAN m THRESHOLD 80 % AS tf_sites;
ENDFOR;
```

---

## Curated Examples Suite (`cql_examples/`)

The repository includes 16 structured `.cql` scripts demonstrating specific language capabilities:

| Script | Description | Primary Features |
| :--- | :--- | :--- |
| `01_extract_genes.cql` | Structural gene extraction | `EXTRACT`, `WHERE LENGTH` |
| `02_length_filter.cql` | Multi-condition range filtering | Logical `AND`, bounded intervals |
| `03_similarity_align.cql` | Paralog discovery via local alignment | Smith-Waterman `SIMILARITY > 70 %` |
| `04_iupac_motifs.cql` | Degenerate promoter motif search | IUPAC translation engine (`TATAWAW`) |
| `05_spatial_promoters.cql` | Upstream regulatory element search | `WITHIN 200 BP UPSTREAM FROM` |
| `06_denovo_orfs.cql` | Unannotated ORF discovery | Regex matching, virtual annotations |
| `07_pwm_scanning.cql` | JASPAR matrix scanning | Log-odds PSSM scoring (`THRESHOLD 80 %`) |
| `08_ctcf_insulators.cql` | Chromatin insulator mapping | `EXCEPT` set subtraction |
| `09_cpg_islands.cql` | Epigenetic CpG island profiling | `ANALYZE CPG_ISLANDS`, `INTERSECT` |
| `10_promoter_union.cql` | Bipartite promoter element merger | Multi-track consolidation via `UNION` |
| `11_strand_search.cql` | Sense vs. antisense motif profiling | `STRAND POSITIVE / NEGATIVE` |
| `12_gc_content.cql` | Sliding-window GC landscape | `ANALYZE GC_CONTENT WINDOW` |
| `13_complex_where.cql` | Multi-property conditional queries | Combined `LENGTH` & `SIMILARITY` |
| `14_integrated_query.cql` | Master genome-wide regulatory map | Complete multi-track pipeline |
| `15_if_else_branching.cql` | Conditional flow execution | `IF / ELSE / ENDIF` |
| `16_foreach_batch_scan.cql` | Batch processing over matrix lists | `FOREACH / DO / ENDFOR` |

---

## Cis-QL Studio Interface

Cis-QL Studio provides a graphical user interface for query development and visualization:

- **Interactive Track Canvas:** Renders ruler coordinates, GC content curves, feature annotation boxes, and nucleotide text at high zoom.
- **Data Export:** Exports execution results to `.cisql_results.json` for external downstream integration.

---

## Compiler Architecture

```
                                  [ .cql Source Query ]
                                            │
                                            ▼
                                   Lexical Analyzer (Lexer)
                                            │
                                            ▼
                                  LL(1) Recursive Parser
                                            │
                                            ▼
                                  Abstract Syntax Tree (AST)
                                            │
                                            ▼
                                    Semantic Analyzer
                                            │
                                            ▼
                               Intermediate Representation (IR)
                                            │
                                            ▼
                       Multithreaded C++11 Execution Engine
                 (FastaReader, GFFReader, Sweep-Line, PWMScanner)
                                            │
                                            ▼
                                [ .cisql_results.json ]
                                            │
                                            ▼
                              Cis-QL Studio GUI Visualizer
```

---

## Current Project Status

Cis-QL is an active academic research project and compiler design implementation.

Current operational components:
- LL(1) Lexical, syntactic, and semantic analyzers.
- Virtual annotation engine for unannotated sequence data.
- IUPAC degeneration engine and regex matching integration.
- Parallel multithreaded PSSM matrix scanner.
- Smith-Waterman pairwise alignment module.
- Sweep-line interval algebra engine (`INTERSECT`, `UNION`, `EXCEPT`).
- Control flow execution engine (`IF/ELSE` and `FOREACH`).
- Desktop GUI workspace (`Cis-QL Studio`).

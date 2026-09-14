# Cis-QL: Cis-Regulatory Query Language & Studio

**Cis-QL** is a specialized Domain-Specific Language (DSL) and desktop environment engineered for bioinformatics, computational genomics, and cis-regulatory element discovery. Built on a native C++11 core, it enables researchers to search, filter, extract, and manipulate genomic regions, regulatory motifs, and sequence features using a human-readable, declarative syntax.

---

## Overview & Purpose

Identifying cis-regulatory elements (such as promoters, enhancers, and transcription factor binding sites) often involves coordinating command-line utilities or custom scripts. Cis-QL explores whether a single declarative syntax can make a subset of these workflows shorter and easier to reproduce.

**Cis-QL provides a unified, declarative interface.** The C++ execution engine currently handles FASTA and GFF3 parsing, motif and matrix scanning, bounded spatial searches, local alignment, and interval arithmetic.

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
3. **Probabilistic PWM Scanning:** Scans sequences using Position Weight Matrices with log-odds PSSM scoring, normalizing thresholds across each matrix's attainable score range.

---

## Key Technical Features

- **Native IUPAC Degeneration Engine:** Translates IUPAC nucleotide ambiguity codes (`R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`, `H`, `V`, `N`) automatically into regular expression search patterns (e.g., `TATAWAW` translates to `TATA[AT]A[AT]`).
- **Multi-Chromosome Scanning:** Uses `std::async` to scan loaded contigs independently. Performance depends on contig count, input size, and the host implementation.
- **Interval Algebra and Spatial Selection:** Supports geometric `INTERSECT`, `UNION`, and `EXCEPT` operations, plus directional `OVERLAPS` selection that retains complete query records and their evidence.
- **Pairwise Smith-Waterman Local Alignment:** Computes a normalized local-alignment score in C++. `SIMILARITY TO alias` requires an explicit one-region reference set. The older implicit-reference form remains available for compatibility.
- **Explicit Dataset Context:** Multiple FASTA and GFF3 datasets can be loaded and selected deterministically with `USE SEQUENCE` and `USE ANNOTATION`.
- **Standard Result Export:** Named region sets can be written as BED, GFF3, or TSV; GC profiles can be written as TSV.
- **Control Flow & Scripting (v2.0):** Supports conditional execution (`IF / ELSE`) based on sequence metrics and batch iteration (`FOREACH`) over matrix collections.
- **Cis-QL Studio (GUI):** Desktop environment built with Electron and React 18 for interactive query authoring, multi-track genomic visualization, and live result inspection.

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
LOAD SEQUENCE "data_examples/ecoli2.fna" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annot;
LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tbp_matrix;
```

The example suite also includes the plant MYB profile `MA0054.1` from
[JASPAR CORE](https://jaspar.elixir.no/matrix/MA0054.1/). The accompanying
`anthocyanin_regulatory_demo` FASTA/GFF3 pair is synthetic and exists only to
make expected regulatory-query results small, deterministic, and inspectable.

When multiple sequence or annotation datasets are loaded, select the active
context explicitly:

```sql
USE SEQUENCE genome;
USE ANNOTATION annot;
```

### 2. Motif Searching & Spatial Conditions (`FIND MOTIF`)

Locate exact motifs, regular expressions, or IUPAC degenerate strings, with optional spatial constraints relative to other features:

```sql
// IUPAC degenerate motif search on the positive strand
FIND MOTIF "TATAWAW" STRAND POSITIVE AS tata_boxes;

// Spatial constraint: motif within 200 BP upstream of coding sequences
FIND MOTIF "TTGACA" WITHIN 200 BP UPSTREAM FROM CDS AS minus35_promoters;

// De novo ORF search downstream of putative promoters
FIND MOTIF "ATG(...)*?(TAA|TAG|TGA)"
    WITHIN 2500 BP DOWNSTREAM FROM tata_boxes
    AS candidate_orfs
    WHERE LENGTH > 600 BP;
```

### 3. Position Weight Matrix Scanning (`SCAN`)

Scan loaded sequences using Position Weight Matrices with log-odds scoring:

```sql
SCAN tbp_matrix STRAND POSITIVE THRESHOLD 80 % AS tbp_sites;

SCAN tbp_matrix BACKGROUND FROM genome
    QVALUE <= 0.05 AS significant_tbp_sites;
```

### 4. Biological & Structural Analysis (`ANALYZE`)

Calculate non-overlapping GC-content windows or identify candidate CpG islands
using 200-bp seed windows with GC and observed/expected CpG thresholds:

```sql
ANALYZE GC_CONTENT WINDOW 1 KB AS gc_profile;
ANALYZE CPG_ISLANDS AS cpg_islands;
```

### 5. Interval Operations (`INTERSECT`, `UNION`, `EXCEPT`, `OVERLAPS`)

Combine or filter interval sets using high-speed interval algebra:

```sql
INTERSECT sp1_sites AND cpg_islands AS supported_sites;
UNION minus35_box AND minus10_box AS promoter_boxes;
EXCEPT ctcf_sites FROM CDS AS noncoding_ctcf_sites;
OVERLAPS myb_sites WITH candidate_promoters AS promoter_myb_sites;
```

`INTERSECT` emits the clipped overlap geometry. `OVERLAPS query WITH reference`
instead performs a directional semi-join: each query interval is retained once
if any reference interval overlaps it. This preserves the query coordinates,
sequence, and motif evidence.

### 6. Feature Extraction & Filtering (`EXTRACT`, `WHERE`)

Filter genomic entities by physical length or alignment similarity:

```sql
EXTRACT GENE AS reference_gene WHERE ID = "geneA";
EXTRACT GENE AS homologous_genes
    WHERE LENGTH > 1 KB
      AND SIMILARITY TO reference_gene > 70 %;
```

Condition properties are result-specific: region sets support `LENGTH`,
`SIMILARITY`, `GC_CONTENT`, and `ID`; motif results support `LENGTH` and
`GC_CONTENT`; GC profiles support `GC_CONTENT`. `IF` currently evaluates the
GC content of the active sequence dataset.

An explicit similarity reference must be a named result set containing exactly
one region. Similarity percentages are normalized local-alignment scores, not
percent identity.

### 7. Result Export (`EXPORT`)

```sql
EXPORT homologous_genes TO "homologous_genes.bed" FORMAT BED;
EXPORT homologous_genes TO "homologous_genes.gff3" FORMAT GFF3;
EXPORT homologous_genes TO "homologous_genes.tsv" FORMAT TSV;
EXPORT gc_profile TO "gc_profile.tsv" FORMAT TSV;
```

Internal and BED coordinates are zero-based and half-open. GFF3 exports convert
the start coordinate to the one-based inclusive convention. Export paths are
relative to the query workspace; absolute paths and parent-directory traversal
are rejected.

### 8. Control Flow (`IF`, `FOREACH`)

Control query execution paths and iterate over collections of loaded matrices:

```sql
// Conditional execution based on sequence properties
IF GC_CONTENT > 50 % THEN
    SCAN sp1 THRESHOLD 80 % AS gc_sites;
ELSE
    SCAN tbp THRESHOLD 80 % AS at_sites;
ENDIF;

// Batch processing over matrix lists
FOREACH m IN [tbp, sp1, ctcf] DO
    SCAN m THRESHOLD 80 % AS tf_sites;
ENDFOR;
```

---

## Formal Syntax & Grammar (CFG)

Cis-QL is formally specified by an LL(1) Context-Free Grammar. Below is the complete EBNF specification matching `GRAMMAR.TXT`:

```ebnf
Program            ::= StatementList
StatementList      ::= Statement StatementList | λ

Statement          ::= LoadStmt | UseStmt | ExportStmt | FindStmt | ExtractStmt
                     | SetOperationStmt | ScanStmt | AnalyzeStmt | IfStmt
                     | ForeachStmt

LoadStmt           ::= LOAD (SEQUENCE | ANNOTATION | MATRIX) STRING AS ID SEMICOLON
UseStmt            ::= USE (SEQUENCE | ANNOTATION) ID SEMICOLON
ExportStmt         ::= EXPORT ID TO STRING FORMAT (BED | GFF3 | TSV) SEMICOLON

AnalyzeStmt        ::= ANALYZE (GC_CONTENT | CPG_ISLANDS) (WINDOW (NUM | FLOAT) Unit)? AliasOpt WhereClause SEMICOLON

FindStmt           ::= FIND MOTIF STRING FindOpts AliasOpt WhereClause SEMICOLON
FindOpts           ::= FindOpt FindOpts | λ
FindOpt            ::= WITHIN (NUM | FLOAT) Unit Direction FROM EntityRef EntityName
                     | STRAND StrandType
                     | CHR STRING

ScanStmt           ::= SCAN ID ScanOpts AliasOpt WhereClause SEMICOLON
ScanOpts           ::= ScanOpt ScanOpts | λ
ScanOpt            ::= STRAND StrandType
                     | THRESHOLD (NUM | FLOAT) PERCENT

SetOperationStmt   ::= (INTERSECT | UNION) EntityRef AND EntityRef AliasOpt WhereClause SEMICOLON
                     | EXCEPT EntityRef FROM EntityRef AliasOpt WhereClause SEMICOLON

ExtractStmt        ::= EXTRACT EntityRef AliasOpt WhereClause SEMICOLON

IfStmt             ::= IF Condition THEN StatementList (ELSE StatementList)? ENDIF (SEMICOLON)?

ForeachStmt        ::= FOREACH ID IN "[" CollectionList "]" DO StatementList ENDFOR (SEMICOLON)?
CollectionList     ::= CollectionItem ("," CollectionItem)* | λ
CollectionItem     ::= ID

WhereClause        ::= WHERE Condition | λ

Condition          ::= Term ConditionPrime
ConditionPrime     ::= OR Term ConditionPrime | λ
Term               ::= Factor TermPrime
TermPrime          ::= AND Factor TermPrime | λ
Factor             ::= NOT Factor | SimpleCondition | "(" Condition ")"

SimpleCondition    ::= Property RelOp Value
                     | SIMILARITY SimilarityRefOpt RelOp Value
Property           ::= LENGTH | GC_CONTENT | ID
SimilarityRefOpt   ::= TO ID | λ
RelOp              ::= ">" | "<" | ">=" | "<=" | "="
Value              ::= (NUM | FLOAT) Unit | (NUM | FLOAT) PERCENT | NUM | FLOAT | STRING

Unit               ::= BP | KB | MB | λ
Direction          ::= UPSTREAM | DOWNSTREAM
Entity             ::= GENE | PROMOTER | ENHANCER | EXON | INTRON | UTR | TSS | CDS | REGION
EntityRef          ::= Entity | ID
EntityName         ::= STRING | λ
AliasOpt           ::= AS ID | λ
StrandType         ::= POSITIVE | NEGATIVE
```

## Correctness Tests and Engineering Benchmarks

Run the automated language and algorithm tests with:

```bash
make test
make validate
```

`make validate` compares motif coordinates, PSSM thresholds, interval
subtraction, and normalized local-alignment filtering with independent
reference implementations. If BEDTools or Biopython are installed, compatible
external checks run as additional optional comparisons.

Deterministic scaling measurements and complete-query timings are available
with:

```bash
make benchmark-core
make benchmark
```

These are engineering measurements of the current implementation, not evidence
of superiority over other tools. See `benchmarks/README.md` for the comparison
requirements needed before reporting external benchmark results.

---

## Curated Examples Suite (`cql_examples/`)

The repository includes 19 structured `.cql` scripts demonstrating specific language capabilities:

The regulatory progression and its expected outputs are described in
[`cql_examples/README.md`](cql_examples/README.md).

| Script | Description | Primary Features |
| :--- | :--- | :--- |
| `01_extract_genes.cql` | Structural gene extraction | `EXTRACT`, `WHERE LENGTH` |
| `02_length_filter.cql` | Multi-condition range filtering | Logical `AND`, bounded intervals |
| `03_similarity_align.cql` | Paralog discovery via local alignment | Smith-Waterman `SIMILARITY > 70 %` |
| `04_iupac_motifs.cql` | Degenerate promoter motif search | IUPAC translation engine (`TATAWAW`) |
| `05_spatial_promoters.cql` | Upstream regulatory element search | `WITHIN 200 BP UPSTREAM FROM` |
| `06_denovo_orfs.cql` | Unannotated ORF discovery | Regex matching, virtual annotations |
| `07_pwm_scanning.cql` | Plant MYB matrix scanning | JASPAR PFM, relative log-odds threshold |
| `08_ctcf_insulators.cql` | Chromatin insulator mapping | `EXCEPT` set subtraction |
| `09_cpg_islands.cql` | Epigenetic CpG island profiling | `ANALYZE CPG_ISLANDS`, `INTERSECT` |
| `10_promoter_union.cql` | Bipartite promoter element merger | Multi-track consolidation via `UNION` |
| `11_strand_search.cql` | Sense vs. antisense motif profiling | `STRAND POSITIVE / NEGATIVE` |
| `12_gc_content.cql` | Sliding-window GC landscape | `ANALYZE GC_CONTENT WINDOW` |
| `13_complex_where.cql` | Multi-property conditional queries | Combined `LENGTH` & `SIMILARITY` |
| `14_integrated_query.cql` | Anthocyanin regulatory evidence map | Promoters, calibrated MYB sites, regulatory overlap |
| `15_if_else_branching.cql` | Conditional flow execution | `IF / ELSE / ENDIF` |
| `16_foreach_batch_scan.cql` | Batch processing over matrix lists | `FOREACH / DO / ENDFOR` |
| `17_explicit_promoters.cql` | TSS-oriented candidate promoters | `DEFINE PROMOTERS`, explicit bounds |
| `18_statistical_pwm_scan.cql` | Calibrated plant MYB-site selection | `BACKGROUND FROM`, `QVALUE` |
| `19_regulatory_overlap.cql` | Promoter-supported MYB evidence | Directional `OVERLAPS`, evidence preservation |

---

## Cis-QL Studio Interface

Cis-QL Studio provides a desktop graphical environment for query development, execution, and visual exploration:

- **Integrated Code Editor & Console:** Write, load, and execute `.cql` queries with real-time terminal output.
- **Multi-Track Visualizer Canvas:** Displays ruler coordinates, GC content profiles, feature annotation tracks, and sequence details.
- **Data Synchronization:** Automatically synchronizes execution results via `.cisql_results.json` for live inspection.

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

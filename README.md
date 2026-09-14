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
- **Regulatory Interval Queries:** Supports geometric interval algebra, directional `OVERLAPS`, evidence-preserving nearest-feature selection, and per-region overlap counts with explicit provenance.
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
./cisql cql_examples/01_define_promoters.cql
```

Use the `--debug` flag to inspect compilation phases, including token stream, Abstract Syntax Tree (AST), Symbol Table, Intermediate Representation (IR), and execution steps:

```bash
./cisql cql_examples/08_integrated_anthocyanin_query.cql --debug
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
LOAD SEQUENCE "data_examples/anthocyanin_regulatory_demo.fasta" AS genome;
LOAD ANNOTATION "data_examples/anthocyanin_regulatory_demo.gff3" AS annotation;
LOAD MATRIX "matrices/MA0054.1_myb.Ph3.pwm" AS myb_matrix;
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
FIND MOTIF "TAACCGTTT" STRAND POSITIVE AS exact_myb_like_sites;

FIND MOTIF "WAACCGTTT"
    WITHIN 100 BP UPSTREAM FROM GENE
    AS degenerate_upstream_sites;
```

`FIND MOTIF` is appropriate for exact strings, regular expressions, or
documented IUPAC patterns. A short consensus is not equivalent to a TF-binding
model; use `SCAN` with a sourced PWM when TF specificity and calibrated scores
matter.

### 3. Position Weight Matrix Scanning (`SCAN`)

Scan loaded sequences using Position Weight Matrices with log-odds scoring:

```sql
SCAN myb_matrix THRESHOLD 90 % AS high_scoring_myb_sites;

SCAN myb_matrix BACKGROUND FROM genome
    QVALUE <= 0.01 AS supported_myb_sites;
```

### 4. Biological & Structural Analysis (`ANALYZE`)

Calculate non-overlapping GC-content windows or identify candidate CpG islands
using 200-bp seed windows with GC and observed/expected CpG thresholds:

```sql
ANALYZE GC_CONTENT WINDOW 1 KB AS gc_profile;
ANALYZE CPG_ISLANDS AS cpg_islands;
```

### 5. Interval Operations (`INTERSECT`, `UNION`, `EXCEPT`, `OVERLAPS`, `NEAR`)

Combine or filter interval sets using high-speed interval algebra:

```sql
OVERLAPS supported_myb_sites WITH candidate_promoters AS promoter_myb_sites;
NEAR supported_myb_sites TO GENE WITHIN 2 KB AS proximal_gene_candidates;
COUNT supported_myb_sites IN candidate_promoters AS promoter_site_counts;
EXTRACT promoter_site_counts AS supported_promoters WHERE COUNT >= 1;
```

`INTERSECT` emits the clipped overlap geometry. `OVERLAPS query WITH reference`
instead performs a directional semi-join: each query interval is retained once
if any reference interval overlaps it. This preserves the query coordinates,
sequence, and motif evidence.

`NEAR query TO reference WITHIN distance` also preserves each complete query
record, but retains it only when its nearest same-chromosome reference is no
farther than the inclusive limit. It records that reference's identity and
coordinates, the interval gap, the requested limit, and whether the intervals
actually overlap. Overlap and directly touching boundaries both have gap 0,
but the `overlaps` field distinguishes them. Equal-distance ties are resolved
deterministically by reference coordinate and metadata. Proximity is a
candidate association—not evidence of regulation by itself.

`COUNT query IN containers` emits every container, including zero-count
regions, and records the number of query intervals that overlap it. The count,
counted set, container set, and overlap relation remain available in JSON,
GFF3, TSV, and Studio. `WHERE COUNT` can then select regions by support. These
are descriptive counts, not motif enrichment: region length, composition,
accessibility, and the statistical background still need explicit treatment.
Counting is strand-agnostic and operates on input records; duplicate records
are counted separately. Restrict or normalize the counted set first when that
is not the intended universe.

### 6. Cis-regulatory modules (`DEFINE MODULE`)

Pair motif or region records under explicit spacing, order, and orientation
constraints:

```sql
DEFINE MODULE
    FROM enhancer_myb_sites WITH enhancer_myb_sites
    SPACING 5 BP TO 30 BP
    ORDER ANY
    ORIENTATION ANY
    AS enhancer_myb_modules;
```

Spacing is the gap between half-open intervals and both bounds are inclusive;
overlapping and adjacent intervals have gap zero. `ORDER AS_WRITTEN` requires a
member from the first set to have a lower reference start than one from the
second set; tied starts do not satisfy it. `ORDER ANY` accepts either order and
records what was observed. `SAME` and `OPPOSITE` require known `+`/`-` strands.
When one alias is used twice, Cis-QL excludes self-pairs and mirror duplicates.
The output span
retains both member identities and PWM evidence in JSON, GFF3, TSV, and Studio.

This constructs candidates satisfying a declared grammar; it does not prove
cooperative binding. Appropriate spacing bounds should come from a stated
hypothesis, reference set, or sensitivity analysis.

### 7. Feature Extraction & Filtering (`EXTRACT`, `WHERE`)

Filter genomic entities by physical length or alignment similarity:

```sql
EXTRACT GENE AS reference_gene WHERE ID = "geneA";
EXTRACT GENE AS homologous_genes
    WHERE LENGTH > 1 KB
      AND SIMILARITY TO reference_gene > 70 %;
```

Condition properties are result-specific: region sets support `LENGTH`,
`SIMILARITY`, `GC_CONTENT`, and `ID`, plus `COUNT` when count evidence is
attached; motif results support `LENGTH` and `GC_CONTENT`; GC profiles support
`GC_CONTENT`. `IF` currently evaluates the GC content of the active sequence
dataset.

An explicit similarity reference must be a named result set containing exactly
one region. Similarity percentages are normalized local-alignment scores, not
percent identity.

### 8. Result Export (`EXPORT`)

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

### 9. Control Flow (`IF`, `FOREACH`)

Control query execution paths and iterate over collections of loaded matrices:

```sql
FOREACH m IN [myb_a, myb_b, myb_c] DO
    SCAN m BACKGROUND FROM genome QVALUE <= 0.01 AS tf_sites;
ENDFOR;
```

`IF` is also supported for execution control. The curated examples avoid using
whole-genome composition to choose a TF model or threshold because that would
hide a scientific decision inside a convenient branch.

---

## Formal Syntax & Grammar (CFG)

Cis-QL is specified by the development CFG in [`grammar.txt`](grammar.txt).
The following compact EBNF lists the main statement forms; the linked file is
the authoritative grammar:

```ebnf
Program            ::= StatementList
StatementList      ::= Statement StatementList | λ

Statement          ::= LoadStmt | UseStmt | ExportStmt | FindStmt | ExtractStmt
                     | DefinePromotersStmt | DefineModuleStmt
                     | SetOperationStmt | CountStmt | ScanStmt | AnalyzeStmt
                     | IfStmt | ForeachStmt

LoadStmt           ::= LOAD (SEQUENCE | ANNOTATION | MATRIX) STRING AS ID SEMICOLON
UseStmt            ::= USE (SEQUENCE | ANNOTATION) ID SEMICOLON
ExportStmt         ::= EXPORT ID TO STRING FORMAT (BED | GFF3 | TSV) SEMICOLON
DefinePromotersStmt ::= DEFINE PROMOTERS OF (GENE | TSS | ID) FROM TSS
                        UPSTREAM (NUM | FLOAT) RequiredUnit
                        DOWNSTREAM (NUM | FLOAT) RequiredUnit AS ID SEMICOLON
DefineModuleStmt   ::= DEFINE MODULE FROM ID WITH ID
                       SPACING (NUM | FLOAT) RequiredUnit TO
                               (NUM | FLOAT) RequiredUnit
                       ORDER (ANY | AS_WRITTEN)
                       ORIENTATION (ANY | SAME | OPPOSITE) AS ID SEMICOLON

AnalyzeStmt        ::= ANALYZE (GC_CONTENT | CPG_ISLANDS) (WINDOW (NUM | FLOAT) Unit)? AliasOpt WhereClause SEMICOLON

FindStmt           ::= FIND MOTIF STRING FindOpts AliasOpt WhereClause SEMICOLON
FindOpts           ::= FindOpt FindOpts | λ
FindOpt            ::= WITHIN (NUM | FLOAT) Unit Direction FROM EntityRef EntityName
                     | STRAND StrandType
                     | CHR STRING

ScanStmt           ::= SCAN ID ScanOpts AliasOpt WhereClause SEMICOLON
ScanOpts           ::= ScanOpt ScanOpts | λ
ScanOpt            ::= IN EntityRef
                     | STRAND StrandType
                     | THRESHOLD (NUM | FLOAT) PERCENT
                     | (PVALUE | QVALUE) (LESS | LESS_EQ) Probability
                     | BACKGROUND (UNIFORM | FROM EntityRef)

SetOperationStmt   ::= (INTERSECT | UNION) EntityRef AND EntityRef AliasOpt WhereClause SEMICOLON
                     | EXCEPT EntityRef FROM EntityRef AliasOpt WhereClause SEMICOLON
                     | OVERLAPS EntityRef WITH EntityRef AliasOpt WhereClause SEMICOLON
                     | NEAR EntityRef TO EntityRef WITHIN
                         (NUM | FLOAT) RequiredUnit AliasOpt WhereClause SEMICOLON

CountStmt          ::= COUNT EntityRef IN EntityRef AS ID WhereClause SEMICOLON

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
Property           ::= LENGTH | GC_CONTENT | COUNT | ID | NAME
SimilarityRefOpt   ::= TO ID | λ
RelOp              ::= ">" | "<" | ">=" | "<=" | "="
Value              ::= (NUM | FLOAT) Unit | (NUM | FLOAT) PERCENT | NUM | FLOAT | STRING
Probability        ::= NUM | FLOAT

Unit               ::= BP | KB | MB | λ
RequiredUnit       ::= BP | KB | MB
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

`make validate` compares motif coordinates, PSSM thresholds and tail
probabilities, interval subtraction, overlap/nearest selection, overlap counts,
constrained homotypic modules, and normalized local-alignment filtering with independent reference
implementations. If BEDTools, Biopython, or FIMO are installed, compatible
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

The repository includes eight `.cql` analyses forming one coherent,
synthetic anthocyanin-regulation path. Compiler feature coverage belongs in
the automated tests; these programs are examples of scientific questions.

The regulatory progression and its expected outputs are described in
[`cql_examples/README.md`](cql_examples/README.md).

| Script | Description | Primary Features |
| :--- | :--- | :--- |
| `01_define_promoters.cql` | Build explicit candidate promoter windows | `DEFINE PROMOTERS`, TSS-relative bounds |
| `02_score_myb_sites.cql` | Inspect high-scoring plant MYB matches | sourced PWM, relative score threshold |
| `03_calibrated_myb_sites.cql` | Select statistically supported MYB sites | genomic background, `QVALUE` |
| `04_promoter_supported_sites.cql` | Retain promoter-overlapping MYB evidence | directional `OVERLAPS` |
| `05_count_promoter_support.cql` | Summarize site support per promoter | zero-preserving `COUNT`, `WHERE COUNT` |
| `06_nearest_gene_candidates.cql` | Form proximity-based gene hypotheses | auditable `NEAR ... WITHIN` |
| `07_enhancer_myb_module.cql` | Detect a constrained homotypic MYB module | spacing, order, orientation, two-member evidence |
| `08_integrated_anthocyanin_query.cql` | Connect the complete evidence path | promoters, calibrated sites, modules, counts, candidate links |

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
- Evidence-preserving regulatory selection, counting, and motif modules.
- Control flow execution engine (`IF/ELSE` and `FOREACH`).
- Desktop GUI workspace (`Cis-QL Studio`).

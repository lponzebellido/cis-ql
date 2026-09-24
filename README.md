# Cis-QL: a query language for cis-regulatory analysis

Cis-QL is an experimental domain-specific language for describing
cis-regulatory analyses as readable, executable queries. A `.cql` program
names its sequence, annotation, motif models, and experimental tracks, then
states how their evidence should be related: inside a promoter, overlapping an
accessible region, supported across replicates, arranged as a motif module, or
near a possible target gene.

The project is motivated by a practical problem. These questions are usually
spread across interval commands, motif-scanning tools, short scripts, and
manually interpreted tables. Cis-QL makes the regulatory relationships part
of the program and carries their provenance into the result. The goal is not
to replace established upstream tools, but to provide a compact language for
the evidence-integration step between their outputs and a testable regulatory
hypothesis.

The command-line interpreter is written in C++11. Cis-QL Studio adds a desktop
editor and result viewer, but `.cql` files do not depend on the graphical
interface.

## What the language does today

The current implementation reads FASTA sequences, hierarchy-preserving GFF3
annotations, JASPAR frequency matrices, and BED or narrowPeak tracks. It can
derive strand-aware promoters, find literal or IUPAC sequence patterns, scan
PWMs with explicit background and multiple-testing metadata, combine genomic
intervals, describe two-site motif modules, count support, select nearby
features, and build anchor-preserving consensus sets across experimental
tracks. Results can be
inspected as structured JSON and exported to BED, GFF3, or TSV.

A typical query has this shape:

```sql
LOAD SEQUENCE "genome.fasta" AS genome;
LOAD ANNOTATION "genes.gff3" AS annotation;
LOAD MATRIX "tf_model.pwm" AS tf_model;
LOAD TRACK "binding_rep1.narrowPeak" FORMAT NARROWPEAK
    EVIDENCE BINDING REPLICATE "R1" AS binding_rep1;
LOAD TRACK "binding_rep2.narrowPeak" FORMAT NARROWPEAK
    EVIDENCE BINDING REPLICATE "R2" AS binding_rep2;

EXTRACT GENE AS candidate_genes WHERE ID = "candidate_1";
DEFINE PROMOTERS OF candidate_genes FROM TSS
    UPSTREAM 1000 BP DOWNSTREAM 100 BP AS candidate_promoters;
SCAN tf_model IN candidate_promoters BACKGROUND FROM genome
    QVALUE <= 0.01 AS promoter_sites;
CONSENSUS FROM [binding_rep1, binding_rep2]
    ANCHOR binding_rep1 MIN_SUPPORT 2
    MIN_RECIPROCAL_OVERLAP 50 %
    MAX_SUMMIT_DISTANCE 20 BP AS reproducible_binding;
OVERLAPS reproducible_binding WITH promoter_sites AS supported_binding;
```

The filenames and biological choices in that fragment are intentionally
generic. Cis-QL does not encode a particular pathway, transcription-factor
family, organism, or assay. The repository contains both a compact eukaryotic
evidence-integration fixture and a bacterial sequence-pattern example. They
exercise different parts of the language and are not its biological boundary.

## Project status and boundaries

Cis-QL is a research prototype, not a complete genomics platform. It does not
align reads, call peaks, infer enhancers, perform IDR, prove a target-gene
relationship, or decide whether a PWM and threshold are biologically suitable.
Those choices remain explicit inputs to the analysis. The value of the
language is that they can be stated, reviewed, rerun, and exported together
instead of disappearing inside an ad hoc pipeline.

The prioritized implementation gaps and their acceptance targets are tracked
in [`docs/RFC_CISQL_V3.md`](docs/RFC_CISQL_V3.md#work-remaining-after-the-regulatory-foundation).

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
./cisql cql_examples/12_regex_denovo.cql --debug
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
LOAD SEQUENCE "genome.fasta" AS genome;
LOAD ANNOTATION "genes.gff3" AS annotation;
LOAD MATRIX "tf_model.pwm" AS tf_model;
LOAD TRACK "accessibility.narrowPeak"
    FORMAT NARROWPEAK
    EVIDENCE ACCESSIBILITY
    ASSAY "ATAC-seq"
    SAMPLE "sample_1"
    CONDITION "treated"
    REPLICATE "R1"
    CONTROL "input"
    AS accessibility_peaks;
```

The example suite also includes the plant MYB profile `MA0054.1` from
[JASPAR CORE](https://jaspar.elixir.no/matrix/MA0054.1/). The accompanying
`anthocyanin_regulatory_demo` FASTA/GFF3 pair is synthetic and exists only to
make expected regulatory-query results small, deterministic, and inspectable.

GFF3 import preserves source, score, phase, `ID`, `Name`, every `Parent`, and
the complete attribute map. Percent-encoded attribute values are decoded for
queries, while GFF3 export retains the original attribute field. Cis-QL does
not guess which Sequence Ontology types constitute a transcript; use
`EXTRACT FEATURE` with explicit type and attribute filters instead.

`LOAD TRACK` currently accepts BED and narrowPeak. Track coordinates are
already zero-based and half-open. When a sequence dataset is active, Cis-QL
validates chromosome names and interval bounds, attaches interval sequence,
and rejects assembly-incompatible records. narrowPeak's `pValue` and `qValue`
columns are preserved according to that format as `-log10(p)` and `-log10(q)`;
they are not confused with the calibrated probabilities produced by `SCAN`.
Every track must declare `EVIDENCE ACCESSIBILITY`, `BINDING`, or `OTHER`.
Optional `ASSAY`, `SAMPLE`, `CONDITION`, `REPLICATE`, and `CONTROL` strings
may appear in any order and travel with derived results and exports. Duplicate
metadata clauses are rejected. These declarations preserve provenance but do
not validate experimental quality or biological interpretation.

When multiple sequence or annotation datasets are loaded, select the active
context explicitly:

```sql
USE SEQUENCE genome;
USE ANNOTATION annot;
```

### 2. Motif Searching & Spatial Conditions (`FIND MOTIF`)

Locate exact motifs, regular expressions, or IUPAC degenerate strings, with optional spatial constraints relative to other features:

```sql
FIND MOTIF "TATA[AT]A[AT]" STRAND POSITIVE AS promoter_like_patterns;

FIND MOTIF "CANNTG"
    WITHIN 100 BP UPSTREAM FROM GENE
    AS degenerate_upstream_patterns;

FIND MOTIF "ATG(?:(?!TAA|TAG|TGA)[ACGT]{3})*(?:TAA|TAG|TGA)"
    STRAND POSITIVE AS start_stop_candidates
    WHERE LENGTH >= 15 BP AND LENGTH MOD 3 = 0;
```

`FIND MOTIF` is appropriate for exact strings, regular expressions, or
documented IUPAC patterns. A short consensus is not equivalent to a TF-binding
model; use `SCAN` with a sourced PWM when TF specificity and calibrated scores
matter.

The start-to-stop expression consumes complete codons and excludes an internal
in-frame stop. The explicit `LENGTH MOD 3 = 0` condition makes the frame
invariant visible to the reader instead of leaving it implicit in the regex.
It reports sequence candidates, not predicted genes; genetic code, alternative
starts, minimum coding length, annotation evidence, and the biological question
remain separate policies.

### 3. Position Weight Matrix Scanning (`SCAN`)

Scan loaded sequences using Position Weight Matrices with log-odds scoring:

```sql
SCAN tf_model THRESHOLD 90 % AS high_scoring_tf_sites;

SCAN tf_model BACKGROUND FROM genome
    QVALUE <= 0.01 AS significant_tf_sites;
```

### 4. Biological & Structural Analysis (`ANALYZE`)

Calculate non-overlapping GC-content windows or identify candidate CpG islands
using 200-bp seed windows with GC and observed/expected CpG thresholds:

```sql
ANALYZE GC_CONTENT WINDOW 1 KB AS gc_profile;
ANALYZE CPG_ISLANDS AS cpg_islands;
```

### 5. Interval Operations (`INTERSECT`, `UNION`, `EXCEPT`, `OVERLAPS`, `NEAR`, `CONSENSUS`)

Combine or filter interval sets using high-speed interval algebra:

```sql
OVERLAPS significant_tf_sites WITH candidate_promoters AS promoter_tf_sites;
NEAR significant_tf_sites TO GENE WITHIN 2 KB AS proximal_gene_candidates;
COUNT significant_tf_sites IN candidate_promoters AS promoter_site_counts;
EXTRACT promoter_site_counts AS supported_promoters WHERE COUNT >= 1;

CONSENSUS FROM [binding_rep1, binding_rep2]
    ANCHOR binding_rep1 MIN_SUPPORT 2
    MIN_RECIPROCAL_OVERLAP 50 %
    MAX_SUMMIT_DISTANCE 20 BP AS reproducible_binding;
```

`INTERSECT` emits the clipped overlap geometry. `OVERLAPS query WITH reference`
instead performs a directional semi-join: each query interval is retained once
if any reference interval overlaps it. This preserves the query coordinates,
sequence, and motif evidence. Every matching reference is recorded in the
result's `overlapEvidence`; when a reference came from `LOAD TRACK`, its
evidence class, experimental metadata, score, signal, and summit remain
attached. If that reference was itself supported by an earlier `OVERLAPS`,
its prior evidence is retained recursively as `supportingEvidence`; this
preserves the provenance path without claiming that a transitive relationship
is a direct overlap.

`CONSENSUS FROM [sets] ANCHOR set MIN_SUPPORT n` retains complete intervals
from the named anchor when they overlap records from enough distinct input
sets. The anchor itself contributes one unit of support; each other set
contributes at most one, regardless of how many of its records overlap. The
optional `MIN_RECIPROCAL_OVERLAP p %` clause requires each accepted match to
cover at least `p` percent of both the anchor interval and the supporting
interval, preventing a marginal one-base overlap from counting as replicate
support. `MAX_SUMMIT_DISTANCE d` additionally requires both matching records
to provide narrowPeak summits no farther apart than `d`; records without a
summit cannot satisfy that criterion. The two optional clauses can be combined
and written in either order. The result records the criteria, requested and
observed support, all input aliases, every matching observation, and the
anchor geometry.
`SUPPORT_COUNT` makes the observed number queryable in `WHERE`. This remains
coordinate-level concordance, not IDR or a statistical reproducibility test.

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
    FROM tf_a_sites WITH tf_b_sites
    SPACING 5 BP TO 30 BP
    ORDER ANY
    ORIENTATION ANY
    AS candidate_regulatory_modules;
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

Filter genomic entities and motif hits by coordinates, orientation, physical
length, alignment similarity, or attached evidence:

```sql
EXTRACT GENE AS reference_gene WHERE ID = "geneA";
EXTRACT GENE AS homologous_genes
    WHERE LENGTH > 1 KB
      AND SIMILARITY TO reference_gene > 70 %;

EXTRACT FEATURE AS coding_features
    WHERE TYPE = "CDS"
      AND PARENT = "transcript_1"
      AND PHASE = "0"
      AND ATTRIBUTE "protein_id" = "protein_1";

EXTRACT accessibility_peaks AS strong_accessibility
    WHERE TRACK_SCORE >= 600
      AND SIGNAL_VALUE >= 10
      AND EVIDENCE_CLASS = "ACCESSIBILITY"
      AND CONDITION = "treated"
      AND REPLICATE = "R1";

FIND MOTIF "ATG" AS phase_zero_starts
    WHERE START MOD 3 = 0
      AND END MOD 3 = 0
      AND STRAND = "+";
```

Condition properties are result-specific. Region sets support `LENGTH`,
`START`, `END`, `STRAND`, `SIMILARITY`, `GC_CONTENT`, and `ID`, plus `COUNT`
when count evidence is attached. Annotation-backed regions also support
`NAME`, `TYPE`, `PARENT`, `SOURCE`, `PHASE`, and arbitrary
`ATTRIBUTE "key" = "value"` filters. `PARENT` matches any member of a
multi-parent GFF3 record. Track-backed regions additionally support
`TRACK_SCORE`,
`SIGNAL_VALUE`, `MINUS_LOG10_PVALUE`, `MINUS_LOG10_QVALUE`,
`EVIDENCE_CLASS`, `ASSAY`, `SAMPLE`, `CONDITION`, `REPLICATE`, and `CONTROL`.
Missing optional narrowPeak values do not satisfy a numeric condition.
Metadata uses exact string equality. These properties filter the primary
track; filter a reference track before combining it with `OVERLAPS`. Motif
results support `LENGTH`, `START`, `END`, `STRAND`, and `GC_CONTENT`; GC
profiles support `GC_CONTENT`.
`IF` currently evaluates the GC content of the active sequence dataset.
Regions emitted by `CONSENSUS` additionally support the non-negative integer
property `SUPPORT_COUNT`.

`MOD` can be applied to numeric properties before comparison, for example
`LENGTH MOD 3 = 0` or `START MOD 3 = 0`. Coordinates are zero-based and
half-open. Modular filtering is a general arithmetic constraint; its biological
meaning comes from the surrounding program.

These filters do not calibrate experimental evidence. BED/narrowPeak scores
and `signalValue` remain upstream-tool-specific, so thresholds require an
assay-appropriate justification and should not be compared across experiments
as if they shared a universal scale.

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
FOREACH m IN [tf_model_a, tf_model_b, tf_model_c] DO
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
                     | SetOperationStmt | ConsensusStmt | CountStmt
                     | ScanStmt | AnalyzeStmt
                     | IfStmt | ForeachStmt

LoadStmt           ::= LOAD (SEQUENCE | ANNOTATION | MATRIX) STRING AS ID SEMICOLON
                     | LOAD TRACK STRING FORMAT (BED | NARROWPEAK)
                       EVIDENCE (ACCESSIBILITY | BINDING | OTHER)
                       TrackMetadata* AS ID SEMICOLON
TrackMetadata      ::= ASSAY STRING | SAMPLE STRING | CONDITION STRING
                     | REPLICATE STRING | CONTROL STRING
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

ConsensusStmt      ::= CONSENSUS FROM "[" ID ("," ID)+ "]"
                       ANCHOR ID MIN_SUPPORT NUM ConsensusCriterion*
                       AS ID WhereClause SEMICOLON
ConsensusCriterion ::= MIN_RECIPROCAL_OVERLAP (NUM | FLOAT) PERCENT
                     | MAX_SUMMIT_DISTANCE (NUM | FLOAT) RequiredUnit

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

SimpleCondition    ::= NumericProperty NumericModifierOpt RelOp Value
                     | SIMILARITY SimilarityRefOpt NumericModifierOpt RelOp Value
                     | StringProperty RelOp STRING
                     | ATTRIBUTE STRING RelOp STRING
NumericModifierOpt ::= MOD (NUM | FLOAT) | λ
NumericProperty    ::= LENGTH | START | END | GC_CONTENT | COUNT
                     | SUPPORT_COUNT | TRACK_SCORE | SIGNAL_VALUE
                     | MINUS_LOG10_PVALUE | MINUS_LOG10_QVALUE
StringProperty     ::= ID | NAME | STRAND | TYPE | PARENT | SOURCE | PHASE
                     | EVIDENCE_CLASS | ASSAY | SAMPLE | CONDITION
                     | REPLICATE | CONTROL
SimilarityRefOpt   ::= TO ID | λ
RelOp              ::= ">" | "<" | ">=" | "<=" | "="
Value              ::= (NUM | FLOAT) Unit | (NUM | FLOAT) PERCENT | NUM | FLOAT | STRING
Probability        ::= NUM | FLOAT

Unit               ::= BP | KB | MB | λ
RequiredUnit       ::= BP | KB | MB
Direction          ::= UPSTREAM | DOWNSTREAM
Entity             ::= GENE | PROMOTER | ENHANCER | EXON | INTRON | UTR | TSS | CDS | REGION | FEATURE
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

The repository includes runnable `.cql` programs grouped by the capability
they demonstrate. Scripts 01-11 form one compact, synthetic eukaryotic
evidence-integration workflow. Scripts 12-13 use an *E. coli* reference to show
sequence-pattern constraints and hierarchy-preserving GFF3 queries. Neither
fixture defines the intended organism, pathway, or regulatory model of Cis-QL.

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
| `09_accessible_myb_evidence.cql` | Combine imported accessibility peaks with motif support | narrowPeak provenance, `COUNT`, `NEAR` |
| `10_accessible_bound_myb_candidates.cql` | Combine filtered accessibility, binding, motif, and proximity | track-evidence `WHERE`, multi-track `overlapEvidence`, `COUNT`, `NEAR` |
| `11_replicate_supported_candidates.cql` | Require substantial coordinate and summit support from two binding replicates | `CONSENSUS`, reciprocal overlap, summit distance, structured experimental provenance |
| `12_regex_denovo.cql` | Find in-frame start-to-stop sequence candidates near TATA-like anchors | regex search, strand, `LENGTH MOD 3` |
| `13_gff3_hierarchy.cql` | Select a real gene and its CDS through preserved annotation identity | `FEATURE`, `ID`, `PARENT`, `PHASE`, `ATTRIBUTE` |

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

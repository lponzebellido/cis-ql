# Cis-QL

Cis-QL is a query language for sequence, annotation, motif, and genomic
interval analysis. A `.cql` file loads FASTA, GFF3, JASPAR matrix, BED, or
narrowPeak data and applies named operations to it. The C++11 interpreter runs
from the command line; Cis-QL Studio provides an optional editor and result
viewer.

## What it does

| Area | Operations |
| :--- | :--- |
| Sequence | literal, IUPAC, and regular-expression search; PWM scanning; GC analysis |
| Annotation | GFF3 validation, filtering, hierarchy traversal, repeated-ID grouping, promoter derivation |
| Intervals | set operations, overlap, counting, distance, two-pattern modules, replicate consensus |
| Results | structured JSON and BED, GFF3, or TSV export |

A query can combine these operations without hiding their parameters:

```sql
LOAD SEQUENCE "data_examples/ecoli2.fna" AS genome;
LOAD ANNOTATION "data_examples/genomic.gff" AS annotation;

EXTRACT GENE AS thra WHERE ID = "gene-b0002";
DEFINE PROMOTERS OF thra FROM TSS
    UPSTREAM 60 BP DOWNSTREAM 20 BP AS promoters;
FIND MOTIF "TATAAT" AS sites;
OVERLAPS sites WITH promoters AS promoter_sites;
EXPORT promoter_sites TO "promoter_sites.tsv" FORMAT TSV;
```

Cis-QL is a research prototype. It does not align reads, call peaks, perform
IDR, infer enhancers, prove gene regulation, or choose a biologically suitable
PWM or threshold. It queries inputs produced or selected elsewhere and records
the operations used on them.

The prioritized implementation gaps and their acceptance targets are tracked
in [`docs/RFC_CISQL_V3.md`](docs/RFC_CISQL_V3.md#work-remaining-after-the-regulatory-foundation).

## Build and run

### Build the interpreter

Cis-QL requires a standard C++11 compiler and `make`:

```bash
make clean
make
```

### Run a query

Execute a `.cql` script using the `cisql` binary:

```bash
./cisql cql_examples/01_promoters.cql
```

Use the `--debug` flag to inspect compilation phases, including token stream, Abstract Syntax Tree (AST), Symbol Table, Intermediate Representation (IR), and execution steps:

```bash
./cisql cql_examples/12_orfs.cql --debug
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

Most examples use the assembly-matched *E. coli* `U00096.3` FASTA and GFF3
files. PWM examples use the plant `MA0054.1` profile from
[JASPAR CORE](https://jaspar.elixir.no/matrix/MA0054.1/) with a short synthetic
FASTA. The matrix demonstrates scanning and p-value syntax; it is not presented
as an *E. coli* model.

GFF3 import preserves source, score, phase, `ID`, `Name`, every `Parent`, and
the complete attribute map. Percent-encoded attribute values are decoded for
queries, while GFF3 export retains the original attribute field. Cis-QL does
not guess which Sequence Ontology types constitute a transcript; use
`EXTRACT FEATURE` with explicit type and attribute filters instead.
Named annotation result sets can be traversed with `EXTRACT CHILDREN OF` for
one downward `Parent` edge or `EXTRACT DESCENDANTS OF` for transitive
relationships. `EXTRACT PARENTS OF` and `EXTRACT ANCESTORS OF` navigate in
the opposite direction.

Validation is an explicit operation, so a partial or imperfect annotation can
still be loaded and queried:

```sql
VALIDATE ANNOTATION annotation AS annotation_report;
EXPORT annotation_report TO "annotation_report.tsv" FORMAT TSV;
```

The report counts records, identities, and parent references; lists unresolved
parents; detects cycles; and flags repeated identities that disagree on
chromosome, feature type, or strand. Repeated IDs with compatible fields are
reported separately as information because GFF3 permits one discontinuous
feature to occupy multiple records. Reports appear in structured JSON, TSV,
and the Validation panel in Cis-QL Studio.

Repeated IDs can be materialized as logical features without treating the
space between their records as part of the feature:

```sql
EXTRACT FEATURE AS coding_parts WHERE TYPE = "CDS";
GROUP coding_parts BY ID AS coding_features;
EXPORT coding_features TO "coding_features.tsv" FORMAT TSV;

EXTRACT MEMBERS OF coding_features AS coding_segments;
EXPORT coding_segments TO "coding_segments.gff3" FORMAT GFF3;
```

Each group retains its original records, their order, phases, scores, parents,
and attributes. Its genomic span and the sum of member lengths are separate
values. Feature groups are not accepted by interval operators; extracting
their members makes the intended interval semantics explicit.

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

EXTRACT GENE AS selected_gene WHERE ID = "gene_1";
EXTRACT CHILDREN OF selected_gene AS direct_components;
EXTRACT DESCENDANTS OF selected_gene AS coding_descendants
    WHERE TYPE = "CDS";

EXTRACT PARENTS OF coding_descendants AS direct_parents;
EXTRACT ANCESTORS OF coding_descendants AS enclosing_genes
    WHERE TYPE = "gene";

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

Hierarchy traversal uses the active annotation dataset and preserved GFF3
relationships. Per-dataset parent and identity indices are built at load time.
`CHILDREN` follows one `Parent` edge;
`DESCENDANTS` follows all reachable edges, excludes the source IDs, preserves
annotation-file order, and terminates even if the input graph contains a
cycle. A source region without `ID` is rejected because it cannot identify a
parent graph node. `PARENTS` and `ANCESTORS` resolve relationships in the
opposite direction; they can start from an annotation record without `ID` if
it declares `Parent`. Unresolved parent identifiers are not inferred from
coordinates. Traversal does not choose a transcript vocabulary or canonical
isoform.

`VALIDATE ANNOTATION` audits the complete named annotation rather than a
derived result set. It reports unresolved parents, incompatible reuse of an
identity, and cycles without silently repairing them. A compatible identity
spread across multiple records is informational, not an error.

`GROUP source BY ID` creates one logical feature per distinct GFF3 identity.
The result is a feature-group collection rather than an interval set, so a
discontinuous CDS, alignment, or transcript cannot accidentally be scanned or
overlapped as one continuous span. `EXTRACT MEMBERS OF` recovers the original
records in group and source order and may be followed by `WHERE`.

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
EXPORT annotation_report TO "annotation_report.tsv" FORMAT TSV;
EXPORT coding_features TO "coding_features.tsv" FORMAT TSV;
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

Statement          ::= LoadStmt | UseStmt | ValidateStmt | GroupStmt | ExportStmt | FindStmt | ExtractStmt
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
ValidateStmt       ::= VALIDATE ANNOTATION ID AS ID SEMICOLON
GroupStmt          ::= GROUP ID BY ID AS ID SEMICOLON
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

ExtractStmt        ::= EXTRACT ExtractSource AliasOpt WhereClause SEMICOLON
ExtractSource      ::= EntityRef | HierarchyRelation OF ID | MEMBERS OF ID
HierarchyRelation  ::= CHILDREN | DESCENDANTS | PARENTS | ANCESTORS

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

## Examples (`cql_examples/`)

Each program focuses on one operation. Most use *E. coli* `U00096.3`; PWM and
track examples use small synthetic fixtures. See
[`cql_examples/README.md`](cql_examples/README.md) for inputs and limits.

| Script | Description | Primary Features |
| :--- | :--- | :--- |
| `01_promoters.cql` | Derive a promoter for `thrA` | `DEFINE PROMOTERS` |
| `02_pwm.cql` | Scan a PWM by relative score | `SCAN`, `THRESHOLD` |
| `03_pvalues.cql` | Filter PWM hits by q-value | `BACKGROUND`, `QVALUE` |
| `04_overlap.cql` | Select motif matches in promoters | `OVERLAPS` |
| `05_count.cql` | Count motif matches per promoter | `COUNT`, `WHERE COUNT` |
| `06_near.cql` | Select motif matches near genes | `NEAR ... WITHIN` |
| `07_modules.cql` | Pair patterns by spacing | `DEFINE MODULE` |
| `08_filters.cql` | Filter annotation fields and attributes | `EXTRACT`, `WHERE` |
| `09_tracks.cql` | Filter narrowPeak fields and metadata | `LOAD TRACK`, `WHERE` |
| `10_track_overlap.cql` | Intersect two tracks | `OVERLAPS` |
| `11_consensus.cql` | Compare two replicate tracks | `CONSENSUS` |
| `12_orfs.cql` | Match complete-codon start-to-stop patterns | regex, `LENGTH MOD 3` |
| `13_gff3.cql` | Validate and traverse GFF3; group repeated IDs | `VALIDATE`, hierarchy queries, `GROUP` |

---

## Cis-QL Studio

Studio edits and runs `.cql` files. It displays sequences, annotations, GC
profiles, validation reports, feature groups, and other result sets written to
`.cisql_results.json`.

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

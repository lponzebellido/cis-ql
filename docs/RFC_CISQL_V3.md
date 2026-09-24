# Cis-QL v3: regulatory genomics direction

Status: incremental implementation. Explicit promoters, calibrated PWM scans,
evidence-preserving interval operations, two-member modules, imported regulatory
tracks, replicate consensus, and modular coordinate filters are implemented.
Later sections are a design contract, not yet accepted syntax.

## Product definition

Cis-QL is a declarative language for expressing and evaluating cis-regulatory
hypotheses. Sequence search and genomic interval operations are its execution
substrate, not its main identity.

The language should make it possible to ask which regulatory evidence supports
a gene or region, how strong that evidence is, and how the result was produced.
It must not equate a motif match or a coexpression edge with direct regulation.

## Non-goals

- Replacing general workflow engines such as Nextflow or Snakemake.
- Becoming a general-purpose DNA editing or sequence utility language.
- Reimplementing every motif, alignment, RNA-seq, or phylogenetics algorithm.
- Reporting regulatory relationships without evidence grades and provenance.

## Work remaining after the regulatory foundation

The next work is ordered by how much it improves real analyses rather than by
how many new statements it adds.

| Priority | Missing capability | Why it matters | Smallest useful acceptance target |
| :--- | :--- | :--- | :--- |
| 1 | Transcript-aware GFF3 selection | Attributes and parentage are now preserved, but transcript graphs and selection policy are not yet first-class | Validate parent graphs; expose explicit `ALL`, canonical, or named-transcript selection |
| 2 | General scalar and oriented-sequence expressions | Regex can describe complex patterns, but programmers need auditable arithmetic, captures, strand-oriented sequence, translation, and genetic-code policy outside the regex itself | Extend the current `MOD`, `START`, `END`, and `STRAND` foundation with reusable expressions and oriented extraction |
| 3 | Enrichment with matched backgrounds | Counts alone cannot distinguish motif enrichment from length or composition effects | Declare foreground/background regions, matching policy, effect size, test, and multiple-testing correction |
| 4 | CRE-to-gene evidence tables | `NEAR` is useful but genomic proximity is only one candidate-linking rule | Import relationship/contact tables and retain typed promoter, distance, contact, expression, and binding evidence per link |
| 5 | General regulatory grammars | Two-site modules cannot express larger heterotypic architectures | Named members, more than two sites, transcript-relative orientation, and grouped aggregation |
| 6 | Dataset identity and catalogs | File provenance alone cannot prevent assembly, annotation-release, or matrix-version mismatches | Assembly identifiers, annotation releases, sequence dictionaries, and catalog-resolved matrix metadata |
| 7 | Comparative regulation | Conservation and motif turnover require orthology and assembly-aware mappings | Import ortholog groups and compare position-aware modules without pretending motif presence proves conserved regulation |

Validation should span several organisms and question types. A bacterial
sequence-pattern fixture tests regex, strand, coordinates, and frame. Yeast can
test compact promoter architecture. A non-MBW plant workflow can test distal
regulation and condition-specific accessibility. A human workflow can test
large annotations, promoter/enhancer evidence, and assembly compatibility.

## Scientific contracts

1. Genomic coordinates are zero-based, half-open internally.
2. Genome assembly, annotation, and imported tracks must be compatible.
3. Promoter boundaries are always explicit; importing a GFF never invents them.
4. Transcript and TSS policy must be visible in the query.
5. Motif hits retain matrix identity, effective background frequencies, motif
   pseudocount, p-value, q-value, and statistical test universe. Database source
   and version become mandatory before catalog-resolved matrices are added.
6. Enrichment requires an explicit or reproducibly generated background.
7. Coexpression, motif presence, accessibility, conservation, and direct
   experimental validation are distinct evidence classes.

## Annotation foundation: preserved GFF3 hierarchy

Implemented syntax:

```cql
EXTRACT FEATURE AS coding_features
  WHERE TYPE = "CDS"
    AND PARENT = "transcript_1"
    AND PHASE = "0"
    AND ATTRIBUTE "protein_id" = "protein_1";
```

GFF3 import retains source, score, phase, `ID`, `Name`, every `Parent`, the
decoded attribute map, and the original attribute field for verbatim
attribute-column re-emission. `ID` and `NAME` are distinct filters. `PARENT`
tests membership rather than comparing the comma-joined field, and `ATTRIBUTE`
permits exact queries on source-specific metadata without turning each key into
a language keyword. JSON and TSV expose the structured annotation evidence.

`FEATURE` means every imported feature type. Cis-QL does not currently decide
that `mRNA`, `transcript`, `lnc_RNA`, or another source vocabulary should be
treated as the canonical transcript concept. Parent-graph validation, child or
descendant traversal, transcript grouping, and explicit transcript-selection
policies remain planned.

## Part 1: explicit TSS-relative promoters

Implemented syntax:

```cql
DEFINE PROMOTERS OF GENE
  FROM TSS
  UPSTREAM 1000 BP
  DOWNSTREAM 200 BP
  AS proximal_promoters;
```

The source can be `GENE`, `TSS`, or an existing result-set alias. The source
features must be stranded. On the positive strand, the TSS anchor is the
zero-based start; on the negative strand it is the half-open end. Intervals are
clamped to the active FASTA chromosome. Both distances are mandatory so a query
cannot silently depend on a biological default.

Transcript-selection policies such as `CANONICAL`, `ALL`, or an explicit
transcript list remain deferred until the language can validate and traverse
the preserved parent graph explicitly.

## Planned language layers

### Part 2A: scoped, evidence-preserving motif scans

Implemented syntax:

```cql
SCAN myb_matrix IN proximal_promoters
  STRAND POSITIVE
  THRESHOLD 85 %
  AS myb_sites;
```

`IN` accepts an annotated biological entity or a named region/motif-hit set.
Coordinates are remapped back to the active genome. Every PWM hit retains the
matrix alias, matrix identifier and name, matrix-file provenance, raw log-odds
score, normalized score percentage, source interval, and a relative start
oriented by the source interval strand. Omitting `IN` deliberately scans the
complete active FASTA, preserving backward compatibility.

Motif-hit aliases are a distinct semantic type and can be exported or reused in
region operations. BED score carries the normalized score on its 0-1000 scale;
GFF3, TSV, JSON, and Cis-QL Studio retain the richer evidence fields. Interval
operations discard hit evidence whenever they change the hit geometry.

### Part 2B1: explicit zero-order background models

Implemented syntax:

```cql
SCAN myb_matrix IN proximal_promoters
  BACKGROUND FROM genome
  THRESHOLD 85 %
  AS myb_sites;

SCAN myb_matrix IN proximal_promoters
  BACKGROUND UNIFORM
  THRESHOLD 85 %
  AS uniform_sites;
```

`BACKGROUND FROM` accepts a loaded sequence dataset, annotated entity, or named
region set. Cis-QL estimates A/C/G/T frequencies after ignoring ambiguous bases
and adds a total pseudocount of 0.1 to prevent zero probabilities. Negative-only
scans complement the frequencies; when both strands are searched, complementary
frequencies are averaged (`A=T`, `C=G`) as in FIMO. The effective frequencies,
source, observed-base count, estimation pseudocount, strand policy, and motif
pseudocount are retained in each hit and its machine-readable exports.

Frequency estimation is incremental and uses constant auxiliary memory, so a
whole plant genome can serve as background without making a second in-memory
copy of its sequences. Region-set backgrounds are counted directly from their
validated intervals in the active genome.

The matrix-to-PSSM conversion now applies the MEME/FIMO convention in which the
total motif pseudocount is distributed according to background frequencies.
Omitting `BACKGROUND` preserves the legacy uniform default, but records it as
such. Because this corrects the former per-letter pseudocount calculation,
percentage thresholds can produce different hit counts than earlier Cis-QL
versions. Future p/q-value syntax will require an explicit statistical policy.

Compatibility basis:

- [FIMO manual](https://meme-suite.org/meme/doc/fimo.html): zero-order
  backgrounds, reverse-complement averaging, motif pseudocounts, and score
  calibration.
- [fasta-get-markov manual](https://meme-suite.org/meme/doc/fasta-get-markov.html):
  zero-order frequency estimation, ambiguity handling, and pseudocount policy.
- [FIMO output contract](https://meme-suite.org/meme/doc/fimo-output-format.html):
  score and future p/q-value meanings.

### Part 2B2: statistically calibrated motif evidence

Implemented for every retained PWM hit:

- a single-window p-value, defined as the probability that a random
  motif-width sequence generated by the selected zero-order background scores
  at least as highly as the observed sequence;
- an exact Benjamini-Hochberg q-value over the complete valid test universe of
  that `SCAN`, including positions below the reporting threshold; and
- the tested-position count, integer score, scale, offset, score range, and
  names of both statistical methods.

The p-value lookup follows FIMO's execution model. Cis-QL maps PSSM cells to
non-negative integer scores using FIMO's range of 1000, then computes
`Pr(score >= x)` by dynamic programming under the same background used for
the log-odds matrix. Windows containing ambiguous bases are neither scored nor
counted as tests. When both strands are selected, each valid position-strand
pair is a test. Separate chromosomes and scoped regions are accumulated before
one multiple-testing correction is applied; overlapping input regions therefore
count as separate tests because they represent separate requested searches.

Unlike an implementation that retains every genomic p-value, Cis-QL keeps a
histogram of integer score levels. This permits exact tied-rank BH correction
for the selected model with memory bounded by motif width and score range,
rather than genome size. The unscaled log-odds score remains available alongside
the discretized score used for calibration.

JSON nests these values under `motifEvidence.statistics`; GFF3 and TSV retain
equivalent fields, and Cis-QL Studio shows p, q, and the test-universe size.
Independent validation exhaustively enumerates the null distribution of a
small motif. The validation suite also compares p-values with the `fimo`
executable when it is locally available.

`THRESHOLD ... %` continues to control relative-score selection, so existing
programs do not silently change semantics.

### Part 2B3: statistical selection in the language

Implemented syntax:

```cql
SCAN myb_matrix IN proximal_promoters
  BACKGROUND FROM genome
  PVALUE <= 1e-4
  AS candidate_myb_sites;

SCAN myb_matrix IN proximal_promoters
  BACKGROUND FROM genome
  QVALUE <= 0.05
  AS significant_myb_sites;
```

Both `<` and `<=` are supported and probability thresholds must lie in `[0, 1]`.
Scientific notation is tokenized as one numeric literal. A `SCAN` may select by
p-value or q-value, but not both. When no relative `THRESHOLD` is written, a
statistical scan evaluates the entire score range instead of inheriting the
legacy 75% default. If both forms are explicit, they are combined with logical
AND.

The statistical universe is accumulated before final selection. To avoid
materializing every genomic window, Cis-QL retains only candidates with
`p <= alpha` while building the complete observed score histogram. This is safe
for q-value filtering because a Benjamini-Hochberg adjusted p-value cannot be
smaller than its raw p-value. Strict-boundary filtering is applied after q-values
have been assigned.

Still planned:

- richer JASPAR metadata such as TF family and matrix release/version;
- installed-FIMO parity fixtures in continuous integration.

### Part 2C: explicit modular coordinate constraints

Implemented syntax:

```cql
FIND MOTIF "ATG(?:(?!TAA|TAG|TGA)[ACGT]{3})*(?:TAA|TAG|TGA)"
  STRAND POSITIVE AS start_stop_candidates
  WHERE LENGTH MOD 3 = 0 AND START MOD 3 = 0;
```

Numeric conditions may apply `MOD` before their relational comparison.
`START` and `END` expose zero-based, half-open genomic coordinates, and
`STRAND` compares with `"+"`, `"-"`, or `"."`. These are general language
properties rather than an ORF-specific operator. A codon-stepping regex is
still necessary when the search engine must skip an out-of-frame stop and keep
looking for the next in-frame stop; filtering a lazy arbitrary-length regex
after it has already chosen a match cannot change that choice.

The current result sequence is stored in reference orientation. Strand-oriented
extraction, capture groups, translation, genetic-code selection, and general
scalar expressions remain planned.

### Part 3: regulatory interval algebra

#### Part 3A: directional overlap selection

Implemented syntax:

```cql
OVERLAPS significant_myb_sites WITH proximal_promoters
  AS promoter_supported_myb_sites;
```

This is an interval semi-join rather than a geometric intersection. It emits
each interval from the left/query set at most once when any interval in the
right/reference set overlaps it. The query interval is not clipped, so its
coordinates, sequence, and motif evidence remain auditable. Every matching
reference is recorded in an `overlapEvidence` array with its source-set alias,
coordinates, identity, and imported `trackEvidence` when present. This allows
an accessibility interval to retain independent binding support without
overwriting either observation. Intervals use the
same zero-based, half-open contract as the rest of the engine; touching
boundaries do not overlap. `INTERSECT` remains available when overlap segments
themselves are the intended result.

The reference side is indexed per chromosome with an interval tree. Matching
references are enumerated in deterministic coordinate order, so the semi-join
remains output-sensitive and emits the query once while retaining all of its
supporting records.

#### Part 3B: bounded nearest-reference selection

Implemented syntax:

```cql
NEAR significant_myb_sites TO GENE WITHIN 2 KB
  AS proximal_gene_candidates;
```

`NEAR query TO reference WITHIN distance` is a directional, bounded nearest
selection. For every valid query interval, Cis-QL finds one nearest reference
on the same chromosome and retains the complete query only when its interval
gap is less than or equal to the explicit limit. Overlaps and directly adjacent
half-open intervals both have gap zero; the stored `overlaps` Boolean keeps
those cases distinguishable. Limits must be finite, non-negative, explicitly
unit-qualified, and resolve to whole base pairs.

Each retained record preserves existing motif evidence and attaches a typed
`spatialRelation` containing the reference set, chosen reference coordinates,
strand, type, name, observed distance, requested maximum distance, and overlap
state. JSON preserves the nested structure; TSV and GFF3 expose equivalent
fields. Reference coordinates embedded in GFF3 attributes remain zero-based,
half-open even though the feature columns follow the GFF3 coordinate system.
BED remains intentionally lossy.

The choice is deterministic: equal-distance references are ordered by start,
end, name, type, then strand. A per-chromosome index of sorted starts and prefix
maximum ends avoids scanning all reference intervals for every query.

This operator creates a proximity-based candidate association. It must not be
reported as proof that the selected feature regulates the recorded gene; later
evidence-integration layers can add accessibility, binding, expression, or
chromatin-contact support.

#### Part 3C: overlap counts per regulatory region

Implemented syntax:

```cql
COUNT significant_myb_sites IN proximal_promoters
  AS promoter_myb_counts;

COUNT significant_myb_sites IN proximal_promoters
  AS supported_promoters
  WHERE COUNT >= 1;
```

The second entity is the container set and therefore determines the emitted
geometry. Every valid container is emitted once, including containers with a
zero count. `countEvidence` records the `OVERLAPS` relation, counted-set alias,
container-set alias, and non-negative integer count. It is preserved by later
selection operations while the container geometry remains unchanged and is
available in JSON, GFF3, TSV, and both Studio inspectors. BED is lossy.

Half-open overlap semantics match `OVERLAPS`, so touching boundaries are not
counted. Per chromosome, separately sorted start and end arrays reduce each
container count to two binary searches instead of comparing every interval
pair. `WHERE COUNT` supports integer thresholds without genomic units.
Counting is strand-agnostic and counts input records, so duplicate intervals
remain separate observations; callers must constrain or normalize the input
set when a distinct-locus interpretation is required.

This is descriptive aggregation. A larger raw count is not automatically motif
enrichment or stronger regulatory evidence because expected counts vary with
region length, nucleotide composition, accessibility, motif model, and the
selected statistical background.

#### Part 3D: constrained two-member cis-regulatory modules

Implemented syntax:

```cql
DEFINE MODULE
  FROM myb_sites WITH bhlh_sites
  SPACING 5 BP TO 30 BP
  ORDER AS_WRITTEN
  ORIENTATION OPPOSITE
  AS myb_bhlh_modules;
```

Both inputs are named region or motif-hit sets. Members must lie on the same
chromosome. `SPACING` is the edge-to-edge gap between their zero-based,
half-open intervals; its explicit lower and upper bounds are inclusive,
non-negative, finite, unit-qualified, and must resolve to whole base pairs.
Overlapping and directly adjacent intervals both have gap zero, while
`observedOrder` distinguishes an overlap from a separated pair.

`ORDER AS_WRITTEN` requires the first-set member to have a lower start
coordinate than the second-set member. Tied starts do not satisfy that policy.
`ORDER ANY` accepts either order and
records `FIRST_BEFORE_SECOND`, `SECOND_BEFORE_FIRST`, or `OVERLAPPING`.
Ordering is deliberately reference-relative, not gene- or transcript-relative.

`ORIENTATION SAME` and `OPPOSITE` compare the members' reference strands and
require both to be `+` or `-`. `ORIENTATION ANY` also accepts unstranded members
and records `UNKNOWN` when their relationship cannot be determined. When the
same alias appears on both sides, Cis-QL produces one canonical unordered pair,
does not pair a record with itself, and does not treat opposite-strand records
at identical coordinates as two independent sites.

Each result spans both members, has type `cis_regulatory_module`, and stores a
typed `moduleEvidence` object containing the requested constraints, observed
spacing/order/orientation, source-set names, member coordinates and metadata,
and each member's PWM evidence when present. JSON retains the nested structure;
GFF3 and TSV expose the constraints, member identity, matrix ID, score, p-value,
and q-value, and both Studio viewers inspect the nested result.
BED carries only the outer span and is intentionally lossy. Spatial and count
operations can be applied to a module without discarding its member evidence;
geometric clipping or merging clears evidence whose span is no longer intact.

The candidate search uses chromosome-local, start-sorted indexes and the
maximum member width to bound the candidate range. Runtime is proportional to
index construction plus the candidates within the requested spatial reach,
rather than the Cartesian product for typical motif-width inputs.

This operator evaluates a declared cis grammar. It does not infer spacing
bounds, establish TF cooperativity, or convert motif co-occurrence into direct
regulatory evidence. Real analyses must justify the grammar from prior
evidence, a benchmark, or a documented sensitivity analysis.

Still planned:

- unbounded `CLOSEST`, standalone `DISTANCE`, and richer grouped aggregation;
- Matched backgrounds and enrichment with multiple-testing correction.
- Modules with more than two members and transcript-relative orientation.

### Part 4: evidence integration

#### Part 4A: imported regulatory tracks

Implemented syntax:

```cql
LOAD TRACK "sample_accessibility.narrowPeak"
  FORMAT NARROWPEAK
  EVIDENCE ACCESSIBILITY
  ASSAY "ATAC-seq"
  SAMPLE "pigmented petal"
  CONDITION "pigmented"
  REPLICATE "A1"
  CONTROL "input"
  AS accessibility_peaks;
```

`LOAD TRACK` accepts BED and narrowPeak files as named interval sets without
replacing the active GFF3 annotation. Coordinates remain zero-based and
half-open. BED name, score, and strand are retained; narrowPeak additionally
retains signalValue, the supplied `-log10(p)` and `-log10(q)` values, and the
summit offset and absolute summit position. The file path, declared format,
query alias, required evidence class (`ACCESSIBILITY`, `BINDING`, or `OTHER`),
and optional assay, sample, condition, replicate, and control labels are
recorded as typed `trackEvidence` in JSON, GFF3, TSV, and Studio. Metadata
clauses may appear in any order; duplicates are rejected. The class is an
explicit user declaration, not an inference from the filename or an assertion
that the experiment is valid.
When `OVERLAPS` combines two tracks, the left/query observation stays in
`trackEvidence` and all matching right/reference observations are stored in
`overlapEvidence`. JSON and Studio expose the structured records directly;
TSV stores the array as JSON in `overlap_evidence_json`, and GFF3 stores a
count plus a percent-encoded JSON attribute. If a reference already has
overlap evidence, the prior observations remain nested under
`supportingEvidence`. This records the real operation path and avoids
misrepresenting transitive support as a direct overlap.

Track evidence is executable in `WHERE` filters:

```cql
EXTRACT accessibility_peaks AS strong_accessibility
  WHERE TRACK_SCORE >= 600
    AND SIGNAL_VALUE >= 10
    AND EVIDENCE_CLASS = "ACCESSIBILITY"
    AND CONDITION = "pigmented"
    AND REPLICATE = "A1";
```

`TRACK_SCORE`, `SIGNAL_VALUE`, `MINUS_LOG10_PVALUE`, and
`MINUS_LOG10_QVALUE` use finite non-negative thresholds without genomic or
percentage units. Missing optional narrowPeak values fail the corresponding
condition. `EVIDENCE_CLASS`, `ASSAY`, `SAMPLE`, `CONDITION`, `REPLICATE`, and
`CONTROL` use exact string equality.
The properties address the region's primary `trackEvidence`. To constrain a
reference track unambiguously when multiple peaks may overlap one query,
filter that track first and pass the filtered alias to `OVERLAPS`.
Filtering does not normalize or calibrate imported measurements: narrowPeak
score and `signalValue` remain defined by the upstream caller and experiment.
Thresholds therefore require assay-specific justification.

When a FASTA is active, every track chromosome and interval bound is checked
against it and interval sequence is attached. A track loaded before its genome
cannot be assembly-validated, so real workflows should load the sequence
first. Imported tracks can participate in `SCAN`, `OVERLAPS`, `COUNT`, `NEAR`,
and module construction. The language preserves upstream statistics; it does
not infer assay type, quality, or biological validity from a filename.
Using `OVERLAPS` between replicate peak sets reports coordinate-level
concordance only. It is not an IDR calculation, does not account for replicate
quality or control design, and is not a formal reproducibility test.

Anchor-preserving coordinate consensus is implemented explicitly:

```cql
CONSENSUS FROM [binding_rep1, binding_rep2]
  ANCHOR binding_rep1
  MIN_SUPPORT 2
  MIN_RECIPROCAL_OVERLAP 50 %
  MAX_SUMMIT_DISTANCE 20 BP
  AS reproducible_binding
  WHERE SUPPORT_COUNT >= 2;
```

The anchor must be one unique member of a list containing at least two named
region sets. It supplies output coordinates and primary evidence and counts as
one supporting set. Every other input contributes at most one support unit per
anchor interval, while all of its overlapping observations remain auditable.
`MIN_RECIPROCAL_OVERLAP` can require each accepted pair to cover a declared
fraction of both intervals. `MAX_SUMMIT_DISTANCE` can additionally require two
narrowPeak summits within an explicit distance; records without summit data do
not satisfy that criterion. `consensusEvidence` records the anchor, input
aliases, thresholds, minimum support, and observed support and remains attached
when the consensus is later used as an `OVERLAPS` reference. `SUPPORT_COUNT`
exposes observed support to `WHERE`. This operation deliberately does not merge
peak geometry, infer replicate quality, or claim IDR-equivalent reproducibility.

Still planned:

- explicit genome-assembly identity and richer experimental design schemas;
- tabular expression/coexpression inputs;
- Accessibility, DAP/ChIP, expression, coexpression, and literature evidence.
- CRE-to-gene linking by promoter, distance, or an imported relationship.
- Evidence tables and auditable candidate ranking.

### Part 5: comparative regulation

- Ortholog groups and genome-version-aware identifiers.
- Conserved modules, motif turnover, and position-aware comparisons.
- Import CoExp/CoExpPhylo results instead of duplicating their pipelines.

## Cross-organism validation portfolio

No single pathway should define the language or its examples. Reference
workflows should be small enough to inspect, use versioned public inputs, and
test a distinct scientific contract:

- *E. coli*: strand-aware sequence patterns, coordinate phase, and candidate
  start-to-stop regions, without presenting regex hits as gene predictions.
- *Saccharomyces cerevisiae*: compact promoter architecture with sourced motif
  models and an independently checkable expected set.
- *Arabidopsis thaliana*: a stress-, hormone-, or light-response workflow using
  annotation, accessibility, and TF evidence rather than an MBW showcase.
- Human: assembly-matched GENCODE-derived annotation and ENCODE-derived tracks
  to exercise transcript policy, large interval sets, and multi-source
  promoter/enhancer support.

Each workflow must distinguish deterministic language tests from biological
validation, document licenses and releases, and state which outputs are direct
observations, coordinate associations, or predictions.

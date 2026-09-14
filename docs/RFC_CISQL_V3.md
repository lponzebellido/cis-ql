# Cis-QL v3: regulatory genomics direction

Status: incremental implementation. Parts 1 through 2B3 and Parts 3A-3C are
implemented; later parts are a design contract, not yet accepted syntax.

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
transcript list are intentionally deferred until the GFF data model preserves
parent/child relationships and all attributes.

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
coordinates, sequence, and motif evidence remain auditable. Intervals use the
same zero-based, half-open contract as the rest of the engine; touching
boundaries do not overlap. `INTERSECT` remains available when overlap segments
themselves are the intended result.

The reference side is indexed per chromosome by sorted starts and prefix
maximum ends. After index construction, each query is evaluated in logarithmic
time without enumerating every overlapping pair.

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

Still planned:

- unbounded `CLOSEST`, standalone `DISTANCE`, and richer grouped aggregation;
- Motif modules with order, orientation, minimum and maximum spacing.
- Matched backgrounds and enrichment with multiple-testing correction.

### Part 4: evidence integration

- BED/narrowPeak and tabular expression/coexpression inputs.
- Accessibility, DAP/ChIP, expression, coexpression, and literature evidence.
- CRE-to-gene linking by promoter, distance, or an imported relationship.
- Evidence tables and auditable candidate ranking.

### Part 5: comparative regulation

- Ortholog groups and genome-version-aware identifiers.
- Conserved modules, motif turnover, and position-aware comparisons.
- Import CoExp/CoExpPhylo results instead of duplicating their pipelines.

## Target anthocyanin use case

The first biological benchmark should use known maize anthocyanin genes and
transport machinery. A direct case such as `ZmMRP3` is a positive control;
`Bz2`, `Wrky33`, transporter families, and ART1-like candidates exercise
different evidence levels. A successful query must reproduce known evidence,
rank plausible new candidates, and clearly identify which relationships remain
predictions.

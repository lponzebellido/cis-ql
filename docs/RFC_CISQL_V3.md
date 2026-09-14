# Cis-QL v3: regulatory genomics direction

Status: incremental implementation. Part 1, Part 2A, and Part 2B1 are
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
   pseudocount, and their available statistical evidence. P-value, q-value,
   database source, and version become mandatory once Part 2B2 introduces
   statistically calibrated scans.
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

### Part 2B2: statistically calibrated motif scans

- JASPAR identifiers, TF names/families, source and version.
- P-values and multiple-testing-corrected q-values.
- A FIMO-compatible execution backend and parity tests.

### Part 3: regulatory interval algebra

- `OVERLAPS`, `NEAR`, `CLOSEST`, `DISTANCE`, `COUNT`, and grouped aggregation.
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

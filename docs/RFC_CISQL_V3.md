# Cis-QL v3: regulatory genomics direction

Status: incremental implementation. Part 1 is implemented; later parts are a
design contract, not yet accepted syntax.

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
5. Motif hits retain matrix identity, score, background model, p-value, q-value,
   database source, and version.
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

### Part 2: evidence-preserving motif scans

- Typed motif sets and motif hits.
- JASPAR identifiers, TF names/families, source and version.
- Raw/log-odds score, configurable background, p-value and q-value.
- Scan a named region set rather than only the complete active genome.
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

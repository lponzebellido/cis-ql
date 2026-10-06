# GRQL v3: regulatory genomics direction

Status: incremental implementation. Explicit promoters, transcript selection,
strand-aware sequence predicates, calibrated PWM scans, evidence-preserving
interval operations, two-member modules, imported regulatory tracks, replicate
consensus, modular coordinate filters, fixed-index sequence slicing, numbered
regular-expression captures, and explicit interval translation are
implemented. Later sections are a design contract, not yet accepted syntax.

## Product definition

GRQL is a declarative language for expressing and evaluating cis-regulatory
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
| 1 | Richer scalar expressions and CDS translation policy | Fixed-index slicing, numbered regex captures, and interval translation are explicit, but programmers still need reusable values, joined CDS assembly, phase handling, and an explicit initiation policy | Add reusable string/numeric expressions; specify a separate CDS mode for joined records and initiator codons |
| 2 | Enrichment with matched backgrounds | Counts alone cannot distinguish motif enrichment from length or composition effects | Declare foreground/background regions, matching policy, effect size, test, and multiple-testing correction |
| 3 | CRE-to-gene evidence tables | `NEAR` is useful but genomic proximity is only one candidate-linking rule | Import relationship/contact tables and retain typed promoter, distance, contact, expression, and binding evidence per link |
| 4 | General regulatory grammars | Two-site modules cannot express larger heterotypic architectures | Named members, more than two sites, transcript-relative orientation, and grouped aggregation |
| 5 | Dataset identity and catalogs | File provenance alone cannot prevent assembly, annotation-release, or matrix-version mismatches | Assembly identifiers, annotation releases, sequence dictionaries, and catalog-resolved matrix metadata |
| 6 | Comparative regulation | Conservation and motif turnover require orthology and assembly-aware mappings | Import ortholog groups and compare position-aware modules without pretending motif presence proves conserved regulation |

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

```grql
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

EXTRACT TRANSCRIPTS OF selected_gene
  TYPES ["mRNA", "transcript"]
  SELECT ALL
  AS transcripts;

EXTRACT TRANSCRIPTS OF selected_gene
  TYPES ["mRNA", "transcript"]
  SELECT ID "transcript_1"
  AS named_transcript;

EXTRACT TRANSCRIPTS OF selected_gene
  TYPES ["mRNA", "transcript"]
  SELECT ATTRIBUTE "tag" = "canonical"
  AS source_labeled_transcript;

VALIDATE ANNOTATION annotation AS hierarchy_report;
EXPORT hierarchy_report TO "hierarchy_report.tsv" FORMAT TSV;

GROUP coding_parts BY ID AS coding_features;
EXTRACT MEMBERS OF coding_features AS coding_segments;
```

GFF3 import retains source, score, phase, `ID`, `Name`, every `Parent`, the
decoded attribute map, and the original attribute field for verbatim
attribute-column re-emission. `ID` and `NAME` are distinct filters. `PARENT`
tests membership rather than comparing the comma-joined field, and `ATTRIBUTE`
permits exact queries on source-specific metadata without turning each key into
a language keyword. JSON and TSV expose the structured annotation evidence.

`CHILDREN` follows one `Parent` edge from every identified source region.
`DESCENDANTS` follows the relation transitively, excludes the source IDs,
preserves annotation-file order, and cannot loop indefinitely on cyclic input.
Downward sources without GFF3 identity are rejected rather than matched by
coordinate.
`PARENTS` and `ANCESTORS` traverse in the opposite direction. They require
annotation evidence but can start from a record without `ID` when it declares
`Parent`. Unresolved identifiers are not replaced by coordinate-based guesses.
Each loaded annotation has reusable parent and identity indices, so traversal
does not rescan the full GFF3 for every query.

`VALIDATE ANNOTATION` produces a typed report over the complete named dataset.
It identifies unresolved `Parent` references, cycles, and repeated IDs whose
records disagree on chromosome, feature type, or strand. Compatible repeated
IDs are listed separately because GFF3 can represent one discontinuous feature
with multiple lines. Validation is explicit rather than a load-time failure so
partial annotations remain usable. JSON retains structured diagnostics, TSV
provides summary and diagnostic rows, and Studio displays the same report.

`GROUP source BY ID` materializes one typed logical feature for each distinct
identity in an annotation dataset or annotation-backed result set. The group
stores the genomic span and summed member length as different quantities and
retains every original record, including member-specific score and phase.
Groups are not intervals and therefore cannot enter interval operators until
`EXTRACT MEMBERS OF` makes the member records explicit. Grouping fails on
missing IDs or incompatible chromosome, type, or strand instead of discarding
records or inventing a continuous feature.
This follows the [Sequence Ontology GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md),
which defines repeated IDs as the records that collectively represent one
discontinuous feature.

`FEATURE` means every imported feature type. `EXTRACT TRANSCRIPTS` makes the
source-specific transcript vocabulary and selection rule part of the query.
It follows the full descendant graph, even through intermediate records whose
types are not listed, then applies exact `TYPES` and `SELECT` matches. `ALL`
keeps every matching record, `ID` names a GFF3 identity, and `ATTRIBUTE` names
a source key and value. GRQL does not treat `mRNA`, `transcript`, `lnc_RNA`,
or an attribute named `canonical` as a universal convention.

## Part 1: explicit TSS-relative promoters

Implemented syntax:

```grql
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

An `EXTRACT TRANSCRIPTS` result can be used as the promoter source. The query
therefore records whether promoters came from every declared transcript, one
named transcript, or records carrying an exact source attribute. GRQL does
not provide a `CANONICAL` keyword because that designation is annotation-source
specific.

## Planned language layers

### Part 2A: scoped, evidence-preserving motif scans

Implemented syntax:

```grql
SCAN tf_matrix IN proximal_promoters
  STRAND POSITIVE
  THRESHOLD 85 %
  AS motif_sites;
```

`IN` accepts an annotated biological entity or a named region/motif-hit set.
Coordinates are remapped back to the active genome. Every PWM hit retains the
matrix alias, matrix identifier and name, matrix-file provenance, raw log-odds
score, normalized score percentage, source interval, and a relative start
oriented by the source interval strand. Omitting `IN` deliberately scans the
complete active FASTA, preserving backward compatibility.

Motif-hit aliases are a distinct semantic type and can be exported or reused in
region operations. BED score carries the normalized score on its 0-1000 scale;
GFF3, TSV, JSON, and GRQL Studio retain the richer evidence fields. Interval
operations discard hit evidence whenever they change the hit geometry.

### Part 2B1: explicit zero-order background models

Implemented syntax:

```grql
SCAN tf_matrix IN proximal_promoters
  BACKGROUND FROM genome
  THRESHOLD 85 %
  AS motif_sites;

SCAN tf_matrix IN proximal_promoters
  BACKGROUND UNIFORM
  THRESHOLD 85 %
  AS uniform_sites;
```

`BACKGROUND FROM` accepts a loaded sequence dataset, annotated entity, or named
region set. GRQL estimates A/C/G/T frequencies after ignoring ambiguous bases
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
percentage thresholds can produce different hit counts than earlier GRQL
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

The p-value lookup follows FIMO's execution model. GRQL maps PSSM cells to
non-negative integer scores using FIMO's range of 1000, then computes
`Pr(score >= x)` by dynamic programming under the same background used for
the log-odds matrix. Windows containing ambiguous bases are neither scored nor
counted as tests. When both strands are selected, each valid position-strand
pair is a test. Separate chromosomes and scoped regions are accumulated before
one multiple-testing correction is applied; overlapping input regions therefore
count as separate tests because they represent separate requested searches.

Unlike an implementation that retains every genomic p-value, GRQL keeps a
histogram of integer score levels. This permits exact tied-rank BH correction
for the selected model with memory bounded by motif width and score range,
rather than genome size. The unscaled log-odds score remains available alongside
the discretized score used for calibration.

JSON nests these values under `motifEvidence.statistics`; GFF3 and TSV retain
equivalent fields, and GRQL Studio shows p, q, and the test-universe size.
Independent validation exhaustively enumerates the null distribution of a
small motif. The validation suite also compares p-values with the `fimo`
executable when it is locally available.

`THRESHOLD ... %` continues to control relative-score selection, so existing
programs do not silently change semantics.

### Part 2B3: statistical selection in the language

Implemented syntax:

```grql
SCAN tf_matrix IN proximal_promoters
  BACKGROUND FROM genome
  PVALUE <= 1e-4
  AS pvalue_sites;

SCAN tf_matrix IN proximal_promoters
  BACKGROUND FROM genome
  QVALUE <= 0.05
  AS qvalue_sites;
```

Both `<` and `<=` are supported and probability thresholds must lie in `[0, 1]`.
Scientific notation is tokenized as one numeric literal. A `SCAN` may select by
p-value or q-value, but not both. When no relative `THRESHOLD` is written, a
statistical scan evaluates the entire score range instead of inheriting the
legacy 75% default. If both forms are explicit, they are combined with logical
AND.

The statistical universe is accumulated before final selection. To avoid
materializing every genomic window, GRQL retains only candidates with
`p <= alpha` while building the complete observed score histogram. This is safe
for q-value filtering because a Benjamini-Hochberg adjusted p-value cannot be
smaller than its raw p-value. Strict-boundary filtering is applied after q-values
have been assigned.

Still planned:

- richer JASPAR metadata such as TF family and matrix release/version;
- installed-FIMO parity fixtures in continuous integration.

### Part 2C1: explicit modular coordinate constraints

Implemented syntax:

```grql
FIND MOTIF "ATG(?:(?!TAA|TAG|TGA)[ACGT]{3})*(?:TAA|TAG|TGA)"
  STRAND POSITIVE AS start_stop_candidates
  WHERE LENGTH MOD 3 = 0 AND START MOD 3 = 0;

FIND MOTIF "ATGNNNTAA" AS bounded_patterns
  WHERE ORIENTED_SEQUENCE STARTS_WITH "ATG"
    AND ORIENTED_SEQUENCE ENDS_WITH "TAA";

EXTRACT bounded_patterns AS alanine_second_codon
  WHERE SLICE ORIENTED_SEQUENCE FROM 3 TO 6 = "GCT";

EXTRACT start_stop_candidates AS terminal_taa
  WHERE SLICE ORIENTED_SEQUENCE FROM -3 TO END = "TAA";

FIND MOTIF "ATG([ACGT]{3})(TAA|TAG|TGA)" AS codon_patterns;
EXTRACT codon_patterns AS tga_middle
  WHERE CAPTURE 1 = "TGA" AND CAPTURE 2 = "TAA";
```

Numeric conditions may apply `MOD` before their relational comparison.
`START` and `END` expose zero-based, half-open genomic coordinates, and
`STRAND` compares with `"+"`, `"-"`, or `"."`. These are general language
properties rather than an ORF-specific operator. A codon-stepping regex is
still necessary when the search engine must skip an out-of-frame stop and keep
looking for the next in-frame stop; filtering a lazy arbitrary-length regex
after it has already chosen a match cannot change that choice.

Result sequences remain stored in reference orientation. `SEQUENCE` queries
that representation, while `ORIENTED_SEQUENCE` reverse-complements a
negative-strand interval for comparison without changing its coordinates or
stored evidence. Both accept literal `=`, `STARTS_WITH`, `ENDS_WITH`, and
`CONTAINS` predicates. Unstranded records do not satisfy
`ORIENTED_SEQUENCE`.

`SLICE property FROM start TO end` applies half-open offsets to `SEQUENCE`,
`ORIENTED_SEQUENCE`, or `PROTEIN_SEQUENCE` before the text comparison.
Non-negative indices count from the start, negative indices count from the end,
and `END` denotes the string length. Thus `FROM -3 TO END` selects the final
three characters regardless of sequence length. An empty, reversed, or
out-of-bounds range evaluates as no match instead of being clamped, so a
program cannot silently compare a shorter string. Slicing composes with the
existing Boolean condition tree and works on motif hits and named region sets.
Materialized sequence expressions and reusable scalar values remain planned.

`CAPTURE n` addresses the nth explicit capturing group in a `FIND MOTIF`
regular expression, beginning at one. The engine's lookahead group is not part
of public numbering. Every match retains the original pattern plus each
group's matched state, value, and zero-based half-open offsets relative to the
oriented full match. Negative-strand captures therefore use the 5'-to-3'
matched text rather than reference-left coordinates. Optional unmatched groups
remain explicit and do not satisfy text predicates; a group index absent from
the pattern also evaluates as no match. JSON preserves typed capture objects,
while TSV and GFF3 carry the same structure as JSON. Geometry-preserving
selection keeps the evidence, whereas clipping and merging discard it.

Named captures, capture slicing, materialized expressions, and reusable scalar
values remain planned.

### Part 2C2: explicit interval translation

Implemented syntax:

```grql
TRANSLATE candidates CODE 1 FRAME 0 AS standard;
TRANSLATE candidates CODE 4 FRAME 0 AS table4;

EXTRACT table4 AS tga_trp
  WHERE PROTEIN_SEQUENCE = "MW*";
```

The source must be a region or motif-hit alias and an active FASTA must be
available. Translation uses each interval's oriented 5'-to-3' sequence, so a
negative-strand interval is reverse-complemented and an unstranded interval is
rejected. `FRAME` is an offset of 0, 1, or 2 from the first oriented base.
Trailing bases that do not form a complete codon are ignored, ambiguous codons
produce `X`, and stop codons produce `*`.

`CODE` accepts the translation-table identifiers currently published in the
[NCBI genetic-code catalog](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/taxonomy-processing/genetic-codes/):
1–6, 9–16, and 21–33. Translation uses the table's amino-acid assignments but
does not apply its separate initiator row. A general interval cannot be assumed
to begin at a biological translation start; applying alternative initiators
there would silently turn a positional query into a CDS policy.

The result is a regular region set with unchanged coordinates and reference
DNA plus typed `translationEvidence`: genetic code, frame, and protein string.
`PROTEIN_SEQUENCE` makes that protein available to case-insensitive literal
`=`, `STARTS_WITH`, `ENDS_WITH`, and `CONTAINS` predicates. JSON, TSV, GFF3,
and Studio retain the evidence. `INTERSECT`, `UNION`, and `EXCEPT` invalidate
it when they alter interval geometry.

This is not yet transcript translation. A later CDS-specific mode must define
how discontinuous members are ordered and concatenated, how GFF3 phase is
applied on both strands, which initiator row is used, and whether terminal stop
or completeness checks are required. Those choices should remain visible in
the program instead of being inferred from feature names.

### Part 3: regulatory interval algebra

#### Part 3A: directional overlap selection

Implemented syntax:

```grql
OVERLAPS qvalue_sites WITH proximal_promoters
  AS promoter_sites;
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

```grql
NEAR qvalue_sites TO GENE WITHIN 2 KB
  AS nearby_gene_links;
```

`NEAR query TO reference WITHIN distance` is a directional, bounded nearest
selection. For every valid query interval, GRQL finds one nearest reference
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

```grql
COUNT qvalue_sites IN proximal_promoters
  AS promoter_site_counts;

COUNT qvalue_sites IN proximal_promoters
  AS nonempty_promoters
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

```grql
DEFINE MODULE
  FROM factor_a_sites WITH factor_b_sites
  SPACING 5 BP TO 30 BP
  ORDER AS_WRITTEN
  ORIENTATION OPPOSITE
  AS paired_sites;
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
same alias appears on both sides, GRQL produces one canonical unordered pair,
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

```grql
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

```grql
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

```grql
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

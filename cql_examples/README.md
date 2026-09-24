# Cis-QL regulatory analysis examples

These eleven runnable programs form a progressive analysis of one synthetic
anthocyanin-regulatory locus. The locus and the MYB model make the examples
concrete; neither is built into Cis-QL. The sequence of operations is the
reusable part: define the regulatory search space, score a TF model, calibrate
and filter its matches, relate them to annotations and experimental tracks,
then retain the evidence used to form each candidate.

## Adapting the workflow

To use the same programs for another regulatory question, replace each input
according to its role rather than copying the anthocyanin interpretation:

| Role in the query | Anthocyanin fixture | Replacement in another study |
| :--- | :--- | :--- |
| Reference sequence | synthetic locus FASTA | the relevant genome, contigs, or locus in the correct assembly |
| Feature coordinates | synthetic genes and enhancer in GFF3 | trusted annotations or explicitly defined candidate regions |
| Binding model | MYB.Ph3 PWM | a justified PWM for the TF or TF family being tested |
| Experimental support | synthetic accessibility and binding peaks | compatible BED/narrowPeak outputs from the assays and conditions of interest |
| Regulatory assumptions | promoter bounds, distance, spacing, support thresholds | values justified for the organism, regulatory system, and question |

For example, the same pattern can test a stress-response TF near induced
genes, a developmental factor in accessible enhancers, or replicate-supported
binding around a microbial promoter. Cis-QL evaluates the declared coordinate
and motif relationships; it does not supply the biological assumptions or
make those examples equivalent across systems.

## Progressive reference workflow

1. `01_define_promoters.cql` derives three strand-aware promoter windows from
   explicit TSS-relative bounds.
2. `02_score_myb_sites.cql` scans a sourced plant MYB matrix with a relative
   score threshold. This helps inspect model scores but is not yet a
   statistically calibrated selection.
3. `03_calibrated_myb_sites.cql` estimates the background from the loaded
   genome and applies `QVALUE <= 0.01`. Four sites remain.
4. `04_promoter_supported_sites.cql` retains the two significant sites that
   overlap a candidate promoter without clipping their coordinates or losing
   PWM evidence.
5. `05_count_promoter_support.cql` counts promoter-scoped sites in all three
   candidate promoters, retaining counts `[1, 1, 0]`, then selects the two
   promoters with support.
6. `06_nearest_gene_candidates.cql` retains two significant sites within 50 bp
   of a gene and records candidate links at interval distances of 41 bp and
   30 bp.
7. `07_enhancer_myb_module.cql` selects the two genome-calibrated sites inside
   the annotated candidate enhancer and pairs them as one homotypic module.
   Their half-open intervals are 400-409 and 420-429, so the observed
   edge-to-edge spacing is 11 bp. Both have the same observed reference
   orientation, but the query uses `ORIENTATION ANY` because the fixture gives
   no biological reason to require that arrangement.
8. `08_integrated_anthocyanin_query.cql` connects the complete path. It keeps
   promoter evidence separate, summarizes it per promoter, builds the enhancer
   module, and records its nearest gene as a deliberately provisional target
   hypothesis.
9. `09_accessible_myb_evidence.cql` imports a synthetic narrowPeak
   accessibility track, counts calibrated MYB sites inside every peak, retains
   peaks with motif support, and records their nearest gene as a candidate
   association. The emitted records preserve both quantitative track evidence
   and the overlap count.
10. `10_accessible_bound_myb_candidates.cql` retains accessible peaks where a
    separately declared, quantitatively filtered MYB-binding track overlaps a
    filtered accessibility track, counts calibrated motif support, and forms
    bounded nearest-gene hypotheses. The accessibility record remains the
    primary `trackEvidence`; every matching binding record is retained in
    `overlapEvidence`.
11. `11_replicate_supported_candidates.cql` attaches condition, replicate,
    and control labels to accessibility and binding tracks, filters each
    replicate explicitly, and uses an anchor-preserving `CONSENSUS` requiring
    two distinct supporting sets, at least 50% reciprocal interval overlap,
    and summits no farther than 5 bp apart before combining binding with
    accessibility and motif support. The result records all three criteria;
    downstream overlap evidence retains the consensus rule and the R2
    observation.

The shared `anthocyanin_regulatory_demo` FASTA and GFF3 files contain three
annotated genes, one candidate enhancer, two promoter-local MYB instances, and
two enhancer-local instances. These deterministic fixtures exercise query
semantics. They do not assert that this sequence, enhancer, or relationship
exists in a plant.

`anthocyanin_accessibility_demo.narrowPeak` is likewise synthetic. Its four
peaks exist to exercise imported experimental-track semantics. The example's
`EVIDENCE`, `ASSAY`, `SAMPLE`, `CONDITION`, `REPLICATE`, and `CONTROL` values
are deliberately marked as synthetic; they do not turn the fixtures into
ATAC-seq, DAP-seq, or ChIP-seq measurements.
`anthocyanin_myb_binding_demo.narrowPeak` and its `rep2` companion are also
synthetic and exist only to exercise multi-track and nested evidence
retention. Their `BINDING` declarations are not experimental claims.

The `MA0054.1 myb.Ph3` frequency matrix is the unmodified JASPAR CORE profile
for *Petunia x hybrida* MYB.Ph3. Applying a related profile to a real
anthocyanin study still requires a documented TF-family rationale, compatible
species and assembly data, and ideally independent accessibility or binding
evidence.

Interpret the outputs conservatively:

- `NEAR` reports genomic proximity, not regulation. Accessibility, binding,
  expression, chromatin contact, or other independent evidence is needed to
  strengthen a candidate link.
- `COUNT` is descriptive overlap aggregation, not motif enrichment. Region
  length, nucleotide composition, accessibility, and the chosen statistical
  background affect expected counts.
- `DEFINE MODULE` reports pairs satisfying the stated spacing, order, and
  orientation grammar. A matching pair is not proof of cooperative binding;
  the constraints need a documented biological or benchmark rationale.
- A narrowPeak interval represents a called enrichment peak from an upstream
  experiment. Its presence supports accessibility or binding only when the
  assay, sample, controls, assembly, and peak-calling provenance are suitable.
  Cis-QL preserves the supplied statistics but does not reinterpret them as
  proof of target-gene regulation.
- Overlap between two replicate peak sets is coordinate-level concordance. It
  is not IDR, does not model replicate quality or controls, and is not by
  itself a formal reproducibility assessment.
- `CONSENSUS` counts distinct input sets, not individual overlapping peaks,
  and retains the declared anchor geometry. Its optional reciprocal-overlap
  threshold rejects marginal intersections, while its summit-distance
  threshold can require closer narrowPeak support. Neither criterion is IDR.
  Choosing an anchor, support threshold, overlap fraction, and summit distance
  remains an analysis decision that must be justified.

The removed historical examples mixed eukaryotic TF models and promoter
assumptions with an *E. coli* fixture. Their language constructs remain covered
by automated tests, but they are intentionally not presented as scientific
workflows.

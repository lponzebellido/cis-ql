# Cis-QL example paths

The examples are intended to be read as small analysis paths, not as a list of
unrelated syntax fragments. Examples 07, 17, 18, 19, and 14 form the current
cis-regulatory path:

1. `07_pwm_scanning.cql` introduces a plant MYB position-frequency matrix and
   shows the difference between a relative score cutoff and statistical
   calibration. Its 90% cutoff retains six strand-specific scores at three
   loci.
2. `17_explicit_promoters.cql` constructs three strand-aware promoter windows
   from explicit TSS-relative boundaries.
3. `18_statistical_pwm_scan.cql` estimates the nucleotide background and uses
   `QVALUE <= 0.01`. It retains three calibrated strand-specific sites.
4. `19_regulatory_overlap.cql` asks which of those sites overlap a candidate
   promoter. Two sites are retained with their original coordinates and motif
   evidence.
5. `14_integrated_query.cql` is the compact end-to-end demonstration. It
   separates two promoter-supported MYB sites from one site supported by the
   annotated candidate enhancer and exports both evidence tracks.

The shared `anthocyanin_regulatory_demo` FASTA and GFF3 files are deliberately
small synthetic fixtures. They contain three annotated genes, one candidate
enhancer, two promoter-local motif instances, and one enhancer-local instance.
They demonstrate query semantics and reproducibility; they are not biological
evidence for a regulatory relationship.

The `MA0054.1 myb.Ph3` frequency matrix is an official plant MYB profile from
[JASPAR CORE](https://jaspar.elixir.no/matrix/MA0054.1/) (*Petunia x hybrida*,
SELEX). The repository stores the unmodified frequency counts shown by JASPAR.
Using a related plant profile in a real anthocyanin study would still require a
documented TF-family rationale, compatible species/assembly data, and
independent accessibility or binding evidence.

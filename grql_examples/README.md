# GRQL examples

Each file demonstrates one language operation with a small query. Run them
from the repository root:

```bash
./grql grql_examples/01_promoters.grql
```

| File | Operation | Input |
| :--- | :--- | :--- |
| `01_promoters.grql` | Select `thrA` and derive its promoter | *E. coli* FASTA and GFF3 |
| `02_pwm.grql` | Scan a PWM with a relative score threshold | synthetic FASTA, JASPAR matrix |
| `03_pvalues.grql` | Estimate a background and filter by q-value | synthetic FASTA, JASPAR matrix |
| `04_overlap.grql` | Keep motif matches that overlap promoters | *E. coli* FASTA and GFF3 |
| `05_count.grql` | Count motif matches per promoter, including zeroes | *E. coli* FASTA and GFF3 |
| `06_near.grql` | Link nearby motif matches to genes by distance | *E. coli* FASTA and GFF3 |
| `07_modules.grql` | Pair two patterns under spacing rules | *E. coli* FASTA |
| `08_filters.grql` | Filter GFF3 records by type, strand, phase, and attribute | *E. coli* GFF3 |
| `09_tracks.grql` | Load narrowPeak data and filter its fields | synthetic narrowPeak track |
| `10_track_overlap.grql` | Overlap accessibility and binding tracks | synthetic narrowPeak tracks |
| `11_consensus.grql` | Require coordinate support across two replicates | synthetic narrowPeak tracks |
| `12_orfs.grql` | Find start-to-stop patterns in complete codons | *E. coli* FASTA |
| `13_gff3.grql` | Validate and traverse GFF3; group repeated IDs | *E. coli* FASTA and GFF3 |
| `14_transcripts.grql` | Select transcript types and derive their promoters | synthetic FASTA and GFF3 |
| `15_strands.grql` | Compare reference and strand-oriented sequence | synthetic FASTA |
| `16_translate.grql` | Translate intervals with two genetic codes | synthetic FASTA |
| `17_captures.grql` | Query numbered regular-expression groups | synthetic FASTA |

`ecoli2.fna` and `genomic.gff` use accession `U00096.3`, so sequence and
annotation coordinates refer to the same assembly. The track and regulatory
FASTA fixtures are synthetic and exist only to keep expected results small.
Examples 02 and 03 use the unmodified JASPAR CORE `MA0054.1` matrix; it is an
input for testing PWM syntax, not a recommended model for *E. coli*.

Example 14 declares both `mRNA` and `transcript` as accepted input types. Its
`canonical` tag is fixture metadata queried with `SELECT ATTRIBUTE`, not a
portable GFF3 convention or an inference made by GRQL.

Example 15 finds the same bounded pattern on both strands. `SEQUENCE` retains
the bases as stored in the FASTA; `ORIENTED_SEQUENCE` reverse-complements the
negative-strand match before applying the text predicates.

Example 16 translates both strand-oriented matches with NCBI tables 1 and 4.
The TGA codon is a stop in table 1 and tryptophan in table 4, so the example
uses a protein slice to show that genetic-code choice is part of the program
rather than an implicit default. It is a translation example, not an ORF or
gene prediction.

Example 17 retains the explicit groups of a regular expression and filters
their values with `CAPTURE`, including a suffix selected with `SLICE`. Capture
offsets are relative to the matched text in its 5'-to-3' orientation, including
matches on the negative strand.

The regular expression in `12_orfs.grql` advances three nucleotides at a time
after `ATG` and rejects an in-frame stop in each repeated codon. It therefore
reaches `TAA`, `TAG`, or `TGA` in the same frame. The `LENGTH MOD 3 = 0`
condition makes that requirement visible, while `SLICE` checks the first and
last codons without assuming a fixed match length. The result is a set of
sequence matches, not a gene prediction.

Interpret interval operations literally. `OVERLAPS` reports coordinate
intersection, `NEAR` reports distance, `COUNT` reports contained matches, and
`CONSENSUS` reports support under the declared overlap and summit thresholds.
None of these operations alone establishes regulation, enrichment, binding,
or reproducibility.

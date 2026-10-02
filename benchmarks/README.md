# Benchmarks and reference checks

Cis-QL has two benchmark levels and one independent correctness suite. None is
a comparison with another language or tool.

## Complete queries

```sh
make check-examples
make benchmark
```

`check-examples` runs every file listed in `examples.json` once. It fails when
an example is missing, its result counts change, or an expected export is not
created. CI runs this target.

`benchmark` performs one warm-up and three measured executions per example.
The JSON report contains wall-clock samples, median, range, input hashes,
binary and manifest hashes, result counts, result hashes, export hashes, system
metadata, and the Git revision when available. Repeated executions must produce
identical JSON and exports.

Run a subset or save the report with:

```sh
python3 benchmarks/run_examples.py \
  --examples 01 07 13 \
  --warmups 1 \
  --repetitions 5 \
  --output benchmark.json
```

Update `examples.json` only when a deliberate language, fixture, or example
change modifies the accepted results.

## Core algorithms

```sh
make benchmark-core
```

This target measures deterministic synthetic workloads for exact motif search,
CpG scanning, PWM scanning, Smith-Waterman alignment, interval intersection,
overlap selection and counting, bounded proximity, promoter construction,
two-pattern modules, and track consensus. Each workload performs one warm-up
and reports the median of five executions as CSV.

## Correctness

```sh
make validate
```

The Python reference suite independently checks motif coordinates, PSSM
thresholds and tail probabilities, interval subtraction, overlap, nearest
selection, counting, consensus, motif modules, reading frame, and normalized
Smith-Waterman scoring. BEDTools, Biopython, and FIMO comparisons run when those
programs are installed.

Before publishing a comparison, use identical records, coordinate conventions,
strand rules, thresholds, and output definitions. Verify coordinates before
timing, pin tool and dataset versions, measure peak memory, and report repeated
runs rather than a single speed ratio.

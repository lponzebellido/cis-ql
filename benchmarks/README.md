# Cis-QL engineering benchmarks

These benchmarks are intended to measure the current implementation
reproducibly. They do not claim superiority over another language or tool.

## Commands

```sh
make benchmark-core
make benchmark
```

`benchmark-core` measures algorithm scaling with deterministic synthetic data.
It performs one warm-up followed by five timed executions and reports the
median:

- exact KMP motif search;
- CpG-island scanning;
- PWM/PSSM scanning;
- Smith-Waterman alignment;
- interval intersection.

`benchmark` executes representative complete `.cql` programs three times and
reports machine metadata, median, minimum, maximum, individual wall-clock
measurements, and result counts as JSON. It aborts if repeated executions
produce different result counts. All generated result files are kept in
temporary directories.

For publication-quality comparisons, reference tools must receive identical
input records, coordinate conventions, strand rules, thresholds, and expected
outputs. Correctness must be checked before timing.

## Correctness before comparison

Run the independent reference suite before recording benchmark results:

```sh
make validate
```

The suite uses separate Python implementations for exact overlapping motif
coordinates, log-odds PSSM thresholds, geometric interval subtraction, and the
normalized Smith-Waterman recurrence. If `bedtools` or Biopython is installed,
the compatible external checks are also executed; otherwise they are reported
as skipped.

FIMO and Cis-QL do not currently expose the same threshold statistic: FIMO
reports p-values and q-values, whereas Cis-QL uses a percentage of the
matrix-specific attainable log-odds range. A fair FIMO comparison therefore
requires a predeclared threshold-mapping protocol and is not yet part of the
automated suite.

## Publication checklist

Before describing a result as a comparison with another tool:

1. Pin tool and dataset versions.
2. Record the exact command lines and environment.
3. Convert all outputs to the same coordinate convention.
4. Compare coordinates and strand assignments before timing.
5. Run warm-ups and repeated measurements on the same machine.
6. Report runtime and peak memory with uncertainty, not only a speed ratio.

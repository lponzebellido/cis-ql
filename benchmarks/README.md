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

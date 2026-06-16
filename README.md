# Cis-QL: Cis-Regulatory Query Language

![Status: In Development](https://img.shields.io/badge/status-in%20development-orange)

**Cis-QL** is a specialized, SQL-like Domain-Specific Language (DSL) designed for bioinformatics and genomic analysis. It enables researchers to intuitively search, filter, extract, and manipulate genomic regions, regulatory motifs, and sequence features using a human-readable, declarative syntax.

---

## What is Cis-QL?

Finding patterns (like promoters, enhancers, or open reading frames) in massive genomic datasets usually requires stringing together multiple command-line tools (e.g., `bedtools`, `grep`, custom Python scripts) or writing complex code.

**Cis-QL solves this by providing a unified, declarative interface.** Instead of writing imperative code to parse files, perform string matching, calculate distances, and run set intersections, you simply *declare what you want to find*.

Cis-QL operates in two primary modes:
1. **Annotation-driven:** Works with established `GFF3` annotation files to query known genes, exons, and regulatory elements.
2. **De novo Discovery:** Works directly on raw `.fasta` files. It can discover regulatory motifs or genes on the fly using literal sequences or **Regular Expressions**, and construct a virtual annotation layer in-memory.

## What Problems Does It Solve?

- **Metagenomic Bioprospecting:** When dealing with newly sequenced environmental samples lacking curated GFF3 files, Cis-QL's *de novo* mode allows researchers to hunt for hypothetical operons or gene clusters on the fly.
- **Spatial Queries:** Easily finding motifs relative to other genomic features (e.g., "Find this motif within 200 base pairs upstream of a gene").
- **Pipeline Simplification:** Reduces complex bioinformatic pipelines into readable, reproducible, and easily shareable query scripts.

---

## How to Use It

### Compilation

Cis-QL is written in modern C++ (C++11). To compile the interpreter:

```bash
make clean
make
```

### Execution

Run a Cis-QL script (`.cql`) by passing it to the interpreter:

```bash
./cisql my_script.cql
```

Use the `--debug` flag to see the underlying AST (Abstract Syntax Tree), Symbol Table, and Intermediate Representation (IR):

```bash
./cisql my_script.cql --debug
```

---

## Syntax & Grammar Overview

Cis-QL reads like SQL but is tailored for genomics.

### 1. Loading Data
You can load sequence data (FASTA) and/or annotation data (GFF3).
```sql
LOAD SEQUENCE "genome.fasta" AS seq_1;
LOAD ANNOTATION "annotations.gff3" AS annot_1;
```

### 2. Finding Motifs (Spatial & Regex)
Find motifs anywhere in the sequence, or spatially relative to other features.
```sql
-- Literal search on the positive strand
FIND MOTIF "TATAAT" STRAND POSITIVE AS tatabox;

-- Regular expression search
FIND MOTIF "ATG(...)*?(TAA|TAG|TGA)" AS putative_genes;

-- Spatial search: Find motif near a known feature or a dynamically discovered feature
FIND MOTIF "ATGCGA" WITHIN 500 BP UPSTREAM FROM tatabox;
```

### 3. Filtering and Extracting
Extract specific biological entities and apply filters using the `WHERE` clause.
Supported properties: `LENGTH`, `SIMILARITY` (computed via Smith-Waterman local alignment).
```sql
EXTRACT GENE WHERE LENGTH > 1.5 KB;
```

### 4. Set Operations
Perform logical intersections, unions, or subtractions between feature sets.
```sql
INTERSECT EXON AND INTRON;
EXCEPT GENE FROM REGION;
```

---

## Examples of Use

### Example 1: Traditional Annotation Query
Find all genes larger than 500 base pairs in an annotated genome.
```sql
LOAD SEQUENCE "ecoli.fasta" AS genome;
LOAD ANNOTATION "ecoli.gff3" AS annotations;

EXTRACT GENE WHERE LENGTH > 500 BP;
```

### Example 2: De Novo Discovery (Regex & Spatial Targeting)
Imagine you have a raw environmental sample (no GFF3 available). You want to find hypothetical promoter sequences and then look for open reading frames (ORFs) immediately downstream.

```sql
/* 1. Load raw sequence */
LOAD SEQUENCE "metagenoma_crudo.fasta" AS sample_1;

/* 2. Discover putative promoters using regex */
FIND MOTIF "TATA[AT]A[AT]" STRAND POSITIVE AS putative_promoters;

/* 3. Look for ORFs specifically downstream of the discovered promoters */
FIND MOTIF "ATG(...)*?(TAA|TAG|TGA)" 
WITHIN 200 BP DOWNSTREAM FROM putative_promoters 
AS putative_genes 
WHERE LENGTH > 900 BP;

/* 4. Extract the final findings */
EXTRACT putative_genes;
```

---

## Current Status

**Cis-QL is actively under development.** 

Current working features:
- Lexical and Semantic analysis engines.
- Virtual `GFF3` annotation building (in-memory).
- KMP string matching and `std::regex` engine integrations.
- Local alignment (Smith-Waterman) for sequence similarity filtering.
- Fully operational set logic (Intersect, Union, Except).

*Note: This language is a research project and compiler design implementation. Expect syntax expansions and further optimizations in future releases.*

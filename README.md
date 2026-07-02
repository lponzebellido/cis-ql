# Cis-QL: Cis-Regulatory Query Language

![Status: In Development](https://img.shields.io/badge/status-in%20development-orange)

**Cis-QL** is a specialized, SQL-like Domain-Specific Language (DSL) designed for bioinformatics and genomic analysis. It enables researchers to intuitively search, filter, extract, and manipulate genomic regions, regulatory motifs, and sequence features using a human-readable, declarative syntax.

---

## What is Cis-QL?

Finding cis-regulatory elements (promoters, enhancers, transcription factor binding sites) in massive genomic datasets usually requires stringing together multiple command-line tools (e.g., `bedtools`, `grep`, custom Python scripts) or writing complex code.

**Cis-QL solves this by providing a unified, declarative interface.** Instead of writing imperative code to parse files, perform string matching, calculate distances, and run set intersections, you simply *declare what you want to find*.

Cis-QL operates in three primary modes:
1. **Annotation-driven:** Works with established `GFF3` annotation files to query known genes, exons, and regulatory elements.
2. **De novo Discovery:** Works directly on raw `.fasta` files. It can discover regulatory motifs or genes on the fly using literal sequences or **Regular Expressions**, and construct a virtual annotation layer in-memory.
3. **PWM Scanning:** Uses **Position Weight Matrices (PWMs)** from standard databases (JASPAR format) to identify transcription factor binding sites with probabilistic scoring — the industry-standard method for cis-regulatory element discovery.

## What Problems Does It Solve?

- **Transcription Factor Binding Sites:** Use pre-built PWMs (JASPAR, TRANSFAC) to find where specific TFs bind, going far beyond exact-match searches.
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

Use the `--debug` flag to see the underlying lexical analysis, AST (Abstract Syntax Tree), Symbol Table, Intermediate Representation (IR), and all intermediate results:

```bash
./cisql my_script.cql --debug
```

---

## Language Reference

### 1. Loading Data

You can load sequence data (FASTA), annotation data (GFF3), and PWM matrices (JASPAR format).

```sql
LOAD SEQUENCE "genome.fasta" AS seq_1;
LOAD ANNOTATION "annotations.gff3" AS annot_1;
LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tata_pwm;
```

### 2. Finding Motifs (Literal & Regex)

Find motifs anywhere in the sequence, or spatially relative to other features.

```sql
-- Literal search on the positive strand
FIND MOTIF "TATAAT" STRAND POSITIVE AS tatabox;

-- Regular expression search
FIND MOTIF "ATG(...)*?(TAA|TAG|TGA)" AS putative_genes;

-- Spatial search: Find motif near a known feature or a dynamically discovered feature
FIND MOTIF "ATGCGA" WITHIN 500 BP UPSTREAM FROM tatabox;
```

### 3. PWM Scanning

Scan a loaded sequence with a Position Weight Matrix to find transcription factor binding sites above a score threshold. The threshold is expressed as a percentage of the theoretical maximum PWM score (0–100%).

```sql
-- Scan both strands, default threshold (75%)
SCAN tata_pwm AS tata_sites;

-- Scan positive strand only with explicit threshold
SCAN tata_pwm STRAND POSITIVE THRESHOLD 80% AS tata_sites;

-- Scan with a more permissive threshold for exploratory analysis
SCAN sp1_pwm THRESHOLD 60% AS gc_boxes;
```

**Options (in any order, all optional):**
| Option | Description | Default |
|--------|-------------|---------|
| `STRAND POSITIVE` / `STRAND NEGATIVE` | Restrict search to one strand | Both strands |
| `THRESHOLD <value>%` | Minimum score as % of max PWM score | `75%` |
| `AS <alias>` | Name the result set for later use | — |

### 4. Filtering and Extracting

Extract specific biological entities and apply filters using the `WHERE` clause.
Supported properties: `LENGTH`, `SIMILARITY` (computed via Smith-Waterman local alignment).

```sql
EXTRACT GENE WHERE LENGTH > 1.5 KB;
EXTRACT tata_sites WHERE LENGTH >= 8 BP;
```

### 5. Set Operations

Perform logical intersections, unions, or subtractions between feature sets.

```sql
INTERSECT tata_sites AND gc_boxes;
EXCEPT GENE FROM REGION;
UNION promoters AND enhancers;
```

---

## Included PWM Matrices (JASPAR)

The `matrices/` directory includes ready-to-use Position Weight Matrices for common transcription factors:

| File | Factor | JASPAR ID | Biological Role |
|------|--------|-----------|-----------------|
| `MA0108.1_TBP.pwm` | TBP | MA0108.1 | TATA-box binding; core promoter element |
| `MA0079.5_SP1.pwm` | SP1 | MA0079.5 | GC-box binding; ubiquitous activator |
| `MA0139.1_CTCF.pwm` | CTCF | MA0139.1 | Insulator / chromatin boundary element |
| `MA0003.4_TFAP2A.pwm` | TFAP2A | MA0003.4 | Enhancer-binding; developmental regulation |
| `MA0002.2_RUNX1.pwm` | RUNX1 | MA0002.2 | Core-binding factor; hematopoiesis |

You can also load any custom PWM in JASPAR format using `LOAD MATRIX`.

---

## Examples

### Example 1: Traditional Annotation Query

Find all genes larger than 500 base pairs in an annotated genome.

```sql
LOAD SEQUENCE "ecoli.fasta" AS genome;
LOAD ANNOTATION "ecoli.gff3" AS annotations;

EXTRACT GENE WHERE LENGTH > 500 BP;
```

### Example 2: De Novo Discovery (Regex & Spatial Targeting)

Find hypothetical promoter sequences in a raw metagenome, then look for open reading frames downstream.

```sql
LOAD SEQUENCE "metagenoma_crudo.fasta" AS sample_1;

-- Discover putative promoters using regex
FIND MOTIF "TATA[AT]A[AT]" STRAND POSITIVE AS putative_promoters;

-- Look for ORFs specifically downstream of the discovered promoters
FIND MOTIF "ATG(...)*?(TAA|TAG|TGA)"
  WITHIN 200 BP DOWNSTREAM FROM putative_promoters
  AS putative_genes
  WHERE LENGTH > 900 BP;

EXTRACT putative_genes;
```

### Example 3: PWM-Based Cis-Regulatory Discovery

Use a JASPAR matrix to find TATA box binding sites with probabilistic scoring.

```sql
LOAD SEQUENCE "data_examples/ecoli.fasta" AS genome;

LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tata_pwm;

-- Find high-confidence TATA box sites (>= 80% of max PWM score)
SCAN tata_pwm STRAND POSITIVE THRESHOLD 80% AS tata_sites;

EXTRACT tata_sites;
```

### Example 4: Combining PWM Scanning with Set Operations

Find putative promoter regions where a TATA box overlaps with a GC-box.

```sql
LOAD SEQUENCE "genome.fasta" AS genome;

LOAD MATRIX "matrices/MA0108.1_TBP.pwm" AS tata_pwm;
LOAD MATRIX "matrices/MA0079.5_SP1.pwm" AS sp1_pwm;

SCAN tata_pwm THRESHOLD 80% AS tata_sites;
SCAN sp1_pwm THRESHOLD 75% AS gc_boxes;

-- Promoters where both elements co-occur
INTERSECT tata_sites AND gc_boxes;
```

---

## Syntax & Grammar Overview

Cis-QL reads like SQL but is tailored for genomics.

### Formal Grammar (CFG)

| Non-Terminal | | Expansion |
| :--- | :---: | :--- |
| `Program` | → | `StatementList` |
| `StatementList` | → | `Statement` `StatementList` <br> \| `λ` (epsilon) |
| `Statement` | → | `LoadStmt` \| `FindStmt` \| `ExtractStmt` \| `SetOperationStmt` \| `ScanStmt` |
| `LoadStmt` | → | **`LOAD`** (**`SEQUENCE`** \| **`ANNOTATION`** \| **`MATRIX`**) *`STRING`* **`AS`** *`ID`* **`;`** |
| `FindStmt` | → | **`FIND MOTIF`** *`STRING`* `FindOpts` `AliasOpt` `WhereClause` **`;`** |
| `FindOpts` | → | `FindOpt` `FindOpts` <br> \| `λ` |
| `FindOpt` | → | **`WITHIN`** *`NUM`* `Unit` `Direction` **`FROM`** `EntityRef` `EntityName` <br> \| **`STRAND`** `StrandType` <br> \| **`CHR`** *`STRING`* |
| `ScanStmt` | → | **`SCAN`** *`ID`* `ScanOpts` `AliasOpt` `WhereClause` **`;`** |
| `ScanOpts` | → | `ScanOpt` `ScanOpts` <br> \| `λ` |
| `ScanOpt` | → | **`STRAND`** `StrandType` <br> \| **`THRESHOLD`** *`NUM`* **`%`** |
| `AliasOpt` | → | **`AS`** *`ID`* <br> \| `λ` |
| `SetOperationStmt`| → | `SetOp` `EntityRef` **`AND`** `EntityRef` `WhereClause` **`;`** |
| `SetOp` | → | **`INTERSECT`** \| **`UNION`** \| **`EXCEPT`** |
| `ExtractStmt` | → | **`EXTRACT`** `EntityRef` `WhereClause` **`;`** |
| `WhereClause` | → | **`WHERE`** `Condition` <br> \| `λ` |
| `Condition` | → | `Term` `ConditionPrime` |
| `ConditionPrime`| → | **`OR`** `Term` `ConditionPrime` <br> \| `λ` |
| `Term` | → | `Factor` `TermPrime` |
| `TermPrime` | → | **`AND`** `Factor` `TermPrime` <br> \| `λ` |
| `Factor` | → | **`NOT`** `Factor` <br> \| `SimpleCondition` <br> \| **`(`** `Condition` **`)`** |
| `SimpleCondition`| → | `Property` `RelOp` `Value` |
| `Property` | → | **`LENGTH`** \| **`SIMILARITY`** |
| `RelOp` | → | **`>`** \| **`<`** \| **`>=`** \| **`<=`** \| **`=`** |
| `Value` | → | *`NUM`* `Unit` <br> \| *`FLOAT`* **`%`** <br> \| *`NUM`* <br> \| *`STRING`* |
| `Unit` | → | **`BP`** \| **`KB`** \| **`MB`** <br> \| `λ` |
| `Direction` | → | **`UPSTREAM`** \| **`DOWNSTREAM`** |
| `Entity` | → | **`GENE`** \| **`PROMOTER`** \| **`ENHANCER`** \| **`EXON`** \| **`INTRON`** \| **`UTR`** \| **`TSS`** \| **`CDS`** \| **`REGION`** |
| `EntityRef` | → | `Entity` \| *`ID`* |
| `EntityName` | → | *`STRING`* <br> \| `λ` |
| `StrandType` | → | **`POSITIVE`** \| **`NEGATIVE`** |

*(Keywords and terminals are highlighted in **`bold code`**. Identifiers and literals are in *`italic code`*. Non-terminals are in `regular code`)*

---

## Current Status

**Cis-QL is actively under development.**

Current working features:
- Lexical, syntactic, and semantic analysis engines.
- Virtual `GFF3` annotation building (in-memory).
- KMP string matching and `std::regex` engine integrations.
- **PWM/PSSM scanning** (JASPAR format) with log-odds scoring and sliding-window search on both strands.
- Local alignment (Smith-Waterman) for sequence similarity filtering.
- Fully operational set logic (Intersect, Union, Except).
- `--debug` mode exposing the full compilation pipeline (tokens, AST, IR, intermediate results).

*Note: This language is a research project and compiler design implementation. Expect syntax expansions and further optimizations in future releases.*

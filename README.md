# mutme

`mutme` is a lightweight Python command-line tool for detecting and annotating
mutations of interest in viral consensus genome sequences.

It runs `Nextclade` against a user-supplied reference genome and GFF3 annotation,
then matches amino-acid substitutions, deletions, insertions, and premature stop
codons against a curated mutation annotation table.

`mutme` helps you answer a simple question at scale:

> “Which mutations of interest are present in my consensus genomes?" 

It can also provide generic annotations associated with each mutation, for example antiviral-resistance mutations, lineage markers, or any other genotype/phenotype you're interested in.

You provide:
- a reference genome (FASTA)
- a genome annotation (GFF3)
- a mutation annotation table (CSV/TSV)
- a query sequence file containing the consensus genomes (FASTA)

> [!IMPORTANT]
`mutme` is intentionally generic and database-agnostic. It does not attempt to
interpret biological meaning beyond what you encode in your annotation table.

## Key Features

- Exact string matching against mutation-of-interest tables
- Supports:
  - Amino-acid substitutions
  - Deletions
  - Insertions
  - Premature stop codons
- Robust handling of multiple annotation rows per mutation
- Long-format output suitable for downstream analysis

## Table of Contents

- [Installation](#installation)
  - [Requirements](#requirements)
  - [Install from source](#install-from-source)
  - [Install with dependencies](#install-with-dependencies)
- [Quick Start](#quick-start)
- [Annotations Table](#annotations-table)
  - [How to write the annotations table](#how-to-write-the-annotations-table)
  - [Mutation encoding (what goes in the mutation column)](#mutation-encoding-what-goes-in-the-mutation-column)
- [Output Table](#output-table)
- [Examples](#examples)
  - [I only care if these mutations are present](#i-only-care-if-these-mutations-are-present)
  - [I want to annotate mutations with mAb susceptibility values](#i-want-to-annotate-mutations-with-mab-susceptibility-values)
- [Command reference](#command-reference)
  - [mutme run](#mutme-run)
  - [mutme subset-gff3](#mutme-subset-gff3)
- [External Dependencies](#external-dependencies)
- [Citation](#citation)

## Installation

### Requirements

- Python >= 3.10
- Nextclade >= v3.18 installation available on `$PATH`

### Install from source

```bash
pip install mutme
```

### Install with dependencies

```bash
mamba create -n mutme -c conda-forge -c bioconda python=3.14 nextclade=3.18
mamba activate mutme
pip install mutme
```

## Quick Start

```bash
mutme run \
  --sequences sequences.fasta \
  --reference reference.fasta \
  --gff reference.gff3 \
  --annotations mutations.csv \
  --output results.csv
```

Notes:

- `--sequences` can contain one or many consensus sequences.
- Output will contain a `seq_name` column so you can see which result belongs to which input sequence.

## Annotations Table

### How to write the annotations table

Your annotation table must have a `mutation` column (case-insensitive). Everything else is optional — add as many extra columns as you like. One special case is a `comment` column, which can contain e.g. links to the mutation phenotype study. 

### Mutation encoding (what goes in the mutation column)

Gene names must be defined in the GFF3 (see below). If CDS are included ensure that the `Name` attribute matches the `gene` name used in the annotation table. Nextclade will prefer the `Name` attribute for amino acid mutation prefixes.

- Substitution: `{gene}:{aa}{pos}{aa}` (`S:N87Y`)
- Deletion: `{gene}:{pos}-` or `{gene}:{aa}{pos}-` (`S:87-` or `S:N87-`)
- Insertion: `{gene}:{pos}{aa-ins}` (`S:214:EPE`)
- Stop codon: `{gene}:{aa}{pos}*` or `{gene}:{pos}` (`S:N87*` or `S:87`)

## Output Table

Default output columns (always present):

- `seq_name` - the FASTA record name (so multiple sequences are kept separate)
- `seq_quality` - Nextclade QC overall status
- `mutation` - the matched mutation from your annotation table

If your annotation table has extra columns, those columns are included too.

## Examples

### “I only care if these mutations are present”

#### Input annotation table (CSV)

```csv
mutation
S:N87Y
S:87-
S:214:EPE
S:87
```

#### Output table (example)

```csv
seq_name,seq_quality,mutation
sample_01,good,S:N87Y
sample_01,good,S:87-
sample_02,mediocre,S:214:EPE
sample_03,good,S:87
```


### "I want to annotate mutations with mAb susceptibility values”

Here’s a more “real” example with multiple numeric columns. The values below are just an illustration of what your table can look like.

#### Input annotation table (CSV)

```csv
mutation,casirivimab_fold,imdevimab_fold,sotrovimab_fold,comment
S:E484K,12.5,8.1,1.2,Reduced neutralization for several mAbs
S:K417N,3.4,1.1,0.9,May reduce some class 1 mAbs
S:214:EPE,,,,Insertion seen in some lineages
S:87,,"",,Premature stop at position 87
S:87-,,"",,Deletion at position 87
```

> [!TIP]
Columns can be floats, ints, strings, or left empty. If you include a comment column, you can choose to carry it into the output with `--include-comments`.

#### Output table (example with `--include-comments`)

```csv
seq_name,seq_quality,mutation,casirivimab_fold,imdevimab_fold,sotrovimab_fold,comment
sample_01,good,S:E484K,12.5,8.1,1.2,Reduced neutralization for several mAbs
sample_01,good,S:214:EPE,,,,Insertion seen in some lineages
sample_02,good,S:K417N,3.4,1.1,0.9,May reduce some class 1 mAbs
sample_03,mediocre,S:87,,"",,Premature stop at position 87
sample_03,mediocre,S:87-,,"",,Deletion at position 87
```

> [!NOTE]
The output keeps your original columns (and their values) attached to each mutation hit. A sequence with no hits won’t appear in the output unless you choose to emit empty rows (not enabled by default).

## Command reference

`mutme` provides two main commands: `run` (the core workflow) and `subset-gff3`
(a helper for trimming GFF3 files based on mutations in an annotation table).

### `mutme run`

Run Nextclade using a custom reference and GFF3, then match detected amino-acid
mutations against an annotation table.

#### Basic usage

```bash
mutme run \
  --sequences sequences.fasta \
  --reference reference.fasta \
  --gff reference.gff3 \
  --annotations mutations.csv \
  --output results.csv
```

> [!NOTE]
Input sequence file `sequences.fasta` may contain **one or many consensus sequences**. Each FASTA record is processed independently and results are distinguished by sequence name in the output table.

---

#### Required options

Either `--reference` + `--gff` or `--nextclade-tsv` must be provided.

| Option | Description |
|------|-------------|
| `--sequences`, `-s` | Input FASTA with one or more consensus sequences |
| `--reference`, `-r` | Reference genome FASTA |
| `--gff`, `-g` | Genome annotation in GFF3 format |
| `--nextclade-tsv`, `-n` | Precomputed Nextclade TSV output (alternative to `--reference` + `--gff`) |
| `--annotations`, `-a` | Mutation annotation table (CSV/TSV) |
| `--output`, `-o` | Output file path (CSV/TSV) |

---

#### Common optional options

| Option | Description |
|------|-------------|
| `--include-comments`, `-c` | Include a `comment` column in output if present in annotation table |
| `--alignment-preset`, `-p` | Nextclade alignment preset (`default`, `high-diversity`, `short-sequences`) |
|`--nextclade-threads`, `-t` | Number of threads to use for Nextclade (default: all available) |
| `--nextclade-extra-args` | Extra arguments passed directly to Nextclade |
| `--nextclade-keep-tsv` | Keep intermediate Nextclade TSV output |
| `--nextclade-bin` | Path or name of Nextclade executable (default `nextclade`) |
| `--annotations-delimiter` | Delimiter used by annotation table (default `,`, use `\t` for TSV) |
| `--output-delimiter` | Delimiter for output table (default `,`) |

---

#### Example: running Nextclade with limited threads

```bash
mutme run \
  -s sequences.fasta \
  -r reference.fasta \
  -g reference.gff3 \
  -a mutations.csv \
  -o results.csv \
  -t 8
```

#### Example: precomputed Nextclade TSV

```bash
mutme run \
  -s sequences.fasta \
  -n results.nextclade.tsv \
  -a mutations.csv \
  -o results.csv
```

#### Example: TSV annotations, CSV output

```bash
mutme run \
  -s sequences.fasta \
  -r reference.fasta \
  -g reference.gff3 \
  -a mutations.tsv \
  -o results.csv \
  --annotations-delimiter '\t'
```

---

#### Example: keep Nextclade output for debugging

```bash
mutme run \
  -s sequences.fasta \
  -r reference.fasta \
  -g reference.gff3 \
  -a mutations.csv \
  -o results.csv \
  --nextclade-keep-tsv
```

---

### `mutme subset-gff3`

Subset a GFF3 file to only genes (and optionally CDS features) referenced by
mutation prefixes in an annotation table.

#### Basic usage

```bash
mutme subset-gff3 reference.gff3 mutations.csv --output subset.gff3
```

---

#### How it works

- Reads the annotation table
- Extracts prefixes from the mutation column (e.g. `S` from `S:E484K`)
- Keeps matching GFF3 records using `gene` or `CDS` features
- Preserves GFF3 headers and comments

---

#### Required arguments

| Argument | Description |
|---------|-------------|
| `gff3_path` | Input GFF3 file |
| `table_path` | Annotation table (CSV/TSV) |

---

#### Common optional options

| Option | Description |
|------|-------------|
| `--output`, `-o` | Output GFF3 path (if omitted, no file is written) |
| `--delimiter`, `-d` | Annotation able delimiter (default `,`, use `\t` for TSV) |
| `--mutation-column` | Name of mutation column (default `mutation`) |
| `--prefix-separator` | Separator used in mutation strings (default `:`) |
| `--no-strip-prefixes` | Do not strip whitespace around prefixes |
| `--gene-feature-type` | GFF3 feature type treated as genes (default `gene`) |
| `--gene-attr-key` | GFF3 attribute used to match gene features to annotation prefix (default `gene`) |
| `--include-cds` | Also include CDS features |
| `--cds-feature-type` | GFF3 feature type for CDS records (default `CDS`) |
| `--cds-gene-attr-key` | GFF3 attribute key used to match CDS features (default `gene`)|
| `--json` | Emit machine-readable JSON summary instead of text |

---

#### Example: subset genes only

```bash
mutme subset-gff3 \
  reference.gff3 \
  mutations.csv \
  --output subset.gff3
```

---

#### Example: include CDS records and emit JSON summary

```bash
mutme subset-gff3 \
  reference.gff3 \
  mutations.tsv \
  --include-cds \
  --output subset.gff3 \
  --json
```

---

#### Example: inspect prefixes without writing a file

```bash
mutme subset-gff3 \
  reference.gff3 \
  mutations.csv
```

This prints summary statistics but does not write an output GFF3.

## External Dependencies

`mutme` relies on `Nextclade` for sequence alignment and mutation calling.

> [Aksamentov et al. (2021)](https://doi.org/10.21105/joss.03773) - Nextclade: clade assignment, mutation calling and quality control for viral genomes - Journal of Open Source Software

## Citation

If you use `mutme` in published work, please cite `Nextclade` (see link above) and acknowledge this repository.
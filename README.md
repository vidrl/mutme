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
> `mutme` is intentionally generic and database-agnostic. It does not attempt to
interpret biological meaning beyond what you encode in your annotation table.

## Key features

- String matching against mutation-of-interest tables
- Supports amino-acid mutations:
  - Substitutions
  - Deletions
  - Insertions
  - Premature stop codons
- Robust handling of multiple annotation rows per mutation
- Long-format output suitable for downstream analysis

## Table of contents

- [Installation](#installation)
  - [Requirements](#requirements)
  - [Install from package index](#install-from-package-index)
  - [Install with dependencies](#install-with-dependencies)
- [Quick start](#quick-start)
- [Annotations table](#annotations-table)
  - [How to write the annotations table](#how-to-write-the-annotations-table)
  - [Mutation encoding (what goes in the mutation column)](#mutation-encoding-what-goes-in-the-mutation-column)
  - [Wildcard matching X (optional)](#wildcard-matching-x-optional)
  - [Deletion ranges (optional)](#deletion-ranges-optional)
  - [Linked mutations (optional)](#linked-mutations-optional)
- [Output table](#output-table)
- [Examples](#examples)
  - [I only care if these mutations are present](#i-only-care-if-these-mutations-are-present)
  - [I want to annotate mutations with mAb susceptibility values](#i-want-to-annotate-mutations-with-mab-susceptibility-values)
- [Command reference](#command-reference)
  - [mutme run](#mutme-run)
  - [mutme subset-gff3](#mutme-subset-gff3)
- [External dependencies](#external-dependencies)
- [Citation](#citation)

## Installation

### Requirements

- Python >= 3.10
- Nextclade >= v3.18 installation available on `$PATH`

### Install from package index

```bash
pip install mutme
```

### Install with dependencies

Latest version available on Anaconda channel with dependencies (`conda/mamba`).

```bash
mamba install -c conda-forge -c bioconda -c esteinig mutme
```

## Quick start

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

## Annotations table

### How to write the annotations table

Your annotation table must have a `mutation` column (case-insensitive). Everything else is optional — add as many extra columns as you like. Two special cases exist:

- a `comment` column, which can contain e.g. links to the mutation phenotype study or database information
- a `linked_mutations` column, which can contain required co-occurring mutations for a mutation-of-interest

### Mutation encoding (what goes in the mutation column)

Gene names must be defined in the GFF3 (see below). If CDS are included, ensure that the `Name` attribute matches the `gene` name used in the annotation table. Nextclade will prefer the `Name` attribute for amino acid mutation prefixes.

| Mutation type | Encoding | Example |
|--------------|----------|---------|
| Substitution | `{gene}:{aa}{pos}{aa}` | `S:N87Y` |
| Deletion | `{gene}:{pos}-` or `{gene}:{aa}{pos}-` or `{gene}:del{start}-{stop}` | `S:87-` or `S:N87-` or `S:del87-99` |
| Insertion | `{gene}:{pos}{aa-ins}` | `S:214:EPE` |
| Stop codon | `{gene}:{aa}{pos}*` or `{gene}:{pos}` | `S:N87*` or `S:87` |


### Wildcard matching `X` (optional)

If enabled, `mutme` can treat `X` in your **annotation table** as a wildcard meaning “any single amino acid” when matching against detected Nextclade mutations.

Supported forms (annotation table only):

- Substitution: `{gene}:{ref}{pos}X`, example: `S:N87X` matches `S:N87Y`, `S:N87F`, etc.
- Insertion: `{gene}:{pos}:{ins}` where `{ins}` may contain one or more `X` characters, example: `S:345:NXY` matches `S:345:NQY`, `S:345:NRY`, etc.

Notes:

- Wildcards apply **only** to substitutions and insertions.
- Deletions and stop codons are not wildcard-matched.
- Wildcards are **disabled by default**; enable with `--allow-x-wildcards`.
- The set of amino-acid characters that `X` can match is configurable via `--x-charset`.

> [!WARNING]
> Wildcard matching applies only to `X` in annotation-table mutations (e.g. `S:N87X`, `S:345:NXY`). `X` in reference amino acids (such as `S:X87N`) will not act as a wildcard.

### Deletion ranges (optional)

In addition to single-position deletions, `mutme` supports a convenience syntax for specifying **contiguous deletion ranges** in the annotation table.

#### Syntax

`{gene}:del{start}-{stop}`

- `{start}` and `{stop}` are inclusive amino-acid indices (1-based).
- The range is expanded internally into individual deletion positions.
- Full deletion ranges are **enabled by default**; disable with `--no-require-full-del-ranges`
- Output rows correspond to individual deletion positions, not a single collapsed range.

#### Annotation table

```csv
mutation,phenotype
S:del87-89,Important deletion region
```

This is internally expanded to:

- `S:87-`
- `S:88-`
- `S:89-`

You do **not** need to write these individually.

#### Full-range requirement

By default, deletion ranges are treated as **linked**: All positions in the range must be deleted for the annotation to match.

For example:

- If a sequence contains only `S:N87-` = **no match**
- If a sequence contains `S:N87-` and `S:N88-` and `S:N89-` all three deletion positions are reported as matches

This behavior can be disabled with `--no-require-full-del-ranges`. When disabled, range rows are still expanded, but each deletion position can match independently. 

### Linked mutations (optional)

The annotation table may optionally contain a column named:

```csv
linked_mutations
```

This column allows you to require that **additional mutations must be present** in the same sequence for a given annotation row to match.

#### Syntax

- Mutations are semicolon-separated (`;`).
- All listed mutations must be present.
- Matching is based on exact mutation strings (as detected by Nextclade and normalized by `mutme`).

#### Example:

```csv
mutation,phenotype,linked_mutations
S:E484K,Escape mutation,S:K417N;S:N501Y
```

This means: `S:E484K` will only be reported if BOTH `S:K417N` AND `S:N501Y` are also present.


#### Important: Multiple rows for the same mutation

If the same mutation appears multiple times in the annotation table, `mutme` merges (unions) all linked mutations across rows.

#### Example

```csv
mutation,phenotype,linked_mutations
S:E484K,Escape,S:K417N
S:E484K,Escape,S:N501Y
```

Internally, the logic becomes: `S:E484K` requires: `S:K417N` **AND** `S:N501Y`

>[!WARNING]:
> Multiple rows for the same mutation do NOT create alternative requirements. Instead, their `linked_mutations` are merged (unioned) into a single combined requirement.

- If only `S:K417N` is present = **no match**
- If only `S:N501Y` is present = **no match**

This means for now:

- You cannot express “either A or B” by adding two rows of the same target mutation with different linked mutations
- All linked mutations across duplicate rows must co-occur for a match to be output.

If you need alternative logic you must encode that outside of `mutme` or restructure your annotation table accordingly. Currently only AND semantics are supported. OR or XOR logic is not implemented. Additional semantics for linked mutation logic will be introduced in the next version.

#### Interaction with deletion ranges

Deletion ranges automatically generate a linked requirement when full-rangemode is enabled (default).

For example:

```csv
mutation,phenotype
S:del87-89,Important deletion region
```

Internally expands to:

- `S:87-`
- `S:88-`
- `S:89-`

and behaves as if `linked_mutations` was `S:87-;S:88-;S:89-`. Thus, partial deletions within the range will not trigger the annotation unless all positions in the range are deleted (unless `--no-require-full-del-ranges` is used). 

Important: 

Only ranges matching `{gene}:del{start}-{stop}` (with permissive whitespace, e.g. `S:del87 - 89`) are treated as ranges. Other `del` strings (`S:delA-B`) are treated as literal mutation keys and will not appear in outpouts as they are not conforming to the standard mutation annotation formnat.

## Output table

Default output columns (always present):

- `seq_name` - the FASTA record name (so multiple sequences are kept separate)
- `seq_quality` - Nextclade overall QC status
- `mutation` - the matched mutation from your annotation table

If your annotation table has extra genotype/phenotype columns, those columns are included. Mutation matching is independent of Nextclade QC status; mutations are reported even for sequences with poor QC.

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
> Columns can be floats, ints, strings, or left empty. If you include a comment column, you can choose to carry it into the output with `--include-comments`.

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
> By default the output keeps all annotations columns (disable with `--no-output-all-columns`), as well as your input sequences (disable with: `--no-output-all-sequences`). 
Sequences without detected mutations will have empty annotation fields. If any output value includes an output delimiter (`--output-delimiter`) the value will be quoted (`""`).

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
> Input sequence file `sequences.fasta` may contain **one or many consensus sequences**. Each FASTA record is processed independently and results are distinguished by sequence name in the output table.

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
| `--nextclade-preset`, `-p` | Nextclade alignment preset (`default`, `high-diversity`, `short-sequences`) |
| `--nextclade-threads`, `-t` | Number of threads to use for Nextclade (default: all available) |
| `--nextclade-extra-args` | Extra arguments passed directly to Nextclade |
| `--nextclade-keep-tsv` | Keep intermediate Nextclade TSV output |
| `--nextclade-bin` | Path or name of Nextclade executable (default `nextclade`) |
| `--annotations-delimiter` | Delimiter used by annotation table (default `,`, use `\t` for TSV) |
| `--output-delimiter` | Delimiter for output table (default `,`) |
| `--allow-x-wildcards` | Treat `X` in annotation-table substitutions/insertions as a wildcard for any single amino acid |
| `--x-charset` | Allowed amino-acid characters that `X` can match when `--allow-x-wildcards` is enabled (default: 20 canonical AAs) |
| `--no-require-full-del-ranges` | Do not require all positions in a deletion range to be present for a match (default: enabled) |
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

#### Example: enable `X` wildcard matching

```bash
mutme run \
  --sequences sequences.fasta \
  --reference reference.fasta \
  --gff reference.gff3 \
  --annotations mutations.csv \
  --output results.csv \
  --allow-x-wildcards
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

## External dependencies

`mutme` relies on `Nextclade` for sequence alignment and mutation calling.

> [Aksamentov et al. (2021)](https://doi.org/10.21105/joss.03773) - Nextclade: clade assignment, mutation calling and quality control for viral genomes - Journal of Open Source Software

## Citation

If you use `mutme` in published work, please cite `Nextclade` and acknowledge this repository.
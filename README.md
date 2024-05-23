# Snakemake workflow: `QAPA_snakemake`


[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)

## Overview

A Snakemake workflow for running QAPA (**TODO CITATION**) for 3'UTR alternative polyA site quantification, plus optional additional steps to perform differential polyA site usage analysis with DEXSeq (as described in paper) and compute condition-wise summary PAS usages (mean, median, deltas between specified contrasts).

The workflow currently has the following functionality: (along with custom extensions to QAPA)

- Build QAPA annotation files and Salmon indices from scratch for given input reference annotation
  - Optionally add genome sequence as decoys to transcriptome index (see doi: 10.1186/s13059-020-02151-8 for)
- Use pre-provided annotations/indices (e.g. from a previous run/QAPA repository) to run Salmon quantification & qapa quant
- Generate isoform count matrices using tximport
- Calculate summary values for QAPA computed polyA usage % (PAU)
  - condition-wise mean and median PAU
  - differences in mean and median PAU between specified contrasts.
- Perform differential polyA site usage analysis with DEXSeq

NOTE: This workflow currently uses a forked version of QAPA available at x. The main changes are:

- Option to add genome sequence as a decoys to transcriptome index built by Salmon.
- Attempt at addressing [issue 49](https://github.com/morrislab/qapa/issues/49)
- Add penultimate exons to transcript models quantified by Salmon (inspired by LABRAT)

The above features are completely optional (with exception of bug fix attempt), so a standard QAPA workflow can be configured if desired.

## Table of Contents

- [Prerequisites and Installation](#prerequisites-and-installation)
- [Example Run](#example-run)
  - [Dry Run](#dry-run)
  - [Local Run](#local-run)
  - [Output](#output)
- [General Configuration](#configuration)
  - [Sample Table](#sample-table)
  - [Differential Usage](#differential-usage)
  - [Re-using previously generated references](#reusing-previously-generated-references) TODO: DECIDE HERE OR BELOW
- [Configuring alternative run modes](#configuring-alternative-run-modes)
  - [Use the annotation files provided by QAPA](#use-annotation-files-provided-with-qapa-repository)
  - [Generate 3'UTR annotation using custom BED file of polyA sites](#generate-3utr-annotation-using-custom-bed-file-of-polya-sites)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)

## Prerequisites and Installation

Before using this pipeline, make sure you have the following software and tools installed:

- [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- [Singularity](https://docs.sylabs.io/guides/latest/user-guide/quick_start.html#quick-installation-steps)


Clone this repository to your local machine: 
**TODO: UPDATE**

```bash
git clone https://github.com/yourusername/your-repo.git
cd your-repo
```

[Return to table of contents](#table-of-contents)

## Example Run

The following steps will work through an example run of the complete workflow of pipeline. This mimics the tutorial in the main repo's README, and extend executing the following steps:

- Construct reference of APA isoforms using Gencode transcript annotations, manual poly(A) site annotations and PolyASite 2.0 poly(A) site database (i.e. the standard workflow with `qapa build`)
- Construct FASTA file of APA isoform sequences with genome decoy sequences appended (custom extended version of `qapa fasta`)
- Quantify APA isoforms using Salmon (with additional bias correction flags enabled)
- Generate QAPA results table with PAS quantification across all samples (`qapa quant`)
- Generate isoform count matrices ready for use with differential usage tools (using `tximport`)
- Perform a simple differential APA analysis between two conditions (using `DEXSeq`)

### Dry Run

First, we perform a 'dry run' to check that the configuration YAML file is specified correctly and to report the execution plan. Run the following command using the pre-provided config file for the test data (`config/config.QAPA.yaml`):

```bash
snakemake -p --configfile config/config.QAPA.yaml --dryrun --use-singularity
```

We should see a successful execution trace (e.g. x)

[Return to table of contents](#table-of-contents)

### Local run

To run the pipeline on your local machine with Singularity containers, use the following command:

```bash

snakemake --use-singularity -p  --configfile path/to/my_config.yaml --cores n
```

- replace `n` with the number of cores for parallel processing.

This takes about x mins (on my bog standard laptop running Windows Subsystem for Linux)

[Return to table of contents](#table-of-contents)

### Output

Pipeline output is stored in subdirectories under the `**TODO**` directory, which is controlled by the `out_dir` flag in the config YAML file:

```bash

TODO: INSERT TREE DIRECTORY STRUCTURE

```

- `annotation` - metadata and annotation files required for downstream steps (computed from input annotation files)
- `qapa_build` - APA isoform reference generated by qapa build
- `qapa_fasta` - FASTA file containing APA isoform sequences (output of `qapa fasta`)
- `salmon` - Outputs generated by salmon index and quant
- `differential_apa` - Output of tximport (count matrices) and DEXSeq (differential usage results table)
  - xx

More complete/detailed description of output files to come...

[Return to table of contents](#table-of-contents)

## Configuration

Below are general instructions for configuring the pipeline to run with your own inputs. Most parameters in the config YAML file are described using comments, but key inputs will be described below:

### Sample Table

Defines input samples and their associated metadata. An example is provided at `config/samplesheet_example.csv` (**TODO: INSERT LINK**). the minimal, mandatory columns are:

- `sample_name` - unique identifier for a given sample
- `condition` - key to group samples belonging to the same experimental condition.
- `fastq1` - path to first mate of paired-end sequencing run in gzipped FASTQ format. If sample is single-end, fill this column with the path to the gzipped FASTQ file
- `fastq2` -  path to second mate of paired-end sequencing run in gzipped FASTQ format. If sample is single-end, leave this column blank/empty

#### Sample Table Notes

- Additional columns are accepted and only used by the pipeline if covariates are specified for differential usage (see differential usage section)
- It is possible to pass a mixture of paired-end and single-end samples to the pipeline

### Differential usage

Differential PAS usage with DEXSeq can be toggled by setting the `run_differential_apa` option in the config file to True/False.

#### 'Contrasts' Table

The '**contrasts**' table is used to specify the two-way comparison(s) of interest. An example is provided at `config/contrasts_example.csv`. The minimal columns are:

- `contrast_name` - unique identifier for the comparison of interest (used to differentiate output tables)
- `column_name` - name of column in sample table from which to extract group identifiers. There is no real need to change this from 'condition' in normal circumstances
- `base_key` - group identifier in 'column_name' column of sample table corresponding to the 'base'/'control' (i.e. denominator) condition in the comparison
- `contrast_key` - group identifier in 'column_name' column of sample table corresponding to the 'contrast'/'treatment' (i.e. numerator) condition in the comparison

#### Adding covariates to model

If you wish to perform differential isoform usage accounting for additional covariates, you must append these columns to the sample table. You should then subsequently update `dexseq_formula_full` and `dexseq_formula_reduced` in the config file to include the additional interaction terms. See DEXSeq documentation for full details on how to do this. A simple example below demonstrates how to construct the formulas if you wanted to add `batch` as an additional covariate in the model.

```yaml
# simple comparison of exon usage with respect to condition

dexseq_formula_full: "~ sample_name + exon + condition:exon"
dexseq_formula_reduced: "~ sample_name + exon"

# additionally account for 'batch' variable in both full and null model, so any differences in exon usage due to batch are considered in both models (and are only testing for differences in exon usage that can be accounted for by 'condition') 
# Note the 'batch' column should be appended to the sample table

dexseq_formula_full: "~ sample_name + exon + batch:exon + condition:exon"
dexseq_formula_reduced: "~ sample_name + exon + batch:exon"
```

#### Differential Usage Notes

- Fold changes, mean and median delta PAU values (WIP) are calculated with contrast/base or contrast - base.
- To specify multiple comparisons of interest, simply add additional rows to the contrasts table!

## Reusing previously generated references

To prevent unnecessary computation on successive runs using the same reference annotation, it is possible to provide pre-computed Salmon indices and reference files for a given run and skip straight to Salmon quantification. This can be achieved by 


## Configuring alternative run modes

### Use annotation files provided with QAPA repository

### Generate 3'UTR annotation using custom BED file of polyA sites



### xx

## Troubleshooting

If you encounter any issues or have questions about the pipeline, please check the Issues section of this repository. If you don't find a solution, feel free to open a new issue.

## Contributing

If you'd like to contribute to this pipeline, please follow the guidelines in CONTRIBUTING.md.

## License
This pipeline is licensed under the MIT License.
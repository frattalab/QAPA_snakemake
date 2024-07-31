# Snakemake workflow: `QAPA_snakemake`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)

## Overview

A Snakemake workflow to run [QAPA](https://doi.org/10.1186/s13059-018-1414-4) for 3'UTR alternative polyA site quantification, plus optional additional steps to perform differential polyA site usage analysis with DEXSeq (as described in QAPA paper) and (**COMING SOON**) compute condition-wise summary PAS usages (mean, median, deltas between specified contrasts).

Inclusive of custom, optional extensions to QAPA, the workflow currently has the following functionality:

- Build QAPA annotation files and Salmon indices from scratch for given input reference annotation
  - Optionally add genome sequence as decoys to transcriptome index ([previously proposed to improve quantification accuracy](https://doi.org/10.1186/s13059-020-02151-8))
- Use pre-provided annotations/indices (e.g. from a previous run/QAPA repository) to run Salmon quantification & qapa quant
- Generate isoform count matrices using tximport
- Perform differential polyA site usage analysis with DEXSeq

**NOTE**: This workflow currently uses a [forked version of QAPA](https://github.com/SamBryce-Smith/qapa). The main changes are:

- Option to add genome sequence as a decoys to transcriptome index built by Salmon.
- Attempt at addressing [issue 49](https://github.com/morrislab/qapa/issues/49)
- (**EXPERIMENTAL**) Add penultimate exons to transcript models quantified by Salmon (inspired by LABRAT)

The above features are completely optional (with exception of bug fix attempt), so a standard QAPA workflow can be configured if desired.

## Table of Contents

- [Prerequisites and Installation](#prerequisites-and-installation)
- [Example Run](#example-run)
  - [Dry Run](#dry-run)
  - [Local Run](#local-run)
  - [Output (TODO)](#output-todo)
- [General Configuration](#configuration)
  - [Sample Table](#sample-table)
  - [Differential Usage](#differential-usage)
  - [Re-using previously generated references](#reusing-previously-generated-references)
- [Configuring alternative run modes (TODO)](#configuring-alternative-run-modes-todo)
  - [Use the annotation files provided by QAPA](#use-annotation-files-provided-with-qapa-repository)
  - [Generate 3'UTR annotation using custom BED file of polyA sites](#generate-3utr-annotation-using-custom-bed-file-of-polya-sites)
- [Troubleshooting (TODO)](#troubleshooting-todo)
- [Contributing (TODO)](#contributing-todo)
- [License](#license)
- [TODOs](#TODOs)

## Prerequisites and Installation

Before using this pipeline, make sure you have the following software and tools installed:

- [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- [Singularity](https://docs.sylabs.io/guides/latest/user-guide/quick_start.html#quick-installation-steps)

Then clone this repository to your local machine:

```bash
git clone https://github.com/frattalab/QAPA_snakemake
cd QAPA_snakemake
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

#### TODO: ADD TEST DATA SO THIS CAN BE RUN

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
- if providing input files not located relative to the execution directory, you may need to mount the corresponding directory to the Singularity container using `--singularity-args "--bind /path/to/dir"`

This should take at most a few mins (on my bog standard laptop running Windows Subsystem for Linux)

[Return to table of contents](#table-of-contents)

### Output (TODO)

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

Defines input samples and their associated metadata. An example is provided at [`config/samplesheet_example.csv`](config/samplesheet_example.csv). the minimal, mandatory columns are:

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

To prevent unnecessary computation on successive runs using the same reference annotation, it is possible to provide pre-computed Salmon indices and reference files and skip straight to Salmon quantification:

- set `use_precomputed_salmon_index` to True
- set `salmon_index_dir` to the directory containing the Salmon index
- (Optionally) set `tx2apa` to the table mapping transcript IDs in the qapa build 3'UTR BED with APA_IDs (output at `<main_output_dir>/annotation/qapa_annotation.tx2apa.tsv`)

## Configuring alternative run modes (TODO)

### Use annotation files provided with QAPA repository

### Generate 3'UTR annotation using custom BED file of polyA sites

## Troubleshooting (TODO)

If you encounter any issues or have questions about the pipeline, please check the Issues section of this repository. If you don't find a solution, feel free to open a new issue.

## Contributing (TODO)

If you'd like to contribute to this pipeline, please follow the guidelines in CONTRIBUTING.md.

## License
This pipeline is licensed under the MIT License.

## TODOs

- Example/test data
- Documentation
  - Example output & directory structure
  - Document output files (specific to pipeline)
  - Alternative run modes
    - Using QAPA provided annotation files
    - Generate 3'UTR annotation with custom BED file
- Add scripts & steps to compute summarised PAU & LABRAT psi values
- Pull request to main QAPA repository where applicable (e.g. genome decoys, bug fix)

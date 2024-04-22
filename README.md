# Snakemake workflow: `QAPA_snakemake`


[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)

## Overview

A Snakemake workflow for running QAPA for 3'UTR alternative polyA site quantification, plus additional steps to generate isoform count matrices & perform differential polyA site usage analysis with DEXSeq.

The workflow currently has the following functionality: (along with custom extensions to QAPA)

- Build QAPA annotation files and Salmon indices from scratch for given input reference annotation
  - Optionally add genome sequence (Citation) as decoys to transcriptome index
- Use pre-provided annotations/indices (e.g. from a previous run/QAPA repository) to run Salmon quantification & qapa quant
- Generate isoform count matrices using tximport
- Perform differential polyA site usage with DEXSeq

NOTE: This workflow currently uses a forked version of QAPA available at x. The main changes are:

- Option to add genome sequence as a decoys to transcriptome index built by Salmon.
- Attempt at addressing [issue 49](https://github.com/morrislab/qapa/issues/49)
- Add penultimate exons to transcript models quantified by Salmon (inspired by LABRAT)

The above features are completely optional (with exception of bug fix attempt), so a standard QAPA workflow can be configured if desired.

## Table of Contents

- [Prerequisites and Installation](#prerequisites-and-installation)
- [Configuration](#configuration)
- [Usage](#usage)
  - [Dry Run](#dry-run)
  - [Running the Pipeline](#running-the-pipeline)
- [Results](#results)
- [Reusing generated references](#reusing-generated-references)
- [Configuring alternative run modes](#configuring-alternative-run-modes)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)

## Prerequisites and Installation

Before using this pipeline, make sure you have the following software and tools installed:

- [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- [Singularity](https://docs.sylabs.io/guides/latest/user-guide/quick_start.html#quick-installation-steps)

Clone this repository to your local machine:

```bash
git clone https://github.com/yourusername/your-repo.git
cd your-repo
```

## Configuration

Edit the config.yaml file to customize the pipeline configuration. This file contains parameters such as input file paths, output directories, and other settings to configure the 'run-mode' of QAPA. Make sure to configure it according to your specific needs.

**Sample Table**: defines input samples and their associated metadata. An example is provided at config/samplesheet_example.csv. the minimal columns are:

- `sample_name` - unique identifier for a given sample
- `condition` - key to group samples belonging to the same experimental condition.
- `fastq1` - path to first mate of paired-end sequencing run in gzipped FASTQ format. If sample is single-end, fill this column with the path to the gzipped FASTQ file
- `fastq2` -  path to second mate of paired-end sequencing run in gzipped FASTQ format. If sample is single-end, leave this column blank/empty

(Note: it is possible to pass a mixture of paired-end and single-end samples to the pipeline)

### Differential usage

Differential PAS usage with DEXSeq can be toggled by setting the `run_differential_apa` option in the config file to True/False.

If performing differential isoform usage using DEXSeq, you will also need to construct '**contrasts**' table, which specifies the group-wise comparisons you'd like to perform. An example is provided at `config/contrasts_example.csv`. The minimal columns are:

- `contrast_name` - unique identifier for the comparison of interest (used to differentiate output tables)
- `column_name` - name of column in sample table from which to extract group identifiers. There is no real need to change this from 'condition' in normal circumstances
- `base_key` - group identifier in 'column_name' column of sample table corresponding to the 'base'/'control' (i.e. denominator) condition in the comparison
- `contrast_key` - group identifier in 'column_name' column of sample table corresponding to the 'contrast'/'treatment' (i.e. numerator) condition in the comparison

Note: Fold changes, mean and median delta PAU values are calculated with contrast/base or contrast - base.

If wish to perform differential isoform usage accounting for additional covariates, you can append these columns to the sample table. You should subsequently update `dexseq_formula_full` and `dexseq_formula_reduced` in the config file to include the additional interaction terms. See DEXSeq documentation for full details on how to do this. A simple example below demonstrates how to construct the formulas if you wanted to add `batch` as an additional covariate in the model.

```
# simple comparison of exon usage with respect to condition

dexseq_formula_full: "~ sample_name + exon + condition:exon"
dexseq_formula_reduced: "~ sample_name + exon"

# additionally account for 'batch' variable in both full and null model, so any differences in exon usage due to batch are considered in both models (and are only testing for differences in exon usage that can be accounted for by 'condition') 
# Note the 'batch' column should be appended to the sample table

dexseq_formula_full: "~ sample_name + exon + batch:exon + condition:exon"
dexseq_formula_reduced: "~ sample_name + exon + batch:exon"
```

## Usage

### Dry Run

Before executing the pipeline, it's a good practice to perform a dry run to ensure everything is set up correctly. This allows you to see the execution plan without actually running the commands.

To perform a dry run with your custom configuration file (e.g., config/config.QAPA.yaml), run the following command:

```bash
snakemake --use-singularity --configfile path/to/my_config.yaml --dryrun
```

Replace path/to/my_config.yaml with the actual path to your custom configuration file. This will generate a list of tasks that would be executed without actually executing them.

### Running the Pipeline

To run the pipeline for real with Singularity containers, use the following command:

```bash

snakemake --use-singularity -p --singularity-args="-b /dir/on/machine" --configfile path/to/my_config.yaml
```

Replace path/to/my_config.yaml with the actual path to your custom configuration file. Similarly, it is highly recommended to change /dir/on/machine/ to a directory 

## Results
The pipeline will generate output files and results in the specified output directories as configured in the config.yaml file.

## Reusing generated references

To prevent unnecessary computation on successive runs using the same reference annotation, it is possible to provide pre-computed Salmon indices and reference files for a given run and skip straight to Salmon quantification. This can be achieved by 


## Configuring alternative run modes

## Troubleshooting

If you encounter any issues or have questions about the pipeline, please check the Issues section of this repository. If you don't find a solution, feel free to open a new issue.

## Contributing (TODO)

If you'd like to contribute to this pipeline, please follow the guidelines in CONTRIBUTING.md.

## License
This pipeline is licensed under the MIT License.
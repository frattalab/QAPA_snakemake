# Snakemake workflow: `QAPA_snakemake`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/<owner>/<repo>/workflows/Tests/badge.svg?branch=main)](https://github.com/<owner>/<repo>/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for running QAPA for 3'UTR alternative polyA site quantification and additional steps to perform differential polyA site usage analysis with DEXSeq.

The workflow currently has the following functionality:
- Build QAPA annotation files and Salmon indices from scratch for given input reference annotation
- Use pre-provided QAPA annotation files and run Salmon quantification & qapa quant
- Use annotations/indices  (e.g. from a previous run) to run Salmon quantification & qapa quant
- Perform differential polyA site usage with DEXSeq


## Usage

This workflow uses singularity to manage dependencies


The workflow was was designed to run all necessary steps to run QAPA on input data. Since it is a reference based tool, these inputs only need to be generated once for a given reference dataabse 


# Use case




The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=<owner>%2F<repo>).




If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) <repo>sitory and its DOI (see above).

# TODO

* Replace `<owner>` and `<repo>` everywhere in the template (also under .github/workflows) with the correct `<repo>` name and owning user or organization.
* Replace `<name>` with the workflow name (can be the same as `<repo>`).
* Replace `<description>` with a description of what the workflow does.
* The workflow will occur in the snakemake-workflow-catalog once it has been made public. Then the link under "Usage" will point to the usage instructions if `<owner>` and `<repo>` were correctly set.
* Remove mentions of conda environments (not used at the mo)
* Add steps to compute mean usages, delta usages (Roza)
* Add tximport steps
* Add saturn
* Stop printing sample table in Snakefile

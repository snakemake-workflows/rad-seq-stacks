# Snakemake workflow: rad-seq-stacks

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.2.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/koesterlab/rad-seq-stacks.svg?branch=master)](https://travis-ci.org/koesterlab/rad-seq-stacks)

This is the template for a new Snakemake workflow. Replace this text with a comprehensive description covering the purpose and domain.
Insert your code into the respective folders, i.e. `scripts`, `rules`, and `envs`. Define the entry point of the workflow in the `Snakefile` and the main configuration in the `config.yaml` file.

## Authors

* Johannes Köster (@johanneskoester)
* Henning Timm (@HenningTimm)

## Usage

### Step 1: Install workflow

If you simply want to use this workflow, download and extract the [latest release](https://github.com/snakemake-workflows/rad-seq-stacks/releases).
If you intend to modify and further develop this workflow, fork this repository. Please consider providing any generally applicable modifications via a pull request.

In any case, if you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this repository and, if available, its DOI (see above).

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`.

### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake -n --use-conda

Execute the workflow locally via

    snakemake --cores $N --use-conda

using `$N` cores or run it in a cluster environment via

    snakemake --cluster qsub --jobs 100 --use-conda

or

    snakemake --drmaa --jobs 100 --use-conda

See the [Snakemake documentation](https://snakemake.readthedocs.io) for further details.

### Step 4: Create report

To enclose all results in a single, portable HTML report, run

    snakemake --report
    
in the working directory.

## Testing

Tests cases are in the subfolder `.test`. They should be executed via continuous integration with Travis CI.

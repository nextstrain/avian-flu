# NCBI Ingest

This workflow ingests public data from NCBI, specifically for the H5N1 outbreak
in the U.S. in 2024.

## Workflow Usage

The workflow can be run from the top level pathogen repo directory:
```
nextstrain build ingest ingest_ncbi --configfile build-configs/ncbi/defaults/config.yaml
```

Alternatively, the workflow can also be run from within the ingest directory:
```
cd ingest
nextstrain build . ingest_ncbi --configfile build-configs/ncbi/defaults/config.yaml
```

This produces the default outputs of the NCBI ingest workflow:

- metadata      = ncbi/results/metadata.tsv
- sequences     = ncbi/results/sequences.fasta


## Defaults

The defaults directory contains all of the default configurations for the NCBI ingest workflow.

[defaults/config.yaml](defaults/config.yaml) contains all of the default configuration parameters
used for the ingest workflow. Use Snakemake's `--configfile`/`--config`
options to override these default values.

## Snakefile and rules

The rules directory contains separate Snakefiles (`*.smk`) as modules of the core ingest workflow.
The modules of the workflow are in separate files to keep the main ingest [Snakefile](Snakefile) succinct and organized.

The `workdir` is hardcoded to be the ingest directory so all filepaths for
inputs/outputs should be relative to the ingest directory.

Modules are all [included](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#includes)
in the main Snakefile in the order that they are expected to run.

# Ingest

This workflow ingests public data from NCBI and outputs curated metadata and
sequences that can be used as input for the phylogenetic workflow.

If you have another data source or private data that needs to be formatted for
the phylogenetic workflow, then you can use a similar workflow to curate your
own data.

## Workflow Usage

The workflow can be run from the top level pathogen repo directory:
```
nextstrain build ingest
```

Alternatively, the workflow can also be run from within the ingest directory:
```
cd ingest
nextstrain build .
```

This produces the default outputs of the ingest workflow:

- metadata      = results/metadata.tsv
- sequences     = results/sequences.fasta

### Dumping the full raw metadata from NCBI Datasets

The workflow has a target for dumping the full raw metadata from NCBI Datasets.

```
nextstrain build ingest dump_ncbi_dataset_report
```

This will produce the file `ingest/data/ncbi_dataset_report_raw.tsv`,
which you can inspect to determine what fields and data to use if you want to
configure the workflow for your pathogen.

## Defaults

The defaults directory contains all of the default configurations for the ingest workflow.

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

### Nextclade

Nextstrain is pushing to standardize ingest workflows with Nextclade runs to include Nextclade outputs in our publicly
hosted data. However, if a Nextclade dataset does not already exist, it requires curated data as input, so we are making
Nextclade steps optional here.

If Nextclade config values are included, the Nextclade rules will create the final metadata TSV by joining the Nextclade
output with the metadata. If Nextclade configs are not included, we rename the subset metadata TSV to the final metadata TSV.

To run Nextclade rules, include the `defaults/nextclade_config.yaml` config file with:

```
nextstrain build ingest --configfile defaults/nextclade_config.yaml
```

> [!TIP]
> If the Nextclade dataset is stable and you always want to run the Nextclade rules as part of ingest, we recommend
moving the Nextclade related config parameters from the `defaults/nextclade_config.yaml` file to the default config file
`defaults/config.yaml`.

## Build configs

The build-configs directory contains custom configs and rules that override and/or
extend the default workflow.

- [nextstrain-automation](build-configs/nextstrain-automation/) - automated internal Nextstrain builds.


## Vendored

This repository uses [`git subrepo`](https://github.com/ingydotnet/git-subrepo)
to manage copies of ingest scripts in [vendored](vendored), from [nextstrain/ingest](https://github.com/nextstrain/ingest).

See [vendored/README.md](vendored/README.md#vendoring) for instructions on how to update
the vendored scripts.

# nextstrain.org/avian-flu

This is the Nextstrain build for avian influenza subtypes A/H5N1, A/H5NX, A/H7N9, and A/H9N2 as well as analysis of the 2024 A/H5N1 cattle-flu outbreak.
The most up-to-date builds of avian influenza can be found [on nextstrain.org](https://nextstrain.org/avian-flu).
Please see [nextstrain.org/docs](https://nextstrain.org/docs) for details about augur and pathogen builds.

The Snakemake pipeline is parameterised by two config files, one for the A/H5N1, A/H5NX, A/H7N9, and A/H9N2 builds and one for the 2024 A/H5N1 cattle-flu outbreak.


## Segment-level GISAID builds

The `config/gisaid.yaml` config builds 32 Auspice datasets (8 segments x 4 subtypes (A/H5N1, A/H5NX, A/H7N9, A/H9N2)) using GISAID data and can be run via

```bash
snakemake --cores 1 -pf --configfile config/gisaid.yaml
```

This pipeline starts by downloading data from a private S3 bucket and the appropriate credentials are required; see below for how to use locally ingested files.
For rapid AWS rebuild run as:


```bash
nextstrain build --aws-batch --aws-batch-cpus 16 --aws-batch-memory 28800 . --jobs 16 --configfile config/gisaid.yaml
```

Please see [nextstrain.org/docs](https://nextstrain.org/docs) for details about augur and pathogen builds.

### Deploying builds

The pipeline can automatically deploy resulting builds within the auspice folder
to nextstrain.org by running:

```bash
nextstrain build . --configfile config/gisaid.yaml -f deploy_all
```

## H5N1 Cattle Outbreak (2024)

We produce per-segment and whole-genome (concatenated segments) builds for the ongoing H5N1 cattle-flu outbreak.
These use NCBI data including consensus genomes and SRA data assembled via the Andersen lab's [avian-influenza repo](https://github.com/andersen-lab/avian-influenza).

> Running this build will overwrite GISAID files in `./data` and thus you can't maintain or run GISAID & NCBI builds in parallel. In most cases this isn't an issue and [we are working on improving this](https://github.com/nextstrain/avian-flu/issues/70). You may want to proactively remove the `./data` directory yourself to make sure everything works as expected.


```bash
snakemake --cores 1 -pf --configfile config/h5n1-cattle-outbreak.yaml
```

This pipeline starts by downloading data from a public S3 bucket.


**Genome builds**

The build is restricted to a set of strains where we think there's no reassortment, with outgroups excluded (`config/dropped_strains_h5n1-cattle-outbreak.txt`).
Output files will be placed in `results/h5n1-cattle-outbreak/genome`.
See `Snakefile.genome` for more details.


**Segment-level builds**

Strains for each segment are chosen by first constructing a general tree for the segment with all strains from 2024 onwards and then taking the clade which contains all strains in the genome build.
This should allow any reassortments to be highlighted and will also include outbreak strains which are missing from the genome build (because they don't have all 8 segments sequenced).

> Note that generating any segment-level build here will necessarily build the genome tree, as it's needed to identify the clade of interest in each segment.

## H5N1 D1.1 Cattle outbreak (2025)

> This build is a work in progress and relies on recent improvements to ingest which add GenoFLU constellations to the metadata TSV

The H5N1-D1.1/genome build uses a similar approach to the H5N1 Cattle Outbreak above however samples all data to a GenoFLU constellation of D1.1.

```bash
snakemake --cores 1 -pf --configfile config/h5n1-d1.1.yaml
```


## Creating a custom build
The easiest way to generate your own, custom avian-flu build is to use the quickstart-build as a starting template. Simply clone the quickstart-build, run with the example data, and edit the Snakefile to customize. This build includes example data and a simplified, heavily annotated Snakefile that goes over the structure of Snakefiles and annotates rules and inputs/outputs that can be modified. This build, with it's own readme, is available [here](https://github.com/nextstrain/avian-flu/tree/master/quickstart-build).

## Features unique to avian flu builds

#### cleavage site annotations
Influenza virus HA is translated as a single peptide (HA0) that is cleaved to form the mature, functional form (HA1 and HA2). In all human strains and many avian strains, the cleavage site is composed of a single, basic amino acid residue. However, some avian influenza subtypes, particularly H5s, have acquired additional basic residues immediately preceding the HA cleavage site. In some cases, this results in addition of a furin cleavage motif, allowing HA to be cleaved by furin, which is ubiquitously expressed, and allows for viral replication across a range of tissues. The addition of this "polybasic cleavage site" is one of the prime determinants of avian influenza virulence. In these builds, we have annotated whether strains contain a furin cleavage motif, defined here as the sequence `R-X-K/R-R` immediately preceding the start of HA2, where `X` can be any amino acid. We have also added a color by for the cleavage site sequence, which we define here as the 4 bases preceding HA2.

#### clade labeling
H5 viruses are classified into clades, which are currently viewable as a color by on [nextstrain.org](https://nextstrain.org/avian-flu/h5n1/ha?c=h5_label_clade). Because clade annotations are not available in all public databases, we annotate sequences with their most likely clade using a tool developed by Samuel S. Shepard at CDC called [LABEL](https://wonder.cdc.gov/amd/flu/label/). The assigned clade for each H5N1 or H5Nx sequence are available as public tsvs [here](https://github.com/nextstrain/avian-flu/tree/master/clade-labeling).

To update the `clades.tsv` file with clade annotations for new sequences, run:

`snakemake -s Snakefile.clades -p --cores 1`

To run the builds on without updating the clades file, run:

`snakemake -p --cores 1`

To string these together and update the `clades.tsv` file for new sequences and then run the builds:

`snakemake -s Snakefile.clades -p --cores 1 && snakemake -p --cores 1`

#### Using the same strains for all segments

By default we subsample data for each segment independently.
Alternatively, you can ask the pipeline to use the same strains for each segment.
This modifies the pipeline in a few ways:
1. An additional metadata processing step is added which counts the number of segments a strain has sequence data for
2. Subsampling is performed as normal for the HA segment, with the additional condition that each sample has sequence data for all segments
3. All other segments are subsampled to contain the same strains as used for HA in step 2

To enable this set the config parameter `same_strains_per_segment` to a truthy value. If you are calling `snakemake` directly you can add

```bash
--config 'same_strains_per_segment=True'
```

If you are using `nextstrain build` then add that to the end of the command (i.e. as a parameter which will be passed through to Snakemake).

Note that you may need to remove any existing data in `results/` in order for snakemake to correctly regenerate the intermediate files.

#### Using locally ingested data (instead of downloading from S3)

If you have run the fauna ingest pipeline locally, you will have files such as
```
ingest/fauna/results/
├── metadata.tsv
├── sequences_ha.fasta
├── sequences_mp.fasta
├── sequences_na.fasta
├── sequences_np.fasta
├── sequences_ns.fasta
├── sequences_pa.fasta
├── sequences_pb1.fasta
└── sequences_pb2.fasta
```

Running with `--config 'inputs=[{"name": "local-ingest", "metadata": "ingest/fauna/results/metadata.tsv", "sequences": "ingest/fauna/results/sequences_{segment}.fasta"}]'` will switch the default S3 source to the local paths listed above.
Alternatively, you can use a config overlay YAML and supply it via `--configfile`:

```yaml
inputs:
  - name: local-ingest
    metadata: ingest/fauna/results/metadata.tsv
    sequences: ingest/fauna/results/sequences_{segment}.fasta
```


For the **h5n1-cattle-outbreak** the approach is similar, but you'll want to replace "fauna" with "joined-ncbi", "ncbi" or "andersen-lab", as appropriate.

## To modify this build to work with your own data
Although the simplest way to generate a custom build is via the quickstart build, you are welcome to clone this repo and use it as a starting point for running your own, local builds if you'd like. The [Nextstrain docs](https://docs.nextstrain.org/en/latest/index.html) are a fantastic resource for getting started with the Nextstrain pipeline, and include some [great tutorials](https://docs.nextstrain.org/en/latest/install.html) to get you started. This build is slightly more complicated than other builds, and has a few custom functions in it to accommodate all the features on [nextstrain.org](https://nextstrain.org/avian-flu), and makes use of wildcards for both subtypes and gene segments. If you'd like to adapt this full, non-simplified pipeline here to your own data (which you may want to do if you also want to annotate clades), you would need to make a few changes and additions:


### Use additional metadata and/or sequences

Using a config YAML provided via `--configfile`, or the same data via `--config` you can supply additional inputs.
For instance, this repo has a small set of example metadata + HA sequences.
We could add these to the default inputs via:

```yaml
additional_inputs:
  - name: example-data
    metadata: example_data/gisaid/metadata.tsv
    sequences:
        ha: example_data/gisaid/sequences_ha.fasta
```

Which will merge `example_data/gisaid/metadata.tsv` with the default metadata, and add the sequences from `example_data/gisaid/sequences_ha.fasta` to the default HA sequences (all other segments will just use the default sequences).

If you had sequences for each segment you could use a wildcard for the segment like so:

```yaml
additional_inputs:
  - name: example-data
    metadata: example_data/gisaid/metadata.tsv
    sequences: example_data/gisaid/sequences_{segment}.fasta
```


> NOTE: These added data will be subject to the same filtering rules as the starting data.
  At a minimum, you'll want to ensure new sequences have metadata which defines their subtype and date, as filter steps will prune out those without valid values here.

> NOTE: Metadata merging is via `augur merge` and we add a column per input source to indicate the origin of data.
  For instance the above example would have a `input_example-data` column with values of `1` for metadata rows which were included in `example_data/gisaid/metadata.tsv` and `0` otherwise.
  You can use this for additional filtering commands as needed.


### Replace the starting metadata and sequences

Using the same approach as above but replacing `additional_inputs` with `inputs` will replace the default inputs.
See "Using locally ingested data (instead of downloading from S3)" (above) for an example of this.


### Run a single H5N1 HA all-time analysis

There is a `test_target` rule which will produce a single H5N1 HA "all-time" tree.
You can combine this with the example data (see above) to run a single build using a small number of sequences:


``` bash
snakemake --cores 2 \
    --configfile config/gisaid.yaml \
    --config 'inputs=[{"name": "example", "metadata": "example_data/gisaid/metadata.tsv", "sequences": "example_data/gisaid/sequences_{segment}.fasta"}]' \
    -pf test_target
```

### clade labeling

If you'd like to run clade labeling, you will need to install [LABEL](https://wonder.cdc.gov/amd/flu/label/) yourself. This pipeline assumes that LABEL is located in `avian-flu/flu-amd/`, and should work if you install it into the `avian-flu` directory. If you do not need to label clades, then you can delete `rule add_h5_clade`, the function `metadata_by_wildcards`. You will need to make sure that all references to `metadata` in the pipeline are referencing `metadata_subtype_segment`, not `metadata-with-clade_subtype_segment`, which is generated by `rule add_h5_clade` and adds a new column to the metadata file with clade information.

# nextstrain.org/avian-flu

This is the Nextstrain build for avian influenza subtypes A/H5N1, A/H5NX, A/H7N9, and A/H9N2 as well as analysis of the 2024 A/H5N1 cattle-flu outbreak.
The most up-to-date builds of avian influenza can be found [on nextstrain.org](https://nextstrain.org/avian-flu).
Please see [nextstrain.org/docs](https://nextstrain.org/docs) for details about augur and pathogen builds.

The Snakemake pipeline is parameterised by two config files, one for the A/H5N1, A/H5NX, A/H7N9, and A/H9N2 builds and one for the 2024 A/H5N1 cattle-flu outbreak.

### Table of Contents

* [Running segment-level GISAID builds](#running-segment-level-gisaid-builds)
* [Running H5N1 Cattle Outbreak (2024) builds](#running-h5n1-cattle-outbreak-2024-builds)
* [Creating a custom build via config overlays](#creating-a-custom-build-via-config-overlays)
* [Running builds in a separate working directory](#running-builds-in-a-separate-working-directory)
* [Adding private data](#adding-private-data)
* [Cleavage Site Annotations](#cleavage-side-annotations)
* [H5 Clade Labeling](#h5-clade-labeling)
* [Other build customisations](#other-build-customisations)


## Running segment-level GISAID builds

The `config/gisaid.yaml` config builds 48 Auspice datasets (8 segments x 4 subtypes (A/H5N1, A/H5NX, A/H7N9, A/H9N2) x 1-2 time resolutions (all-time, 2y)) using GISAID data and can be run via

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

## Running H5N1 Cattle Outbreak (2024) builds

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


## Creating a custom build via config overlays

The two config files introduced above are designed to be extended (customised) by overlaying your own YAML configurations.
The aim is to allow easy customisation of the workflow outputs (e.g. which subtypes / segments / time-frames to run) and the stopping points (Auspice JSONs or particular intermediate files) via such overlays.
Additionally, modification of parameters (e.g. clock rates, minimum lengths, subsampling criteria, DTA metadata) is possible through overlays without needing to modify the underlying base config or Snakemake pipeline itself.

Config overlays allow you to essentially maintain one or more modifications of the workflow for your own purposes.
For instance you may want a way to easily run only H5N1 HA+NA 2y builds, and a config overlay can achieve this.
Using an overlay keeps this change separate from the config used for Nextstrain automation and also avoids `git` conflicts emerging over time.
When combined with running in a separate working directory (explained below) this becomes even more powerful.

As an example we'll create a custom config `config_local.yaml` which you can then add to the build command described above, e.g. `nextstrain build . --configfile config/gisaid.yaml config_local.yaml`.
We'll use the GISAID config as the example we're extending, however the concepts are the similar for the H5N1 cattle outbreak config.


### Restrict which builds we want to produce

By default we produce 48 Auspice JSONs (4 subtypes * 8 segments * 1-2 time resolutions). We can restrict these by redefining the builds in our config overlay. For instance the following will produce only 2 datasets, `h5n1/ha/2y` and `h5n1/na/2y`:
```yaml
builds:
  - subtype: h5n1
    segment:
      - ha
      - na
    time: 2y
```
Here `builds` is an array of sub-configs, each of which define a combination of subtype, segment and time parameters. Each subtype, segment and time can be a single string (as subtype and time are, above) or an array (as segment is, above).

### Only produce certain intermediate files, not Auspice datasets

By default the "target" for each build is the Auspice JSON (`auspice/avian-flu_h5n1_ha_2y.json` etc) however we can change this if we just want certain intermediate files. Adding the following to the config overlay will stop once we've filtered to the metadata & sequences we would use for each tree

```yaml
target_patterns:
  - "results/{subtype}/{segment}/{time}/metadata.tsv"
  - "results/{subtype}/{segment}/{time}/sequences.tsv"
```

(This works in combination with the custom `builds` definition, above).

### Change the target number of sequences per build

In the `config/gisaid.yaml` this parameter is defined as

```yaml
filter: 
  target_sequences_per_tree:
    "*/*/*": 3000
```

Where the `"*/*/*"` syntax is slash-separated matchers for subtype, segment and time-resolution, and the `*` character means "match everything".
So here we're saying "for every subtype, for every segment, for every time-resolution target 3000 sequences".

If we want to change our h5n1 builds to instead have 5000 sequences (whilst keeping the rest at 3000) we could add the following to our config overlay:
```yaml
filter: 
  target_sequences_per_tree:
    "h5n1/*/*": 5000
```
And since `"h5n1/*/*"` is more specific than `"*/*/*"` it'll take precedence when `subtype="h5n1"`.
(Internally, Snakemake merges these configs together resulting in `'target_sequences_per_tree': {'h5n1/*/*': 5000, '*/*/*': 3000}`. For each combination of subtype/segment/time values within the workflow we consult these dictionaries and pick the most specific match.)


This syntax is concise but powerful, for instance we can parameterise the builds like so:
```yaml
filter: 
  target_sequences_per_tree:
    "*/*/all-time": 5000 # target 5k sequences for the all-time builds
    "h5n1/*/all-time": 10000 # but for h5n1 all-time builds target 10k sequences (this is more specific as it only includes one '*' character)
    "*/*/*": 1000 # target 1k sequences for any other builds (i.e. the 2y builds)
```



### Change other parameters

The pattern introduced in the preceeding section generally applies for all parameters in the workflow.
By reading through the base config YAML (e.g. `config/gisaid.yaml`) you should be able to learn the config structure and add then modify that within your overlay config as needed.

If the parameter is not exposed via the config YAML and you find yourself modifying the underlying Snakefile consider exposing it via the config so that this customisation then becomes available to config overlays.


## Running builds in a separate working directory

> Note: This section doesn't yet work using AWS / docker runtimes but support is forthcoming.

We can run the worfklow from an analysis directory separate to this repo. Imagine the following directory structure:

```
├── avian-flu
│   ├── README.md
│   ├── Snakefile
│   └── ... (lots more files here)
└── analysis_dir
    └── config.yaml
```

Where `config.yaml` is the your config overlay, introduced in the previous section as `config_local.yaml`.

> TODO XXX following uses a different invocation - this is all in flux

From within the `analysis_dir` we can run `snakemake --snakefile ../avian-flu/gisaid/Snakefile --cores 1` and we'll automatically run the workflow using your `config.yaml` config overlay.
This keeps the workflow outputs separate and isolated from the workflow itself.
Depending on how you run builds this can be very liberating; for instance if you are switching between versions of the workflow you can maintain different analysis directories for each version, or if you are running the workflow at different times you can separate the analyses for easy before/after comparisons.

> You don't have to name it `config.yaml`, but if you use a different filename you'll have to specify it via `--configfile <filename>`.

## Adding private data

Private metadata and/or sequences can be supplied by defining `additional_inputs` in a config overlay YAML.

```yaml
additional_inputs:
  - name: secret
    metadata: secret.tsv
    sequences: 
      ha: secret_ha.fasta
```

The filenames here can be S3 URIs (ensure you have credentials set in your environment) or local files.
In this case local files should be specified relative to the analysis directory (typically where you run the command from).

If you have data for all segments you can use a slightly different and more concise syntax:
```yaml
additional_inputs:
  - name: secret
    metadata: secret.tsv
    sequencs: secret_{segment}.fasta
```

> These added data will be subject to the same filtering rules as the starting data.
  At a minimum, you'll want to ensure new sequences have metadata which defines their subtype and date, as filter steps will prune out those without valid values here.

> Metadata merging is via `augur merge` and we add one-hot columns for each input source with the column name `input_{NAME}`, for instance the above example would have a `input_secret` column with values of `1` for metadata rows which were included in `secret.tsv` and `0` otherwise.
  You can use this for additional filtering commands as needed.

By default each workflow config defines a single metadata input and one FASTA per segment.


## Cleavage Site Annotations

Influenza virus HA is translated as a single peptide (HA0) that is cleaved to form the mature, functional form (HA1 and HA2). In all human strains and many avian strains, the cleavage site is composed of a single, basic amino acid residue. However, some avian influenza subtypes, particularly H5s, have acquired additional basic residues immediately preceding the HA cleavage site. In some cases, this results in addition of a furin cleavage motif, allowing HA to be cleaved by furin, which is ubiquitously expressed, and allows for viral replication across a range of tissues. The addition of this "polybasic cleavage site" is one of the prime determinants of avian influenza virulence. In these builds, we have annotated whether strains contain a furin cleavage motif, defined here as the sequence `R-X-K/R-R` immediately preceding the start of HA2, where `X` can be any amino acid. We have also added a color by for the cleavage site sequence, which we define here as the 4 bases preceding HA2.

## H5 Clade Labeling
H5 viruses are classified into clades, which are currently viewable as a color by on [nextstrain.org](https://nextstrain.org/avian-flu/h5n1/ha?c=h5_label_clade). Because clade annotations are not available in all public databases, we annotate sequences with their most likely clade using a tool developed by Samuel S. Shepard at CDC called [LABEL](https://wonder.cdc.gov/amd/flu/label/). The assigned clade for each H5N1 or H5Nx sequence are available as public tsvs [here](https://github.com/nextstrain/avian-flu/tree/master/clade-labeling).

To update the `clades.tsv` file with clade annotations for new sequences, run:

`snakemake -s Snakefile.clades -p --cores 1`

To run the builds on without updating the clades file, run:

`snakemake -p --cores 1`

To string these together and update the `clades.tsv` file for new sequences and then run the builds:

`snakemake -s Snakefile.clades -p --cores 1 && snakemake -p --cores 1`

## Other build customisations

### Using the same strains for all segments

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

### Using locally ingested data (instead of downloading from S3)

Each workflow defines an input via `config['inputs']`, e.g. the GISAID build uses:
```yaml
inputs:
  - name: gisaid
    metadata: s3://nextstrain-data-private/files/workflows/avian-flu/metadata.tsv.zst
    sequences: s3://nextstrain-data-private/files/workflows/avian-flu/{segment}/sequences.fasta.zst
```

You can use an approach analogous to the [addition of private data](#adding-private-data) above to replace this array in your config overlay to point to local files instead.
If you have run the default ingest pipeline a config overlay of

```yaml
inputs:
  - name: ingest
    metadata: ingest/fauna/results/metadata.tsv
    sequences: ingest/fauna/results/sequences_{segment}.fasta
```

Will result in the default inputs being replaced by paths to these local ingest files.
The search order for local files is:
  1. Relative to the analysis directory
  2. Relative to the entry snakefile
  3. Relative to the `avian-flu` directory


### To modify this build to work with your own data
Although the simplest way to generate a custom build is via the quickstart build, you are welcome to clone this repo and use it as a starting point for running your own, local builds if you'd like. The [Nextstrain docs](https://docs.nextstrain.org/en/latest/index.html) are a fantastic resource for getting started with the Nextstrain pipeline, and include some [great tutorials](https://docs.nextstrain.org/en/latest/install.html) to get you started. This build is slightly more complicated than other builds, and has a few custom functions in it to accommodate all the features on [nextstrain.org](https://nextstrain.org/avian-flu), and makes use of wildcards for both subtypes and gene segments. If you'd like to adapt this full, non-simplified pipeline here to your own data (which you may want to do if you also want to annotate clades), you would need to make a few changes and additions:


#### 1. Data is stored on a private S3 bucket

The phylogenetics pipeline starts by downloading data from a private S3 bucket.
If you don't have credentials for this bucket you can run the build using the example data provided in this repository.
Before running the build, copy the example sequences and metadata into the `data/` directory like so:

```bash
mkdir -p data/
cp -r -v example_data/* data/
```

Then run the build with the test target, a H5N1 HA "all time" tree:

``` bash
snakemake test_target
```

If you'd like to consistently run your own data, then you can place your fasta file in `data`. Alternatively, you can alter the `Snakefile` to remove references to S3 and add paths to your own files (see rules `download_sequences` and `download_metadata`).

#### 2. clade labeling
If you'd like to run clade labeling, you will need to install [LABEL](https://wonder.cdc.gov/amd/flu/label/) yourself. This pipeline assumes that LABEL is located in `avian-flu/flu-amd/`, and should work if you install it into the `avian-flu` directory. If you do not need to label clades, then you can delete `rule add_h5_clade`, the function `metadata_by_wildcards`. You will need to make sure that all references to `metadata` in the pipeline are referencing `metadata_subtype_segment`, not `metadata-with-clade_subtype_segment`, which is generated by `rule add_h5_clade` and adds a new column to the metadata file with clade information.

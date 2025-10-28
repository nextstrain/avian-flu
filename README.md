# nextstrain.org/avian-flu

This is the Nextstrain build for avian influenza analyses.
Currently there are two distinct workflows:

**Segment-focused builds** produce segment-level analyses of subtypes A/H5N1, A/H5NX, A/H7N9, and A/H9N2 using GISAID data.

**Genome-focused builds** produce genome-focused analyses of the (ongoing) 2024 A/H5N1 cattle-flu outbreak and the (ongoing) D1.1 outbreaks using NCBI data

The most up-to-date builds of avian influenza can be found [on nextstrain.org](https://nextstrain.org/avian-flu).
Please see [nextstrain.org/docs](https://nextstrain.org/docs) for details about augur and pathogen builds.


## Segment-focused builds (from GISAID data)

**This pipeline starts by downloading data from a private S3 bucket and the appropriate credentials are required;** see below for how to use locally ingested files or NCBI inputs.

The `segment-focused/Snakefile` entrypoint and associated `segment-focused/config.yaml` config builds 48 Auspice datasets (8 segments x 4 subtypes (A/H5N1, A/H5NX, A/H7N9, A/H9N2) with 2 time-ranges (all-time, 2y) for A/H5N1 and A/H5NX) using GISAID data and can be run via

```bash
snakemake --cores 1 -pf --snakefile segment-focused/Snakefile
```

For rapid AWS rebuild run as:

```bash
nextstrain build --aws-batch --aws-batch-cpus 16 --aws-batch-memory 28800 . --jobs 16 --snakefile segment-focused/Snakefile
```

This uses the `segment-focused/config.yaml` config, but you can amend this with additional configs (see below).

### Deploying builds

The pipeline can automatically deploy resulting builds within the auspice folder
to nextstrain.org by running:

```bash
nextstrain build . --snakefile segment-focused/Snakefile -f deploy_all
```

## Genome-focused builds (from NCBI data)

We produce whole-genome focused (and, sometimes, individual segment builds) for the ongoing H5N1 cattle-flu outbreak.
These use NCBI data including consensus genomes and SRA data assembled via the Andersen lab's [avian-influenza repo](https://github.com/andersen-lab/avian-influenza) by default.

Currently we build two sets of analyses, both of which use GenoFLU to partition data into non-reassorting sets.
The `h5n1-cattle-flu` datasets are GenoFLU constellation B3.13 and represent the 2024-ongoing North American cattle outbreak.
The `h5n1-d1.1` datasets are GenoFLU constellation D1.1 and represent the separate but contemporaneous (2025) zoonoses into cattle; as of Feb 2025 there are at least two separate D1.1 jumps into cattle.


```bash
snakemake --cores 1 -pf --snakefile genome-focused/Snakefile
```

This uses the `genome-focused/config.yaml` config, but you can amend this with additional configs (see below).
The pipeline starts by downloading data from a public S3 bucket before filtering the data based on the aforementioned GenoFLU constellations.
(The annotations themselves are added in the separate ingest pipeline.)

The B3.13 datasets also include per-segment analyses, which are also limited to B3.13 strains.
We're aiming to improve the diversity and significance of these datasets soon.

## Running analyses in a separate working directory

We can run the avian-flu worfklows from analysis directory/directories separate to this repo. Imagine the following directory structure:

```
├── avian-flu
│   ├── README.md
│   ├── segment-focused
│   │   ├── Snakefile
│   │   └── config.yaml
│   └── ... (lots more files)
└── analysis_dir
```

From within the `analysis_dir` we can run `snakemake --cores 1 -pf --snakefile ../avian-flu/segment-focused/Snakefile` and we'll automatically run the workflow and write all output files (`data/`, `logs/`, `results/`, `auspice/`) within the `analysis_dir`. This keeps the workflow outputs separate and isolated from the workflow itself.

Depending on how you run builds this can be very liberating; for instance if you are running the workflow at different times you can separate the analyses for easy before/after comparisons. Similarly if you make changes to the workflow you can maintain independent sets of results, similar to git worktrees.

> You can use Snakemake's `--directory` argument to define the analysis directory (working directory) if that's easier

## Creating a custom build via config overlays

Each of the workflows described above is designed to be extended (customised) by overlaying your own YAML configurations.
The aim is to allow easy customisation of the workflow outputs (e.g. which subtypes / segments / time-frames to run) and the stopping points (Auspice JSONs or particular intermediate files) via such overlays.
Additionally, modification of parameters (e.g. clock rates, minimum lengths, subsampling criteria, DTA metadata) is possible through overlays without needing to modify the underlying base config or Snakemake pipeline itself.

Config overlays allow you to essentially maintain one or more modifications of the workflow for your own purposes.
For instance you may want a way to easily run only H5N1 HA+NA 2y builds, and a config overlay can achieve this.
Using an overlay keeps this change separate from the config used for Nextstrain automation and also avoids `git` conflicts emerging over time.
When combined with running in a separate working directory (explained below) this becomes even more powerful.

For the following examples we presume you are running within a separate analysis directory as introduced in the previous section.

We'll start by creating a config (overlay) YAML, `config.yaml` in our analysis directory, and use the segment-focused workflow described above, however the concepts are the similar for the genome-focused workflow.

> You can choose a different name for the file, but if you do you'll have to supply it to the command via `--configfile <filename>` as only `config.yaml` will be automatically detected.

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

In the `segment-focused/config.yaml` this parameter is defined as

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
By reading through the relevant base config YAML (e.g. `segment-focused/config.yaml` or `genome-focused/config.yaml`) you should be able to learn the config structure and add then modify that within your overlay config as needed.

If the parameter is not exposed via the config YAML and you find yourself modifying the underlying Snakefile consider exposing it via the config so that this customisation then becomes available to config overlays.


## Creating a custom build via quickstart pipeline
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

## Changing the input data

The following sections explain how to:

* run the segment-level builds with NCBI (public) data (or the genome builds with GISAID data)
* use locally ingested data (instead of downloading from S3)
* use your own input data either on its own or merged with other inputs

#### Running the segment-level builds using NCBI data

We can use a second config YAML to changes the inputs.
For instance, create a `config_local.yaml` with:
```yaml
inputs:
  - name: ncbi
    metadata: s3://nextstrain-data/files/workflows/avian-flu/h5n1/metadata.tsv.zst
    sequences: s3://nextstrain-data/files/workflows/avian-flu/h5n1/{segment}/sequences.fasta.zst
```
and run the workflow (described above) with `--configfile config_local.yaml`, e.g. `snakemake --cores 1 -pf --snakefile segment-focused/Snakefile --configfile config_local.yaml`

To run the genome-focused builds with GISAID data use the same approach, but using the (private) data endpoints defined in `segment-focused/config.yaml`.

#### Using locally ingested data (instead of downloading from S3)

If you have run one of the ingest pipeline locally, you will have files such as
```
ingest
└── <data-source-name>
   └── results
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

Where `<data-source-name>` may be "fauna", "joined-ncbi" etc. 

Similar to before we can use a config-overlay to switch the starting inputs of the workflow, or append to the existing inputs, as desired.
Create a config overlay YAML and supply it via `--configfile`:

```yaml
inputs:
  - name: local-ingest
    metadata: ingest/fauna/results/metadata.tsv
    sequences: ingest/fauna/results/sequences_{segment}.fasta
```

(replace "fauna" with "joined-ncbi" etc, as appropriate).


#### Using your own input data

If you have your own data you can use the same approach, but rather than the paths being `ingest/fauna/...` etc have them point to your sequences & metadata files.

If you instead want to merge / combine your data with the existing data, we can specify the `additional_inputs` config parameter.
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


### Run a single H5N1 HA all-time analysis

There is an example build which will produce a single H5N1 HA "all-time" tree
using the example data (see above):


``` bash
snakemake --cores 2 \
    --snakefile segment-focused/Snakefile \
    --configfile build-configs/ci/config.yaml
```

### clade labeling

If you'd like to run clade labeling, you will need to install [LABEL](https://wonder.cdc.gov/amd/flu/label/) yourself. This pipeline assumes that LABEL is located in `avian-flu/flu-amd/`, and should work if you install it into the `avian-flu` directory. If you do not need to label clades, then you can delete `rule add_h5_clade`, the function `metadata_by_wildcards`. You will need to make sure that all references to `metadata` in the pipeline are referencing `metadata_subtype_segment`, not `metadata-with-clade_subtype_segment`, which is generated by `rule add_h5_clade` and adds a new column to the metadata file with clade information.

# nextstrain.org/flu/avian

This is the Nextstrain build for avian influenza subtypes A/H5N1, A/H5NX, A/H7N9, and A/H9N2.
The most up-to-date builds of avian influenza can be found [on nextstrain.org](https://nextstrain.org/flu/avian).
Please see [nextstrain.org/docs](https://nextstrain.org/docs) for details about augur and pathogen builds.

## Building

All 32 builds (4 subtypes x 8 segments) can be built by running `snakemake`. For rapid AWS rebuild run as:

```bash
nextstrain build --aws-batch --aws-batch-cpus 16 --aws-batch-memory 28800 . --jobs 16
```

Please see [nextstrain.org/docs](https://nextstrain.org/docs) for details about augur and pathogen builds.

## Creating a custom build 
The easiest way to generate your own, custom avian-flu build is to use the  quickstart-build as a starting template. Simply clone the quickstart-build, run with the example data, and edit the Snakefile to customize. This build includes example data and a simplified, heavily annotated Snakefile that goes over the structure of Snakefiles and annotates rules and inputs/outputs that can be modified. This build, with it's own readme, is available [here](https://github.com/nextstrain/avian-flu/tree/master/quickstart-build).

### Features unique to avian flu builds

#### cleavage site annotations 
Influenza virus HA is translated as a single peptide (HA0) that is cleaved to form the mature, functional form (HA1 and HA2). In all human strains and many avian strains, the cleavage site is composed of a single, basic amino acid residue. However, some avian influenza subtypes, particularly H5s, have acquired additional basic residues immediately preceding the HA cleavage site. In some cases, this results in addition of a furin cleavage motif, allowing HA to be cleaved by furin, which is ubiquitously expressed, and allows for viral replication across a range of tissues. The addition of this "polybasic cleavage site" is one of the prime determinants of avian influenza virulence. In these builds, we have annotated whether strains contain a furin cleavage motif, defined here as the sequence `R-X-K/R-R` immediately preceding the start of HA2, where `X` can be any amino acid. We have also added a color by for the cleavage site sequence, which we define here as the 4 bases preceding HA2. 

#### clade labeling
H5 viruses are classified into clades, which are currently viewable as a color by on [nextstrain.org](https://nextstrain.org/flu/avian/h5n1/ha?c=h5_label_clade). Because clade annotations are not available in all public databases, we annotate sequences with their most likely clade using a tool developed by Samuel S. Shepard at CDC called [LABEL](https://wonder.cdc.gov/amd/flu/label/). The assigned clade for each H5N1 or H5Nx sequence are available as public tsvs [here](https://github.com/nextstrain/avian-flu/tree/master/clade-labeling).

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

### To modify this build to work with your own data
Although the simplest way to generate a custom build is via the quickstart build, you are welcome to clone this repo and use it as a starting point for running your own, local builds if you'd like. The [Nextstrain docs](https://docs.nextstrain.org/en/latest/index.html) are a fantastic resource for getting started with the Nextstrain pipeline, and include some [great tutorials](https://docs.nextstrain.org/en/latest/install.html) to get you started. This build is slightly more complicated than other builds, and has a few custom functions in it to accommodate all the features on [nextstrain.org](https://nextstrain.org/flu/avian), and makes use of wildcards for both subtypes and gene segments. If you'd like to adapt this full, non-simplified pipeline here to your own data (which you may want to do if you also want to annotate clades), you would need to make a few changes and additions:

#### 1. fauna / RethinkDB credentials
This build starts by pulling sequences from our live [fauna](https://github.com/nextstrain/fauna) database (a RethinkDB
instance). This requires environment variables `RETHINK_HOST` and `RETHINK_AUTH_KEY` to be
set.

If you don't have access to our database, you can run the build using the example data
provided in this repository. Before running the build, copy the example sequences into the
`data/` directory like so:

```bash
mkdir data/
cp example_data/* data/
```

Then run the the build. If you'd like to consistently run your own data, then you can place your fasta file in `data`. Alternatively, you can alter the `Snakefile` to remove references to our database and add paths to your own files. To do this, remove `rule download`, add paths to your input data (sequences and metadata) in `rule files`, and add those paths as the input to `rule parse`. 

#### 2. clade labeling 
If you'd like to run clade labeling, you will need to install [LABEL](https://wonder.cdc.gov/amd/flu/label/) yourself. This pipeline assumes that LABEL is located in `avian-flu/flu-amd/`, and should work if you install it into the `avian-flu` directory. If you do not need to label clades, then you can delete `rule add_h5_clade`, the function `metadata_by_wildcards`. You will need to make sure that all references to `metadata` in the pipeline are referencing `metadata_subtype_segment`, not `metadata-with-clade_subtype_segment`, which is generated by `rule add_h5_clade` and adds a new column to the metadata file with clade information. 

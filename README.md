# nextstrain.org/flu/avian

This is the Nextstrain build for avian influenza subtypes A/H5N1, A/H5NX, A/H7N9, and A/H9N2.
The most up-to-date builds of avian influenza can be found [on nextstrain.org](https://nextstrain.org/flu/avian).
Please see [nextstrain.org/docs](https://nextstrain.org/docs) for details about augur and pathogen builds.

# Building

All 32 builds (4 subtypes x 8 segments) can be build by running `snakemake`. For rapid AWS rebuild run as:

    nextstrain build --aws-batch --aws-batch-cpus 16 --aws-batch-memory 28800 . --jobs 16

Please see [nextstrain.org/docs](https://nextstrain.org/docs) for details about augur and pathogen builds.

# A note on clade labeling
H5 viruses are classified into clades, which are currently viewable as a color by on [nextstrain.org](https://nextstrain.org/flu/avian/h5n1/ha?c=h5_label_clade). Because clade annotations are not available in all public databases, we annotate sequences with their most likely clade using a tool developed by Samuel S. Shepard at CDC called [LABEL](https://wonder.cdc.gov/amd/flu/label/). The assigned clade for each H5N1 or H5Nx sequence are available as public tsvs [here](https://github.com/nextstrain/avian-flu/tree/master/clade-labeling).

To update the `clades.tsv` file with clade annotations for new sequences, run: 

`snakemake -s -p Snakefile.clades --cores 1`

To run the builds on without updating the clades file, run: 

`snakemake -p --cores 1`

To string these together and update the `clades.tsv` file for new sequences and then run the builds: 

`snakemake -s -p Snakefile.clades --cores 1 && snakemake -p --cores 1`
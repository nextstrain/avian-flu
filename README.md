# nextstrain.org/flu/avian

This is the Nextstrain build for avian influenza subtypes A/H5N1, A/H7N9, and A/H9N2.
The most up-to-date builds of avian influenza can be found [on nextstrain.org](https://nextstrain.org/flu/avian).
Please see [nextstrain.org/docs](https://nextstrain.org/docs) for details about augur and pathogen builds.

# Building

All 24 builds (3 subtypes x 8 segments) can be build by running `snakemake`. For rapid AWS rebuild run as:

    nextstrain build --aws-batch --aws-batch-cpus 16 --aws-batch-memory 28800 . --jobs 24
Please see [nextstrain.org/docs](https://nextstrain.org/docs) for details about augur and pathogen builds.

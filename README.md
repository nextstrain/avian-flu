# nextstrain.org/flu/avian

This is the Nextstrain build for avian influenza subtypes A/H5N1 and A/H7N9.
The most up-to-date build of H7N9 can be found at [nextstrain.org](https://nextstrain.org/flu/avian/h7n9/ha).
Please see [nextstrain.org/docs](https://nextstrain.org/docs) for details about augur and pathogen builds.

# Building

All 16 builds (2 subtypes x 8 segments) can be build by running `snakemake`. For rapid AWS rebuild run as:

    nextstrain build --aws-batch --aws-batch-cpus 16 --aws-batch-memory 28800 . --jobs 16

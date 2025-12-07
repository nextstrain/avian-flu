This directory represents the workflow which takes curated data (TSV + FASTAs) from our all-influenza curation pipeline and adds in GenoFLU metadata, and then uploads the (new) metadata and (unchanged) sequences to S3. The all-influenza curation pipeline is part of the [ingest workflow in the seasonal-flu repo](https://github.com/nextstrain/seasonal-flu/tree/master/ingest).

See the sister `config.yaml` for the S3 addresses.

**GitHub action**

The `genoflu-gisaid` GitHub action runs this workflow.
The intention is for the seasonal-flu repo to trigger it when newly curated data are available.

**Maual usage**

> You may wish to edit the `config.yaml` to source locally curated data

Working directory: `avian-flu/ingest`

Command: `snakemake --cores 1 --snakefile gisaid/Snakefile -npf`

Add `upload_all` to the end of that rule if you also want to upload files.



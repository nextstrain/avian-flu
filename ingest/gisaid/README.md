This directory represents the workflow which takes curated data (TSV + FASTAs) from our all-influenza curation pipeline (part of the seasonal-flu repo) and adds in GenoFLU metadata, and then uploads the (new) metadata and (unchanged) sequences to S3.

See the sister `config.yaml` for the S3 addresses.

**GitHub action**

The `genoflu-gisaid` GitHub action runs this workflow.
The intention is for the seasonal-flu repo to trigger it when newly curated data are available.

**Maual usage**

Working directory: `avian-flu/ingest`

Command: `snakemake --cores 1 --snakefile gisaid/Snakefile -npf`

Add `upload_all` to the end of that rule if you also want to upload files.



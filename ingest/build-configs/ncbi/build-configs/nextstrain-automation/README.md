# Nextstrain automation

> [!NOTE]
> External users can ignore this directory!
> This build config/customization is tailored for the internal Nextstrain team
> to extend the core ingest workflow for automated workflows.

## Update the config

Update the [config.yaml](config.yaml) for your pathogen:

1. Edit the `s3_dst` param to add the pathogen repository name.
2. Edit the `files_to_upload` param to a mapping of files you need to upload for your pathogen.
The default includes suggested files for uploading curated data and Nextclade outputs.

## Run the workflow

Provide the additional config file to the Snakemake options in order to
include the custom rules from [upload.smk](upload.smk) in the workflow.
Specify the `upload_all` target in order to run the additional upload rules.

The upload rules will require AWS credentials for a user that has permissions
to upload to the Nextstrain data bucket.

The customized workflow can be run from the top level pathogen repo directory with:
```
nextstrain build \
    --env AWS_ACCESS_KEY_ID \
    --env AWS_SECRET_ACCESS_KEY \
    ingest \
        upload_all \
        --configfile build-configs/nextstrain-automation/config.yaml
```

## Automated GitHub Action workflows

Additional instructions on how to use this with the shared `pathogen-repo-build`
GitHub Action workflow to come!

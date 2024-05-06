# nextstrain.org/avian-flu/ingest

This is the ingest pipeline for avian virus sequences.

## Software requirements

Follow the [standard installation instructions](https://docs.nextstrain.org/en/latest/install.html) for Nextstrain's suite of software tools.
This workflow requires the Nextstrain CLI's Docker runtime which includes [fauna](https://github.com/nextstrain/fauna) as a sibling directory of the workflow directory and includes the Python RethinkDB bindings required for downloading from fauna.

## Usage

> NOTE: All command examples assume you are within the `ingest` directory.
> If running commands from the outer `avian-flu` directory, replace the `.` with `ingest`.

### Ingest and upload data from fauna to S3

The ingest pipeline downloads all sequences and metadata from fauna and uploads those data to S3.
Run the complete ingest pipeline and upload results to AWS S3 with the following command.

```sh
nextstrain build \
    --docker \
    --env RETHINK_HOST \
    --env RETHINK_AUTH_KEY \
    --env AWS_ACCESS_KEY_ID \
    --env AWS_SECRET_ACCESS_KEY \
    .
```

Locally, this workflow produces one metadata file, `results/metadata.tsv`, and one sequences file per gene segment like `results/sequences_ha.fasta`.
The workflow compresses and uploads these files to S3 to corresponding paths like `s3://nextstrain-data-private/files/workflows/avian-flu/metadata.tsv.zst` and `s3://nextstrain-data-private/files/workflows/avian-flu/ha/sequences.fasta.zst`.
Each file represents all available subtypes.


### Ingest data without uploading to S3

Run the above command however add specify the rule `all_local`, e.g.

```sh
nextstrain build \
    ... \
    . all_local
```


## Configuration

### Environment Variables

The complete ingest pipeline with AWS S3 uploads uses the following environment variables:

#### Required

- `RETHINK_HOST`
- `RETHINK_AUTH_KEY`
- `AWS_ACCESS_KEY_ID`
- `AWS_SECRET_ACCESS_KEY`

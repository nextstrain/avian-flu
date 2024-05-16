# nextstrain.org/avian-flu/ingest

This is the ingest pipeline for avian virus sequences.

## Software requirements

Follow the [standard installation instructions](https://docs.nextstrain.org/en/latest/install.html) for Nextstrain's suite of software tools.
This workflow requires the Nextstrain CLI's Docker runtime which includes [fauna](https://github.com/nextstrain/fauna) as a sibling directory of the workflow directory and includes the Python RethinkDB bindings required for downloading from fauna.

## Usage

> NOTE: All command examples assume you are within the `ingest` directory.
> If running commands from the outer `avian-flu` directory, replace the `.` with `ingest`.

### Ingest and upload data from fauna to S3

The ingest pipeline supports downloading all sequences and metadata from fauna and uploading those data to S3.

To download all data locally from fauna, run the ingest pipeline with the following command.

```sh
nextstrain build \
    --docker \
    --env RETHINK_HOST \
    --env RETHINK_AUTH_KEY \
    --env AWS_ACCESS_KEY_ID \
    --env AWS_SECRET_ACCESS_KEY \
    .
```

This command produces one metadata file, `results/fauna/metadata.tsv`, and one sequences file per gene segment like `results/fauna/sequences_ha.fasta`.
Each file represents all available subtypes.

Add the `upload_all` target to the command above to run the complete ingest pipeline _and_ upload results to AWS S3.
The workflow compresses and uploads the local files to S3 to corresponding paths like `s3://nextstrain-data-private/files/workflows/avian-flu/metadata.tsv.zst` and `s3://nextstrain-data-private/files/workflows/avian-flu/ha/sequences.fasta.zst`.

```sh
nextstrain build \
    --docker \
    --env RETHINK_HOST \
    --env RETHINK_AUTH_KEY \
    --env AWS_ACCESS_KEY_ID \
    --env AWS_SECRET_ACCESS_KEY \
    . upload_all
```

### Ingest from Andersen lab's avian-influenza repo

Ingest publicly available consensus sequences and metadata from Andersen lab's [avian-influenza repo](https://github.com/andersen-lab/avian-influenza).
Only run this workflow as needed to see the latest available data in the repo.
It does not merge or deduplicate the data with the fauna data used in the default ingest workflow.

```sh
nextstrain build . merge_andersen_segment_metadata
```

The results will be available in `results/andersen-lab/`.

## Configuration

### Environment Variables

The complete ingest pipeline with AWS S3 uploads uses the following environment variables:

#### Required

- `RETHINK_HOST`
- `RETHINK_AUTH_KEY`
- `AWS_ACCESS_KEY_ID`
- `AWS_SECRET_ACCESS_KEY`

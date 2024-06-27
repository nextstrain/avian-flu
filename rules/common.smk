S3_SRC = config.get('s3_src', "s3://nextstrain-data-private/files/workflows/avian-flu")
LOCAL_INGEST = bool(config.get('local_ingest', False))
INGEST_SOURCE = config.get('ingest_source', 'fauna')


def subtypes_by_subtype_wildcard(wildcards):
    db = {
        'h5nx': ['h5n1', 'h5n2', 'h5n3', 'h5n4', 'h5n5', 'h5n6', 'h5n7', 'h5n8', 'h5n9'],
        'h5n1': ['h5n1'],
        'h7n9': ['h7n9'],
        'h9n2': ['h9n2'],
    }
    db['h5n1-cattle-outbreak'] = [*db['h5nx']]
    return(db[wildcards.subtype])

if LOCAL_INGEST:
    rule copy_sequences_from_ingest:
        output:
            sequences = "data/{segment}/sequences.fasta",
        params:
            sequences = lambda w: f"ingest/{INGEST_SOURCE}/results/sequences_{w.segment}.fasta"
        shell:
            """
            cp {params.sequences} {output.sequences}
            """

    rule copy_metadata_from_ingest:
        output:
            metadata = "data/metadata.tsv",
        params:
            metadata = f"ingest/{INGEST_SOURCE}/results/metadata.tsv",
        shell:
            """
            cp {params.metadata} {output.metadata}
            """

else:
    rule download_sequences:
        output:
            sequences = "data/{segment}/sequences.fasta",
        params:
            s3_src=S3_SRC,
        shell:
            """
            aws s3 cp {params.s3_src:q}/{wildcards.segment}/sequences.fasta.zst - | zstd -d > {output.sequences}
            """

    rule download_metadata:
        output:
            metadata = "data/metadata.tsv",
        params:
            s3_src=S3_SRC,
        shell:
            """
            aws s3 cp {params.s3_src:q}/metadata.tsv.zst - | zstd -d > {output.metadata}
            """

rule filter_sequences_by_subtype:
    input:
        sequences = "data/{segment}/sequences.fasta",
        metadata = "data/metadata.tsv",
    output:
        sequences = "data/sequences_{subtype}_{segment}.fasta",
    params:
        subtypes=subtypes_by_subtype_wildcard,
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --query "subtype in {params.subtypes!r}" \
            --output-sequences {output.sequences}
        """

rule filter_metadata_by_subtype:
    input:
        metadata = "data/metadata.tsv",
    output:
        metadata = "data/metadata_{subtype}.tsv",
    params:
        subtypes=subtypes_by_subtype_wildcard,
    shell:
        """
        augur filter \
            --metadata {input.metadata} \
            --query "subtype in {params.subtypes!r}" \
            --output-metadata {output.metadata}
        """

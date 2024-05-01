path_to_fauna = '../fauna'
config["s3_dst"] = "s3://nextstrain-data-private/files/workflows/avian-flu"
config["segments"] = ["pb2", "pb1", "pa", "ha","np", "na", "mp", "ns"]

rule all:
    input:
        sequences=expand("upload/s3/sequences_{segment}.done", segment=config["segments"]),
        metadata="upload/s3/metadata.done",

rule download_segment:
    output:
        sequences = "upload/data/{segment}.fasta",
    params:
        fasta_fields = "strain virus accession collection_date region country division location host domestic_status subtype originating_lab submitting_lab authors PMID gisaid_clade h5_clade",
    benchmark:
        "benchmarks/download_segment_{segment}.txt"
    shell:
        """
        python3 {path_to_fauna}/vdb/download.py \
            --database vdb \
            --virus avian_flu \
            --fasta_fields {params.fasta_fields} \
            --select  locus:{wildcards.segment} \
            --path upload/data \
            --fstem {wildcards.segment}
        """

rule parse_segment:
    input:
        sequences = "upload/data/{segment}.fasta",
    output:
        sequences = "upload/results/sequences_{segment}.fasta",
        metadata = "upload/results/metadata_{segment}.tsv",
    params:
        fasta_fields =  "strain virus isolate_id date region country division location host domestic_status subtype originating_lab submitting_lab authors PMID gisaid_clade h5_clade",
        prettify_fields = "region country division location host originating_lab submitting_lab authors PMID"
    shell:
        """
        augur parse \
            --sequences {input.sequences} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --fields {params.fasta_fields} \
            --prettify-fields {params.prettify_fields}
        """

rule merge_segment_metadata:
    """
    For each subtype's HA metadata file add a column "n_segments" which reports
    how many segments have sequence data (no QC performed). This will force the
    download & parsing of all segments for a given subtype. Note that this does
    not currently consider the prescribed min lengths (see min_length function)
    for each segment, but that would be a nice improvement.
    """
    input:
        segments = expand("upload/results/metadata_{segment}.tsv", segment=config["segments"]),
        metadata = "upload/results/metadata_ha.tsv",
    output:
        metadata = "upload/results/metadata.tsv",
    shell:
        """
        python scripts/add_segment_counts.py \
            --segments {input.segments} \
            --metadata {input.metadata} \
            --output {output.metadata}
        """

rule upload_sequences:
    input:
        sequences="upload/results/sequences_{segment}.fasta",
    output:
        flag=touch("upload/s3/sequences_{segment}.done"),
    params:
        s3_dst=config["s3_dst"],
    shell:
        """
        zstd -c {input.sequences:q} \
            | aws s3 cp \
                  - \
                  {params.s3_dst:q}/{wildcards.segment}/sequences.fasta.zst
        """

rule upload_metadata:
    input:
        metadata="upload/results/metadata.tsv",
    output:
        flag=touch("upload/s3/metadata.done"),
    params:
        s3_dst=config["s3_dst"],
    shell:
        """
        zstd -c {input.metadata:q} \
            | aws s3 cp \
                  - \
                  {params.s3_dst:q}/metadata.tsv.zst
        """

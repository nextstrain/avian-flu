path_to_fauna = '../fauna'
config["s3_dst"] = "s3://nextstrain-data-private/files/workflows/avian-flu"
config["segments"] = ["pb2", "pb1", "pa", "ha","np", "na", "mp", "ns"]

rule all:
    input:
        sequences=expand("upload/s3/sequences_{segment}.done", segment=config["segments"]),
        metadata=expand("upload/s3/metadata_{segment}.done", segment=config["segments"]),

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
        metadata="upload/results/metadata_{segment}.tsv",
    output:
        flag=touch("upload/s3/metadata_{segment}.done"),
    params:
        s3_dst=config["s3_dst"],
    shell:
        """
        zstd -c {input.metadata:q} \
            | aws s3 cp \
                  - \
                  {params.s3_dst:q}/{wildcards.segment}/metadata.tsv.zst
        """

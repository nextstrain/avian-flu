"""
This part of the workflow handles uploading files to AWS S3.
"""


rule upload_sequences:
    input:
        sequences="{data_source}/results/sequences_{segment}.fasta",
    output:
        flag=touch("{data_source}/s3/sequences_{segment}.done"),
    params:
        s3_dst=lambda wildcards: config["s3_dst"][wildcards.data_source],
    shell:
        """
        zstd -c {input.sequences:q} \
            | aws s3 cp \
                  - \
                  {params.s3_dst:q}/{wildcards.segment}/sequences.fasta.zst
        """


rule upload_metadata:
    input:
        metadata="{data_source}/results/metadata.tsv",
    output:
        flag=touch("{data_source}/s3/metadata.done"),
    params:
        s3_dst=lambda wildcards: config["s3_dst"][wildcards.data_source],
    shell:
        """
        zstd -c {input.metadata:q} \
            | aws s3 cp \
                  - \
                  {params.s3_dst:q}/metadata.tsv.zst
        """

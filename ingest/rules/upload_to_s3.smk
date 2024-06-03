"""
This part of the workflow handles uploading files to AWS S3.
"""


rule upload_sequences:
    input:
        sequences="{data_source}/results/sequences_{segment}.fasta",
    output:
        flag="{data_source}/s3/sequences_{segment}.done",
    params:
        s3_dst=lambda wildcards: config["s3_dst"][wildcards.data_source],
        cloudfront_domain=config.get("cloudfront_domain", ""),
    shell:
        """
        ./vendored/upload-to-s3 \
            --quiet \
            {input.sequences:q} \
            {params.s3_dst:q}/{wildcards.segment}/sequences.fasta.zst \
            {params.cloudfront_domain} 2>&1 | tee {output.flag}
        """


rule upload_metadata:
    input:
        metadata="{data_source}/results/metadata.tsv",
    output:
        flag="{data_source}/s3/metadata.done",
    params:
        s3_dst=lambda wildcards: config["s3_dst"][wildcards.data_source],
        cloudfront_domain=config.get("cloudfront_domain", ""),
    shell:
        """
        ./vendored/upload-to-s3 \
            --quiet \
            {input.metadata:q} \
            {params.s3_dst:q}/metadata.tsv.zst \
            {params.cloudfront_domain} 2>&1 | tee {output.flag}
        """

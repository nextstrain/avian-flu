"""
This part of the workflow handles uploading files to AWS S3.

Files to upload must be defined in the `files_to_upload` config param, where
the keys are the remote files and the values are the local filepaths
relative to the ingest directory.

Produces a single file for each uploaded file:
    "results/upload/{remote_file}.upload"

The rule `upload_all` can be used as a target to upload all files.
"""
import os

slack_envvars_defined = "SLACK_CHANNELS" in os.environ and "SLACK_TOKEN" in os.environ
send_notifications = (
    config.get("send_slack_notifications", False) and slack_envvars_defined
)


rule upload_to_s3:
    input:
        file_to_upload=lambda wildcards: config["files_to_upload"][wildcards.remote_file],
    output:
        "results/upload/{remote_file}.upload",
    params:
        quiet="" if send_notifications else "--quiet",
        s3_dst=config["s3_dst"],
        cloudfront_domain=config["cloudfront_domain"],
    shell:
        """
        ./vendored/upload-to-s3 \
            {params.quiet} \
            {input.file_to_upload:q} \
            {params.s3_dst:q}/{wildcards.remote_file:q} \
            {params.cloudfront_domain} 2>&1 | tee {output}
        """


rule upload_all:
    input:
        uploads=[
            f"results/upload/{remote_file}.upload"
            for remote_file in config["files_to_upload"].keys()
        ],
    output:
        touch("results/upload_all.done")

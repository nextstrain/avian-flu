DEPLOY_URL = config.get('deploy_url', "s3://nextstrain-data")


rule deploy_all:
    """
    Upload all builds to AWS S3
    Depends on indendent Snakemake workflow's defined `all` rule
    """
    input: rules.all.input
    params:
        s3_dst = DEPLOY_URL
    shell:
        """
        nextstrain remote upload {params.s3_dst:q} {input}
        """

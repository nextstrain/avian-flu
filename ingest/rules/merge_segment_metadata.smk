"""
This part of the workflow handles how we merge the metadata for each segment
into a central metadata file.
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
        segments = expand("{{data_source}}/data/metadata_{segment}.tsv", segment=config["segments"]),
        metadata = "{data_source}/data/metadata_ha.tsv",
    output:
        metadata = "{data_source}/results/metadata.tsv",
    shell:
        """
        python scripts/add_segment_counts.py \
            --segments {input.segments} \
            --metadata {input.metadata} \
            --output {output.metadata}
        """

"""
This part of the workflow handles the curation of data from NCBI

REQUIRED INPUTS:

    ndjson      = ncbi/data/ncbi.ndjson

OUTPUTS:

    metadata    = ncbi/results/metadata.tsv
    seuqences   = ncbi/results/sequences_{segment}.fasta

"""


rule fetch_general_geolocation_rules:
    output:
        general_geolocation_rules="ncbi/data/general-geolocation-rules.tsv",
    params:
        geolocation_rules_url=config["curate"]["geolocation_rules_url"],
    shell:
        """
        curl {params.geolocation_rules_url} > {output.general_geolocation_rules}
        """


rule concat_geolocation_rules:
    input:
        general_geolocation_rules="ncbi/data/general-geolocation-rules.tsv",
        local_geolocation_rules=config["curate"]["local_geolocation_rules"],
    output:
        all_geolocation_rules="ncbi/data/all-geolocation-rules.tsv",
    shell:
        """
        cat {input.general_geolocation_rules} {input.local_geolocation_rules} >> {output.all_geolocation_rules}
        """


def format_field_map(field_map: dict[str, str]) -> str:
    """
    Format dict to `"key1"="value1" "key2"="value2"...` for use in shell commands.
    """
    return " ".join([f'"{key}"="{value}"' for key, value in field_map.items()])


rule curate:
    input:
        sequences_ndjson="ncbi/data/ncbi.ndjson",
        all_geolocation_rules="ncbi/data/all-geolocation-rules.tsv",
        annotations=config["curate"]["annotations"],
    output:
        curated_ndjson="ncbi/data/curated.ndjson",
    log:
        "ncbi/logs/curate.txt",
    benchmark:
        "ncbi/benchmarks/curate.txt"
    params:
        field_map=format_field_map(config["curate"]["field_map"]),
        strain_regex=config["curate"]["strain_regex"],
        strain_backup_fields=config["curate"]["strain_backup_fields"],
        date_fields=config["curate"]["date_fields"],
        expected_date_formats=config["curate"]["expected_date_formats"],
        articles=config["curate"]["titlecase"]["articles"],
        abbreviations=config["curate"]["titlecase"]["abbreviations"],
        titlecase_fields=config["curate"]["titlecase"]["fields"],
        authors_field=config["curate"]["authors_field"],
        authors_default_value=config["curate"]["authors_default_value"],
        abbr_authors_field=config["curate"]["abbr_authors_field"],
        segments=format_field_map(config["ncbi_segments"]),
        host_map=config["curate"]["host_map"],
        annotations_id=config["curate"]["annotations_id"],
    shell:
        """
        (cat {input.sequences_ndjson} \
            | ./vendored/transform-field-names \
                --field-map {params.field_map} \
            | augur curate normalize-strings \
            | ./vendored/transform-strain-names \
                --strain-regex {params.strain_regex} \
                --backup-fields {params.strain_backup_fields} \
            | ./build-configs/ncbi/bin/parse-metadata-from-strain \
            | augur curate format-dates \
                --date-fields {params.date_fields} \
                --expected-date-formats {params.expected_date_formats} \
            | ./vendored/transform-genbank-location \
            | augur curate titlecase \
                --titlecase-fields {params.titlecase_fields} \
                --articles {params.articles} \
                --abbreviations {params.abbreviations} \
            | ./vendored/transform-authors \
                --authors-field {params.authors_field} \
                --default-value {params.authors_default_value} \
                --abbr-authors-field {params.abbr_authors_field} \
            | ./vendored/apply-geolocation-rules \
                --geolocation-rules {input.all_geolocation_rules} \
            | ./build-configs/ncbi/bin/transform-segment \
                --segments {params.segments} \
            | ./build-configs/ncbi/bin/transform-host \
                --host-map {params.host_map} \
            | ./build-configs/ncbi/bin/transform-to-match-fauna \
            | ./vendored/merge-user-metadata \
                --annotations {input.annotations} \
                --id-field {params.annotations_id} ) 2>> {log} > {output.curated_ndjson}
        """


rule split_curated_ndjson_by_segment:
    """
    Split out the full curate NDJSON by segment, then we can deduplicate
    records by strain name within each segment
    """
    input:
        curated_ndjson="ncbi/data/curated.ndjson",
    output:
        metadata="ncbi/data/all_metadata_{segment}.tsv",
        sequences="ncbi/results/sequences_{segment}.fasta",
    params:
        id_field=config["curate"]["output_id_field"],
        sequence_field=config["curate"]["output_sequence_field"],
    log:
        "ncbi/logs/{segment}/split_curated_ndjson_by_segment.txt",
    benchmark:
        "ncbi/benchmarks/{segment}/split_curated_ndjson_by_segment.txt"
    shell:
        """
        (cat {input.curated_ndjson} \
            | ./build-configs/ncbi/bin/filter-ndjson-by-segment \
                --segment {wildcards.segment} \
            | ./build-configs/ncbi/bin/dedup-by-strain \
            | augur curate passthru \
                --output-metadata {output.metadata} \
                --output-fasta {output.sequences} \
                --output-id-field {params.id_field} \
                --output-seq-field {params.sequence_field} ) 2>> {log}
        """


rule subset_metadata:
    input:
        metadata="ncbi/data/all_metadata_{segment}.tsv",
    output:
        subset_metadata="ncbi/data/metadata_{segment}.tsv",
    params:
        metadata_fields=",".join(config["curate"]["metadata_columns"]),
    shell:
        """
        tsv-select -H -f {params.metadata_fields} \
            {input.metadata} > {output.subset_metadata}
        """


rule merge_ncbi_segment_metadata:
    """
    Add a column "n_segments" which reports how many segments
    have sequence data (no QC performed).
    """
    input:
        segments = expand("ncbi/data/metadata_{segment}.tsv", segment=config["ncbi_segments"]),
        metadata = "ncbi/data/metadata_ha.tsv",
    output:
        metadata = "ncbi/results/metadata.tsv",
    shell:
        """
        python scripts/add_segment_counts.py \
            --segments {input.segments} \
            --metadata {input.metadata} \
            --output {output.metadata}
        """

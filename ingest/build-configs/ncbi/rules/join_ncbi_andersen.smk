"""
This part of the workflow handles the joining and deduplication of NCBI and
Andersen lab data.
"""


rule select_missing_metadata:
    """
    Uses tsv-join --exclude flag to exclude matching records so we are
    left with any missing metadata from the Andersen lab.
    """
    input:
        ncbi_metadata = "ncbi/results/metadata.tsv",
        andersen_metadata = "andersen-lab/results/metadata.tsv",
    output:
        missing_metadata = "joined-ncbi/data/missing_metadata.tsv",
    params:
        match_field = config["join_ncbi_andersen"]["match_field"],
    shell:
        """
        tsv-join -H \
            --exclude \
            --filter-file {input.ncbi_metadata} \
            --key-fields {params.match_field} \
            {input.andersen_metadata} > {output.missing_metadata}
        """


rule append_missing_metadata_to_ncbi:
    input:
        ncbi_metadata = "ncbi/results/metadata.tsv",
        missing_metadata = "joined-ncbi/data/missing_metadata.tsv",
    output:
        joined_metadata = "joined-ncbi/data/all_metadata.tsv",
    params:
        source_column_name = config["join_ncbi_andersen"]["source_column_name"],
        ncbi_source = config["join_ncbi_andersen"]["ncbi_source"],
        andersen_source = config["join_ncbi_andersen"]["andersen_source"],
        dedup_column = config["join_ncbi_andersen"]["dedup_column"],
    shell:
        """
        tsv-append \
            --source-header {params.source_column_name} \
            --file {params.ncbi_source}={input.ncbi_metadata} \
            --file {params.andersen_source}={input.missing_metadata} \
            | tsv-uniq -H -f {params.dedup_column} \
            > {output.joined_metadata}
        """


rule dedup_metadata_by_sample_id:
    input:
        all_joined_metadata="joined-ncbi/data/all_metadata.tsv",
    output:
        joined_metadata="joined-ncbi/results/metadata.tsv",
    params:
        id_field=config["curate"]["output_id_field"],
    log:
        "logs/joined-ncbi/dedup_metadata_by_sample_id.txt",
    shell:
        r"""
        (augur curate passthru \
            --id-column {params.id_field:q} \
            --metadata {input.all_joined_metadata:q} \
            | ./build-configs/ncbi/bin/dedup-by-sample-id \
            | augur curate passthru \
                --output-metadata {output.joined_metadata}) 2>> {log}
        """


rule select_final_strain_names:
    input:
        joined_metadata="joined-ncbi/results/metadata.tsv",
    output:
        strain_names = "joined-ncbi/data/final_strain_names.txt",
    params:
        sequence_id_column = config["curate"]["output_id_field"],
    shell:
        r"""
        tsv-select -H -f {params.sequence_id_column} \
            {input.joined_metadata} \
            > {output.strain_names}
        """


rule select_final_sequences:
    input:
        strain_names = "joined-ncbi/data/final_strain_names.txt",
        ncbi_sequences = "ncbi/results/sequences_{segment}.fasta",
        andersen_sequences = "andersen-lab/results/sequences_{segment}.fasta",
    output:
        joined_sequences = "joined-ncbi/results/sequences_{segment}.fasta",
    shell:
        r"""
        cat {input.ncbi_sequences} {input.andersen_sequences} \
            | seqkit grep -f {input.strain_names} \
                > {output.joined_sequences}
        """

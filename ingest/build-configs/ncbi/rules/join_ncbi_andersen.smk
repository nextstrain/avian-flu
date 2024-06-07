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


rule select_missing_strain_names:
    input:
        missing_metadata = "joined-ncbi/data/missing_metadata.tsv",
    output:
        missing_sequence_ids = "joined-ncbi/data/missing_sequence_ids.txt",
    params:
        sequence_id_column = config["curate"]["output_id_field"],
    shell:
        """
        tsv-select -H -f {params.sequence_id_column} \
            {input.missing_metadata} \
            > {output.missing_sequence_ids}
        """


rule select_missing_sequences:
    input:
        missing_sequence_ids = "joined-ncbi/data/missing_sequence_ids.txt",
        andersen_sequences = "andersen-lab/results/sequences_{segment}.fasta",
    output:
        missing_sequences = "joined-ncbi/data/missing_sequences_{segment}.fasta",
    shell:
        """
        seqkit grep -f {input.missing_sequence_ids} \
            {input.andersen_sequences} \
            > {output.missing_sequences}
        """


rule append_missing_metadata_to_ncbi:
    input:
        ncbi_metadata = "ncbi/results/metadata.tsv",
        missing_metadata = "joined-ncbi/data/missing_metadata.tsv",
    output:
        joined_metadata = "joined-ncbi/results/metadata.tsv",
    params:
        source_column_name = config["join_ncbi_andersen"]["source_column_name"],
        ncbi_source = config["join_ncbi_andersen"]["ncbi_source"],
        andersen_source = config["join_ncbi_andersen"]["andersen_source"],
    shell:
        """
        tsv-append \
            --source-header {params.source_column_name} \
            --file {params.ncbi_source}={input.ncbi_metadata} \
            --file {params.andersen_source}={input.missing_metadata} \
            > {output.joined_metadata}
        """


rule append_missing_sequences_to_ncbi:
    input:
        ncbi_sequences = "ncbi/results/sequences_{segment}.fasta",
        missing_sequences = "joined-ncbi/data/missing_sequences_{segment}.fasta",
    output:
        joined_sequences = "joined-ncbi/results/sequences_{segment}.fasta",
    shell:
        """
        cat {input.ncbi_sequences} {input.missing_sequences} > {output.joined_sequences}
        """

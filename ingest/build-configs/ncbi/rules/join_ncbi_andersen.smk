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

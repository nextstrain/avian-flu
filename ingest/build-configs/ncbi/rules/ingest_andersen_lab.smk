"""
This part of the workflow handles ingest of metadata and consensus sequences
from the Andersen Lab's avian-influenza repo
<https://github.com/andersen-lab/avian-influenza>
"""

rule fetch_andersen_lab_repo:
    output:
        andersen_lab_repo = temp("andersen-lab/data/avian-influenza.tar.gz")
    shell:
        """
        curl -fsSL \
            -H "Accept: application/vnd.github+json" \
            -H "X-GitHub-Api-Version: 2022-11-28" \
            https://api.github.com/repos/andersen-lab/avian-influenza/tarball \
            > {output.andersen_lab_repo}
        """

rule extract_metadata:
    input:
        andersen_lab_repo = "andersen-lab/data/avian-influenza.tar.gz"
    output:
        metadata = "andersen-lab/data/PRJNA1102327_metadata.csv"
    params:
        metadata_file_path = "metadata/SraRunTable_PRJNA1102327_automated.csv",
    shell:
        """
        tar xz -O --file={input.andersen_lab_repo} \
            --wildcards \
            "*/{params.metadata_file_path:q}" \
            > {output.metadata}
        """

rule extract_consensus_sequences:
    input:
        andersen_lab_repo = "andersen-lab/data/avian-influenza.tar.gz"
    output:
        fasta = directory("andersen-lab/data/fasta"),
        output_flag = touch("andersen-lab/data/extract_consensus_sequences.done")
    params:
        output_dir = lambda wildcards, output: Path(output.fasta).parent
    shell:
        """
        tar xz --file={input.andersen_lab_repo} \
            --strip-components=1 \
            -C {params.output_dir} \
            --wildcards \
            "*/fasta"
        """

rule rename_and_concatenate_segment_fastas:
    """
    Truncate the FASTA headers to just the SRA run accessions
    and concatenate FASTAs of the same segment
    """
    input:
        extract_consensus_sequences_flag = "andersen-lab/data/extract_consensus_sequences.done"
    output:
        fasta = "andersen-lab/data/{segment}.fasta"
    params:
        segment = lambda wildcards: wildcards.segment.upper()
    shell:
        """
        for fasta in andersen-lab/data/fasta/*_{params.segment}_cns.fa; do
            seqkit replace \
                -p "Consensus_(SRR[0-9]+)_.*" \
                -r '$1' \
                "$fasta" \
                >> {output.fasta}
        done
        """

rule curate_metadata:
    input:
        metadata = "andersen-lab/data/PRJNA1102327_metadata.csv",
        geolocation_rules = "defaults/geolocation_rules.tsv"
    output:
        metadata = "andersen-lab/data/metadata.tsv"
    log:
        "andersen-lab/logs/curate_metadata.txt",
    shell:
        """
        augur curate normalize-strings \
            --metadata {input.metadata} \
            | ./build-configs/ncbi/bin/curate-andersen-lab-data \
            | ./vendored/apply-geolocation-rules \
                --geolocation-rules {input.geolocation_rules} \
            | augur curate passthru \
                --output-metadata {output.metadata} 2>> {log}
        """

rule match_metadata_and_segment_fasta:
    """
    Matches the full metadata with the corresponding segment sequence FASTAs
    and outputs the matching metadata TSV and sequence FASTAs per segment.
    """
    input:
        metadata = "andersen-lab/data/metadata.tsv",
        fasta = "andersen-lab/data/{segment}.fasta"
    output:
        metadata = "andersen-lab/data/matched_metadata_{segment}.tsv",
        fasta = "andersen-lab/results/sequences_{segment}.fasta",
    params:
        input_id_field="isolate_id",
        sequence_field="sequence",
        output_id_field=config["curate"]["output_id_field"],
    log:
        "andersen-lab/logs/match_segment_metadata_and_fasta/{segment}.txt",
    shell:
        """
        augur curate passthru \
            --metadata {input.metadata} \
            --fasta {input.fasta} \
            --seq-id-column {params.input_id_field} \
            --seq-field {params.sequence_field} \
            --unmatched-reporting warn \
            --duplicate-reporting warn \
            --output-metadata {output.metadata} \
            --output-fasta {output.fasta} \
            --output-id-field {params.output_id_field} \
            --output-seq-field {params.sequence_field} \
            2> {log}
        """


rule reorder_metadata_columns:
    """
    Using tsv-select to reorder the columns of the Andersen lab metadata to
    exactly match the NCBI metadata columns. Ensures that we don't accidently
    append the wrong columns in joining steps.

    tsv-select will exit with error if the column does not exist.
    """
    input:
        metadata = "andersen-lab/data/matched_metadata_{segment}.tsv",
    output:
        reordered_metadata = "andersen-lab/data/metadata_{segment}.tsv"
    params:
        metadata_fields=",".join(config["curate"]["metadata_columns"]),
    shell:
        """
        tsv-select -H -f {params.metadata_fields} \
            {input.metadata} > {output.reordered_metadata}
        """

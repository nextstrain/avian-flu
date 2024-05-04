"""
This part of the workflow handles ingest of metadata and consensus sequences
from the Andersen Lab's avian-influenza repo
<https://github.com/andersen-lab/avian-influenza>
"""

rule fetch_andersen_lab_repo:
    output:
        andersen_lab_repo = temp("data/andersen-lab-avian-influenza.tar.gz")
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
        andersen_lab_repo = "data/andersen-lab-avian-influenza.tar.gz"
    output:
        metadata = "data/andersen-lab/PRJNA1102327_metadata.csv"
    shell:
        """
        tar xz --file={input.andersen_lab_repo} \
            --strip-components=2 \
            -C data/andersen-lab \
            --wildcards \
            "*/metadata/PRJNA1102327_metadata.csv"
        """

rule extract_consensus_sequences:
    input:
        andersen_lab_repo = "data/andersen-lab-avian-influenza.tar.gz"
    output:
        fasta = directory("data/andersen-lab/fasta"),
        output_flag = touch("data/andersen-lab/extract_consensus_sequences.done")
    shell:
        """
        tar xz --file={input.andersen_lab_repo} \
            --strip-components=1 \
            -C data/andersen-lab \
            --wildcards \
            "*/fasta"
        """

rule rename_and_concatenate_segment_fastas:
    """
    Truncate the FASTA headers to just the SRA run accessions
    and concatenate FASTAs of the same segment
    """
    input:
        extract_consensus_sequences_flag = "data/andersen-lab/extract_consensus_sequences.done"
    output:
        fasta = "data/andersen-lab/{segment}.fasta"
    params:
        segment = lambda wildcards: wildcards.segment.upper()
    shell:
        """
        for fasta in data/andersen-lab/fasta/*_{params.segment}_cns.fa; do
            seqkit replace \
                -p "Consensus_(SRR[0-9]+)_.*" \
                -r '$1' \
                "$fasta" \
                >> {output.fasta}
        done
        """

rule curate_metadata:
    input:
        metadata = "data/andersen-lab/PRJNA1102327_metadata.csv",
        geolocation_rules = "defaults/geolocation_rules.tsv"
    output:
        metadata = "data/andersen-lab/metadata.tsv"
    log:
        "logs/curate_metadata.txt",
    shell:
        """
        augur curate normalize-strings \
            --metadata {input.metadata} \
            | python3 ./scripts/curate_andersen_lab_data.py \
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
        metadata = "data/andersen-lab/metadata.tsv",
        fasta = "data/andersen-lab/{segment}.fasta"
    output:
        metadata = "results/andersen-lab/metadata_{segment}.tsv",
        fasta = "results/andersen-lab/sequences_{segment}.fasta"
    log:
        "logs/match_segment_metadata_and_fasta/{segment}.txt",
    shell:
        """
        augur curate passthru \
            --metadata {input.metadata} \
            --fasta {input.fasta} \
            --seq-id-column isolate_id \
            --seq-field sequence \
            --unmatched-reporting warn \
            --duplicate-reporting warn \
            --output-metadata {output.metadata} \
            --output-fasta {output.fasta} \
            --output-id-field strain \
            --output-seq-field sequence \
            2> {log}
        """

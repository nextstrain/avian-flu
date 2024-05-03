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
        fasta = directory("data/andersen-lab/fasta")
    shell:
        """
        tar xz --file={input.andersen_lab_repo} \
            --strip-components=1 \
            -C data/andersen-lab \
            --wildcards \
            "*/fasta"
        """

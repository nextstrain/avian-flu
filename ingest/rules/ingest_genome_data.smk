
## This ruleset is in flux and will change often
## Currently it merges hardcoded (committed) source-data with fauna-derived data
## No checking is done for duplicate strains

rule merge_genome_metadata:
    input:
        fauna = "results/fauna/metadata.tsv"
    params:
        source_data = "source-data/metadata.tsv"
    output:
        metadata = "results/genome/metadata.tsv"
    shell:
        """
        diff <(head -n 1 {params.source_data}) <(head -n 1 {input.fauna}) &&
            cp {params.source_data} {output.metadata} && \
            tail -n +2 {input.fauna} >> {output.metadata}
        """

rule merge_genome_sequences:
    input:
        fauna = "results/fauna/sequences_{segment}.fasta"
    params:
        source_data = "source-data/sequences_{segment}.fasta"
    output:
        metadata = "results/genome/sequences_{segment}.fasta"
    shell:
        """
        cat {params.source_data} {input.fauna} > {output.metadata}
        """


"""
This part of the workflow handles assigning lineages and genotype with GenoFLU.

We are using a vendored version of <https://github.com/moncla-lab/GenoFLU-multi>,
which is built on top of USDA's GenoFLU <https://github.com/USDA-VS/GenoFLU>.
"""

rule run_genoflu:
    input:
        sequences=expand("{{data_source}}/results/sequences_{segment}.fasta", segment=config["segments"]),
    output:
        # Hardcoded output file from GenoFLU is {fasta_dir}/results/results.tsv
        genoflu="{data_source}/results/results/results.tsv"
    params:
        fasta_dir="{data_source}/results/",
    log:
        "{data_source}/logs/run_genoflu.txt",
    threads: 4
    shell:
        """
        python ./vendored-GenoFLU-multi/bin/genoflu-multi.py \
            -f {params.fasta_dir:q} \
            -n {threads} > {log}
        """


rule subset_genoflu:
    input:
        genoflu="{data_source}/results/results/results.tsv",
    output:
        genotypes="{data_source}/results/genoflu_genotypes.tsv",
    shell:
        """
        csvtk cut -t \
            -f Strain,Genotype \
            {input.genoflu} \
            | csvtk rename -t \
                -f Strain,Genotype \
                -n strain,genoflu_genotype > {output.genotypes}
        """


rule merge_metadata_with_genotypes:
    input:
        metadata="{data_source}/results/metadata.tsv",
        genotypes="{data_source}/results/genoflu_genotypes.tsv",
    output:
        final_metadata="{data_source}/results/final_metadata.tsv",
    shell:
        """
        augur merge \
            --metadata metadata={input.metadata} genotypes={input.genotypes} \
            --metadata-id-columns strain \
            --no-source-columns \
            --output-metadata {output.final_metadata}
        """

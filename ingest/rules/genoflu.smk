"""
This part of the workflow handles assigning lineages and genotype with GenoFLU.

We are using a vendored version of <https://github.com/moncla-lab/GenoFLU-multi>,
which is built on top of USDA's GenoFLU <https://github.com/USDA-VS/GenoFLU>.
"""

rule provision_genoflu_sequences:
    """
    GenoFLU will consume all the FASTA files in the provided directory, so we set up a
    new directory with FASTA files we want to call. (We use the final ("results") sequences
    here because the sequences themselves aren't modified by GenoFLU.)

    The current implementation is a simple file copy however we may wish to use `augur filter`
    in the future to restrict the samples we process.
    """
    input:
        sequences = "{data_source}/results/sequences_{segment}.fasta",
    output:
        sequences = temp("{data_source}/data/genoflu/sequences_{segment}.fasta"),
    threads: 1
    shell:
        """
        cp {input.sequences} {output.sequences}
        """

rule run_genoflu:
    input:
        sequences=expand("{{data_source}}/data/genoflu/sequences_{segment}.fasta", segment=config["segments"]),
    output:
        # Hardcoded output file from GenoFLU is {fasta_dir}/results/results.tsv
        genoflu="{data_source}/data/genoflu/results/results.tsv"
    params:
        fasta_dir="{data_source}/data/genoflu/",
    log:
        "{data_source}/logs/run_genoflu.txt",
    threads: 12
    shell:
        """
        python ./vendored-GenoFLU-multi/bin/genoflu-multi.py \
            -f {params.fasta_dir:q} \
            -n {threads} > {log}
        """


rule parse_genoflu:
    """
    Parses the genoflu TSV to produce a TSV with 10 columns:
    * strain - ID used for matching
    * genoflu - the "genotype" or "constellation"
    * genoflu_<SEGMENT> - the individual segment genoflu calls
    """
    input:
        genoflu="{data_source}/data/genoflu/results/results.tsv"
    output:
        genotypes="{data_source}/data/genoflu/genoflu_genotypes.tsv",
    shell:
        r"""
        cat {input.genoflu} | \
            csvtk cut -t -F -f Strain,Genotype,'Genotype List Used*' | \
            csvtk rename -t -F -f Strain,Genotype,'Genotype List Used*' -n strain,genoflu,details | \
            python scripts/parse_genoflu.py \
            > {output.genotypes}
        """


rule merge_metadata_with_genotypes:
    input:
        metadata = "{data_source}/data/metadata_combined.tsv",
        genotypes = "{data_source}/data/genoflu/genoflu_genotypes.tsv",
    output:
        metadata = "{data_source}/results/metadata.tsv",
    wildcard_constraints:
        # Note - we need a 'placeholder' (or somesuch) to prevent an empty string constraint which
        data_source = '|'.join([*[k for k,v in config.get('genoflu', {}).items() if v], 'placeholder'])
    shell:
        """
        augur merge \
            --metadata metadata={input.metadata} genotypes={input.genotypes} \
            --metadata-id-columns strain \
            --no-source-columns \
            --output-metadata {output.metadata}
        """

rule skip_genoflu:
    input:
        metadata = "{data_source}/data/metadata_combined.tsv",
    output:
        metadata = "{data_source}/results/metadata.tsv",
    shell:
        """
        cp {input.metadata} {output.metadata}
        """

# This ruleorder effectively determines whether GenoFLU will be run or not, and in turn this is
# controlled by the `wildcard_constraints` on `merge_metadata_with_genotypes` which are config-
# controlled.
ruleorder: merge_metadata_with_genotypes > skip_genoflu
"""
This part of the workflow handles assigning lineages and genotype with GenoFLU.

We are using a vendored version of <https://github.com/moncla-lab/GenoFLU-multi>,
which is built on top of USDA's GenoFLU <https://github.com/USDA-VS/GenoFLU>.
"""

def genoflu_filter_args(wildcards):
    # NOTE: it's crucial to get the quoting right here, of the following three strings on the command line
    # "gisaid_clade=='2.3.4.4b'" gisaid_clade=='2.3.4.4b' and 'gisaid_clade==2.3.4.4b'
    # only the first works.
    # NOTE 2: This filtering may not be correct - see <https://github.com/nextstrain/avian-flu/pull/127#issuecomment-2669102995>
    if wildcards.data_source=='fauna':
        return "--query \"gisaid_clade=='2.3.4.4b'\""
    return ""


rule provision_genoflu_sequences:
    """
    GenoFLU will consume all the FASTA files in the provided directory, so we set up a
    new directory with (filtered) FASTA files we want to call. Note that we use the
    final/results sequences as _inputs_ here, because GenoFLU isn't going to modify
    those in any way, and as such they are marked as temporary.
    """
    input:
        sequences = "{data_source}/results/sequences_{segment}.fasta",
        metadata = "{data_source}/data/metadata_combined.tsv",
    output:
        sequences = temp("{data_source}/data/genoflu/sequences_{segment}.fasta"),
    threads: 1
    params:
        query = genoflu_filter_args
    shell:
        """
        augur filter \
            {params.query} \
            --metadata {input.metadata} --sequences {input.sequences} \
            --output-sequences {output.sequences}
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


rule subset_genoflu:
    input:
        genoflu="{data_source}/data/genoflu/results/results.tsv"
    output:
        genotypes="{data_source}/data/genoflu/genoflu_genotypes.tsv",
    shell:
        """
        csvtk cut -t \
            -f Strain,Genotype \
            {input.genoflu} \
            | csvtk rename -t \
                -f Strain,Genotype \
                -n strain,genoflu > {output.genotypes}
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
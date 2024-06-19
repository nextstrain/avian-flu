"""
This part of the workflow handles running Nextclade on the curated metadata
and sequences.
"""


DATASET_NAME = config["nextclade"]["dataset_name"]


rule get_nextclade_dataset:
    """Download Nextclade dataset"""
    output:
        dataset=f"data/nextclade/{DATASET_NAME}.zip",
    benchmark:
        "benchmarks/get_nextclade_dataset.txt"
    params:
        dataset_name=DATASET_NAME
    shell:
        """
        nextclade3 dataset get \
            --name={params.dataset_name:q} \
            --output-zip={output.dataset} \
            --verbose
        """


rule run_nextclade:
    input:
        dataset=f"data/nextclade/{DATASET_NAME}.zip",
        # The H5NX datasets should only be for the HA segment
        sequences="{data_source}/results/sequences_ha.fasta",
    output:
        nextclade="{data_source}/results/nextclade.tsv",
        alignment="{data_source}/results/alignment.fasta",
    benchmark:
        "{data_source}/benchmarks/run_nextclade.txt"
    shell:
        """
        nextclade3 run \
            {input.sequences} \
            --input-dataset {input.dataset} \
            --output-tsv {output.nextclade} \
            --output-fasta {output.alignment}
        """

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


rule join_metadata_and_nextclade:
    input:
        nextclade="{data_source}/results/nextclade.tsv",
        metadata="{data_source}/data/merged_segment_metadata.tsv",
        nextclade_field_map=config["nextclade"]["field_map"],
    output:
        metadata="{data_source}/results/metadata.tsv",
    params:
        # Making this param optional because we don't have curate pipeline for fauna data
        metadata_id_field=config.get("curate", {}).get("output_id_field", "strain"),
        nextclade_id_field=config["nextclade"]["id_field"],
    shell:
        """
        export SUBSET_FIELDS=`grep -v '^#' {input.nextclade_field_map} | awk '{{print $1}}' | tr '\n' ',' | sed 's/,$//g'`

        csvtk fix-quotes -t {input.nextclade} \
        | csvtk -t cut -f $SUBSET_FIELDS \
        | csvtk -t rename2 \
            -F \
            -f '*' \
            -p '(.+)' \
            -r '{{kv}}' \
            -k {input.nextclade_field_map} \
        | csvtk del-quotes -t \
        | tsv-join -H \
            --filter-file - \
            --key-fields {params.nextclade_id_field} \
            --data-fields {params.metadata_id_field} \
            --append-fields '*' \
            --write-all ? \
            {input.metadata} \
        | tsv-select -H --exclude {params.nextclade_id_field} \
            > {output.metadata}
        """

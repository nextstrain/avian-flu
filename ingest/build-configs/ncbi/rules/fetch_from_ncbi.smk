"""
This part of the workflow handles fetching sequences and metadata from NCBI.

REQUIRED INPUTS:

    None

OUTPUTS:

    ndjson = ncbi/data/ncbi.ndjson

"""


rule fetch_ncbi_dataset_package:
    params:
        ncbi_taxon_id=config["ncbi_taxon_id"],
        released_after=config["ncbi_released_after"],
        geo_location=config["ncbi_geo_location"],
    output:
        dataset_package=temp("ncbi/data/ncbi_dataset.zip"),
    # Allow retries in case of network errors
    retries: 5
    benchmark:
        "ncbi/benchmarks/fetch_ncbi_dataset_package.txt"
    shell:
        """
        datasets download virus genome taxon {params.ncbi_taxon_id:q} \
            --released-after {params.released_after:q} \
            --geo-location {params.geo_location:q} \
            --no-progressbar \
            --filename {output.dataset_package}
        """


rule extract_ncbi_dataset_sequences:
    input:
        dataset_package="ncbi/data/ncbi_dataset.zip",
    output:
        ncbi_dataset_sequences=temp("ncbi/data/ncbi_dataset_sequences.fasta"),
    benchmark:
        "ncbi/benchmarks/extract_ncbi_dataset_sequences.txt"
    shell:
        """
        unzip -jp {input.dataset_package} \
            ncbi_dataset/data/genomic.fna > {output.ncbi_dataset_sequences}
        """


rule format_ncbi_dataset_report:
    input:
        dataset_package="ncbi/data/ncbi_dataset.zip",
    output:
        ncbi_dataset_tsv=temp("ncbi/data/ncbi_dataset_report.tsv"),
    params:
        ncbi_datasets_fields=",".join(config["ncbi_datasets_fields"]),
    benchmark:
        "ncbi/benchmarks/format_ncbi_dataset_report.txt"
    shell:
        """
        dataformat tsv virus-genome \
            --package {input.dataset_package} \
            --fields {params.ncbi_datasets_fields:q} \
            --elide-header \
            | csvtk fix-quotes -Ht \
            | csvtk add-header -t -l -n {params.ncbi_datasets_fields:q} \
            | csvtk rename -t -f accession -n accession_version \
            | csvtk -t mutate -f accession_version -n accession -p "^(.+?)\." \
            | csvtk del-quotes -t \
            | tsv-select -H -f accession --rest last \
            > {output.ncbi_dataset_tsv}
        """


rule select_accessions_from_dataset_report:
    input:
        ncbi_dataset_tsv="ncbi/data/ncbi_dataset_report.tsv",
    output:
        genbank_accessions=temp("ncbi/data/genbank_accessions.txt"),
    params:
        id_column_name=config["ncbi_join_field"],
    shell:
        """
        tsv-select -H \
            -f {params.id_column_name} \
            {input.ncbi_dataset_tsv} \
            > {output.genbank_accessions}
        """


rule fetch_from_ncbi_entrez_with_accessions:
    input:
        genbank_accessions="ncbi/data/genbank_accessions.txt",
    output:
        entrez_metadata=temp("ncbi/data/entrez_metadata.tsv"),
    params:
        source_fields=" ".join(config["entrez_source_fields"]),
        id_column_name=config["ncbi_join_field"],
    log:
        "ncbi/logs/fetch_from_ncbi_entrez_with_accessions.txt",
    benchmark:
        "ncbi/benchmarks/fetch_from_ncbi_entrez_with_accessions.txt"
    shell:
        """
        ./build-configs/ncbi/bin/fetch_from_ncbi_entrez_with_accessions \
            --genbank-accessions {input.genbank_accessions} \
            --source-fields {params.source_fields} \
            --id-column-name {params.id_column_name} \
            --output {output.entrez_metadata} 2>> {log}
        """


rule join_entrez_with_datasets_metadata:
    input:
        entrez_metadata="ncbi/data/entrez_metadata.tsv",
        datasets_metadata="ncbi/data/ncbi_dataset_report.tsv",
    output:
        joined_metadata="ncbi/data/ncbi_metadata.tsv",
    params:
        ncbi_join_field=config["ncbi_join_field"]
    benchmark:
        "ncbi/benchmarks/join_entrez_with_datasets_metadata"
    shell:
        """
        tsv-join -H \
            --filter-file {input.entrez_metadata} \
            --key-fields {params.ncbi_join_field} \
            --append-fields '*' \
            --write-all '?' \
            {input.datasets_metadata} \
            > {output.joined_metadata}
        """


rule format_ncbi_datasets_ndjson:
    input:
        ncbi_dataset_sequences="ncbi/data/ncbi_dataset_sequences.fasta",
        ncbi_metadata="ncbi/data/ncbi_metadata.tsv",
    output:
        ndjson="ncbi/data/ncbi.ndjson",
    log:
        "ncbi/logs/format_ncbi_datasets_ndjson.txt",
    benchmark:
        "ncbi/benchmarks/format_ncbi_datasets_ndjson.txt"
    shell:
        """
        augur curate passthru \
            --metadata {input.ncbi_metadata} \
            --fasta {input.ncbi_dataset_sequences} \
            --seq-id-column accession_version \
            --seq-field sequence \
            --unmatched-reporting warn \
            --duplicate-reporting warn \
            2> {log} > {output.ndjson}
        """

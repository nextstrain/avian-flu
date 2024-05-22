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
    output:
        dataset_package=temp("ncbi/data/ncbi_dataset.zip"),
    # Allow retries in case of network errors
    retries: 5
    benchmark:
        "ncbi/benchmarks/fetch_ncbi_dataset_package.txt"
    shell:
        """
        datasets download virus genome taxon {params.ncbi_taxon_id:q} \
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


rule format_ncbi_datasets_ndjson:
    input:
        ncbi_dataset_sequences="ncbi/data/ncbi_dataset_sequences.fasta",
        ncbi_dataset_tsv="ncbi/data/ncbi_dataset_report.tsv",
    output:
        ndjson="ncbi/data/ncbi.ndjson",
    log:
        "ncbi/logs/format_ncbi_datasets_ndjson.txt",
    benchmark:
        "ncbi/benchmarks/format_ncbi_datasets_ndjson.txt"
    shell:
        """
        augur curate passthru \
            --metadata {input.ncbi_dataset_tsv} \
            --fasta {input.ncbi_dataset_sequences} \
            --seq-id-column accession_version \
            --seq-field sequence \
            --unmatched-reporting warn \
            --duplicate-reporting warn \
            2> {log} > {output.ndjson}
        """

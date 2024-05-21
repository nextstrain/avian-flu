"""
This part of the workflow handles fetching sequences and metadata from NCBI.

REQUIRED INPUTS:

    None

OUTPUTS:

    ndjson = ncbi/data/ncbi.ndjson

"""
import datetime


def _get_date_filter():
    """
    Construct the NCBI Virus collection date filter to use
    today as the max date.
    """
    min_date = config["ncbi_virus_min_collection_date"]
    max_date = datetime.datetime.today().strftime("%Y-%m-%d")
    return '{!tag=CollectionDate_dr}CollectionDate_dr:' + f'[{min_date} TO {max_date}]'


rule fetch_from_ncbi_virus:
    output:
        ncbi_virus_csv="ncbi/data/ncbi_virus.csv",
    params:
        github_repo="nextstrain/avian-flu",
        ncbi_taxon_id=config["ncbi_taxon_id"],
        ncbi_collection_date_filter=_get_date_filter(),
        ncbi_virus_filters=" ".join(f"{filter!r}" for filter in config["ncbi_virus_filters"]),
    shell:
        """
        ./build-configs/ncbi/bin/fetch-from-ncbi-virus \
            {params.ncbi_taxon_id} \
            {params.github_repo} \
            --filters {params.ncbi_collection_date_filter:q} {params.ncbi_virus_filters} \
            > {output.ncbi_virus_csv}
        """


rule select_accessions_from_ncbi_virus:
    input:
        ncbi_virus_csv="ncbi/data/ncbi_virus.csv",
    output:
        genbank_accessions=temp("ncbi/data/genbank_accessions.txt"),
    params:
        accession_column_name="accession_version",
    shell:
        """
        cat {input.ncbi_virus_csv} \
            | csvtk cut -U -f {params.accession_column_name} \
            > {output.genbank_accessions}
        """


rule fetch_ncbi_dataset_package:
    input:
        genbank_accessions="ncbi/data/genbank_accessions.txt",
    params:
        ncbi_taxon_id=config["ncbi_taxon_id"]
    output:
        dataset_package=temp("ncbi/data/ncbi_dataset.zip"),
    # Allow retries in case of network errors
    retries: 5
    benchmark:
        "ncbi/benchmarks/fetch_ncbi_dataset_package.txt"
    shell:
        """
        datasets download virus genome accession \
            --inputfile {input.genbank_accessions} \
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


rule format_ncbi_datasets_ndjson:
    input:
        ncbi_dataset_sequences="ncbi/data/ncbi_dataset_sequences.fasta",
        ncbi_virus_csv="ncbi/data/ncbi_virus.csv",
    output:
        ndjson="ncbi/data/ncbi.ndjson",
    log:
        "ncbi/logs/format_ncbi_datasets_ndjson.txt",
    benchmark:
        "ncbi/benchmarks/format_ncbi_datasets_ndjson.txt"
    shell:
        """
        augur curate passthru \
            --metadata {input.ncbi_virus_csv} \
            --fasta {input.ncbi_dataset_sequences} \
            --seq-id-column accession_version \
            --seq-field sequence \
            --unmatched-reporting warn \
            --duplicate-reporting warn \
            2> {log} > {output.ndjson}
        """

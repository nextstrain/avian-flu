from pathlib import Path


rule download_segment:
    output:
        sequences = "fauna/data/{segment}.fasta",
    params:
        fasta_fields = "strain virus accession collection_date region country division location host domestic_status subtype originating_lab submitting_lab authors PMID gisaid_clade h5_clade",
        output_dir = lambda wildcards, output: Path(output.sequences).parent,
        output_fstem = lambda wildcards, output: Path(output.sequences).stem,
    benchmark:
        "fauna/benchmarks/download_segment_{segment}.txt"
    shell:
        """
        python3 {path_to_fauna}/vdb/download.py \
            --database vdb \
            --virus avian_flu \
            --fasta_fields {params.fasta_fields} \
            --select  locus:{wildcards.segment} \
            --path {params.output_dir} \
            --fstem {params.output_fstem}
        """

rule parse_segment:
    input:
        sequences = "fauna/data/{segment}.fasta",
    output:
        sequences = "fauna/results/sequences_{segment}.fasta",
        metadata = "fauna/data/metadata_{segment}.tsv",
    params:
        fasta_fields =  "strain virus isolate_id date region country division location host domestic_status subtype originating_lab submitting_lab authors PMID gisaid_clade h5_clade",
        prettify_fields = "region country division location host originating_lab submitting_lab authors PMID"
    shell:
        """
        augur parse \
            --sequences {input.sequences} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --fields {params.fasta_fields} \
            --prettify-fields {params.prettify_fields}
        """

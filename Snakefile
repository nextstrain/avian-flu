SEGMENTS = ["pb2","pb1","pa","ha","np","na","mp","ns"]

path_to_fauna = '../fauna'

rule all:
    input:
        auspice_tree = expand("auspice/flu_avian_h5n1_{segment}_tree.json", segment=SEGMENTS),
        auspice_meta = expand("auspice/flu_avian_h5n1_{segment}_meta.json", segment=SEGMENTS)

rule files:
    params:
        dropped_strains = "config/dropped_strains.txt",
        reference = "config/h5n1_{segment}.ref.gb",
        colors = "config/colors.tsv",
        auspice_config = "config/auspice_config.json"

files = rules.files.params

def _get_min_length_by_wildcards(wildcards):
    len_dict = {"pb2": 2100, "pb1": 2100, "pa": 2000, "ha":1600, "np":1400,"na":1270,"mp":900, "ns":800}
    length = len_dict[wildcards.segment]
    return(length)

rule download:
    message: "Downloading sequences from fauna"
    output:
        sequences = "data/h5n1_{segment}.fasta"
    params:
        fasta_fields = "strain virus accession collection_date region country division location passage_category submitting_lab age gender"
    shell:
        """
        python3 {path_to_fauna}/vdb/download.py \
            --database vdb \
            --virus avian_flu \
            --fasta_fields {params.fasta_fields} \
            --select locus:{wildcards.segment} subtype:h5n1 \
            --path data \
            --fstem h5n1_{wildcards.segment}
        """

rule parse:
    message: "Parsing fasta into sequences and metadata"
    input:
        sequences = rules.download.output.sequences
    output:
        sequences = "results/sequences_h5n1_{segment}.fasta",
        metadata = "results/metadata_h5n1_{segment}.tsv"
    params:
        fasta_fields =  "strain virus isolate_id date region country division location passage authors age gender"
    shell:
        """
        augur parse \
            --sequences {input.sequences} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --fields {params.fasta_fields}
        """

rule filter:
    message:
        """
        Filtering to
          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
          - excluding strains in {input.exclude}
          - samples with missing region and country metadata
        """
    input:
        sequences = rules.parse.output.sequences,
        metadata = rules.parse.output.metadata,
        exclude = files.dropped_strains
    output:
        sequences = "results/filtered_{segment}.fasta"
    params:
        group_by = "country year month",
        sequences_per_group = 5,
        min_date = 1996,
        min_length = _get_min_length_by_wildcards,
        exclude_where = "host_species=Ferret country=? region=?"

    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --output {output.sequences} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-date {params.min_date} \
            --exclude-where {params.exclude_where} \
            --min-length {params.min_length}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = files.reference
    output:
        alignment = "results/aligned_{segment}.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree-raw_{segment}.nwk"
    params:
        method = "iqtree"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method {params.method}
        """

rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.align.output,
        metadata = rules.parse.output.metadata
    output:
        tree = "results/tree_{segment}.nwk",
        node_data = "results/branch_lengths_{segment}.json"
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_filter_iqd = 4
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "results/nt_muts_{segment}.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output {output.node_data} \
            --inference {params.inference}
        """

def _get_genes_by_wildcards(wildcards):
    if wildcards.segment == "ns":
        genes = ["NS1", "NS2"]
    elif wildcards.segment == "mp":
        genes = ["M1","M2"]
    else:
        genes = [wildcards.segment.upper()]
    return(genes)

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "results/aa_muts_{segment}.json"
    params:
        genes = _get_genes_by_wildcards
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data} \
            --genes {params.genes} \
        """

rule traits:
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata
    output:
        node_data = "results/traits_{segment}.json",
    params:
        columns = "country"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        colors = files.colors,
        auspice_config = files.auspice_config
    output:
        auspice_tree = "auspice/flu_avian_h5n1_{segment}_tree.json",
        auspice_meta = "auspice/flu_avian_h5n1_{segment}_meta.json"
    shell:
        """
        augur export \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} \
            --colors {input.colors} \
            --auspice-config {input.auspice_config} \
            --output-tree {output.auspice_tree} \
            --output-meta {output.auspice_meta}
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"

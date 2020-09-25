SUBTYPES = ["h5nx"] #,"h5n1", "h7n9", "h9n2"
SEGMENTS = ["pb2", "pb1", "pa", "ha", "np", "na", "mp", "ns"]

path_to_fauna = '../fauna'

rule all:
    input:
        auspice_json = expand("auspice/flu_avian_{subtype}_{segment}.json", subtype=SUBTYPES, segment=SEGMENTS)

rule files:
    params:
        dropped_strains = "config/dropped_strains_{subtype}.txt",
        reference = "config/reference_{subtype}_{segment}.gb",
        colors = "config/colors_{subtype}.tsv",
        lat_longs = "config/lat_longs_{subtype}.tsv",
        auspice_config = "config/auspice_config_{subtype}.json"

files = rules.files.params

def download_by(w):
    db = {'h5nx': '', 'h5n1': 'subtype:h5n1', 'h7n9': 'subtype:h7n9', 'h9n2': 'subtype:h9n2'}
    return(db[w.subtype])

def group_by(w):
    gb = {'h5nx': 'subtype country year','h5n1': 'country year', 'h7n9': 'division year', 'h9n2': 'country year'}
    return gb[w.subtype]

def sequences_per_group(w):
    spg = {'h5nx':20,'h5n1': '10', 'h7n9': '70', 'h9n2': '10'}
    return spg[w.subtype]

def min_length(w):
    len_dict = {"pb2": 2100, "pb1": 2100, "pa": 2000, "ha":1600, "np":1400, "na":1270, "mp":900, "ns":800}
    length = len_dict[w.segment]
    return(length)

def min_date(w):
    date = {'h5nx':'1996','h5n1': '1996', 'h7n9': '2013', 'h9n2': '1966'}
    return date[w.subtype]

def traits_columns(w):
    traits = {'h5nx':'region country','h5n1': 'region country', 'h7n9': 'country division', 'h9n2': 'region country'}
    return traits[w.subtype]

rule download:
    message: "Downloading sequences from fauna"
    output:
        sequences = "data/{subtype}_{segment}.fasta"
    params:
        fasta_fields = "strain virus accession collection_date region country division location host subtype originating_lab submitting_lab",
        download_by = download_by
    shell:
        """
        python3 {path_to_fauna}/vdb/download.py \
            --database test_vdb \
            --virus avian_flu \
            --fasta_fields {params.fasta_fields} \
            --select  {params.download_by} locus:{wildcards.segment} \
            --path data \
            --fstem {wildcards.subtype}_{wildcards.segment}
        """

rule parse:
    message: "Parsing fasta into sequences and metadata"
    input:
        sequences = rules.download.output.sequences
    output:
        sequences = "results/sequences_{subtype}_{segment}.fasta",
        metadata = "results/metadata_{subtype}_{segment}.tsv"
    params:
        fasta_fields =  "strain virus isolate_id date region country division location host subtype originating_lab submitting_lab",
        prettify_fields = "region country division location host originating_lab submitting_lab"
    shell:
        """
        augur parse \
            --sequences {input.sequences} \
            --output-sequences {output.sequences} \
            --output-metadata {output.metadata} \
            --fields {params.fasta_fields} \
            --prettify-fields {params.prettify_fields}
        """

rule filter:
    message:
        """
        Filtering to
          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
          - excluding strains in {input.exclude}
          - samples with missing region and country metadata
          - excluding strains prior to {params.min_date}
        """
    input:
        sequences = rules.parse.output.sequences,
        metadata = rules.parse.output.metadata,
        exclude = files.dropped_strains
    output:
        sequences = "results/filtered_{subtype}_{segment}.fasta"
    params:
        group_by = group_by,
        sequences_per_group = sequences_per_group,
        min_date = min_date,
        min_length = min_length,
        exclude_where = "host=laboratoryderived host=ferret host=unknown host=other country=? region=?"

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
            --min-length {params.min_length} \
            --non-nucleotide
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
        alignment = "results/aligned_{subtype}_{segment}.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps \
            --remove-reference \
            --nthreads 1
        """

rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree-raw_{subtype}_{segment}.nwk"
    params:
        method = "iqtree"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method {params.method} \
            --nthreads 1
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
        tree = "results/tree_{subtype}_{segment}.nwk",
        node_data = "results/branch-lengths_{subtype}_{segment}.json"
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
        node_data = "results/nt-muts_{subtype}_{segment}.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "results/aa-muts_{subtype}_{segment}.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data}
        """

rule traits:
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.parse.output.metadata
    output:
        node_data = "results/traits_{subtype}_{segment}.json",
    params:
        columns = traits_columns,
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence
        """

rule reconstruct_translations:
    message: "Reconstructing translations for DMS data view"
    input:
        tree = rules.refine.output.tree,
        node_data = "results/aa-muts_{subtype}_{segment}.json",
    output:
        aa_alignment = "results/aa-alignment_{subtype}_{segment}.fasta"
    shell:
        """
        augur reconstruct-sequences \
            --tree {input.tree} \
            --mutations {input.node_data} \
            --gene PB2 \
            --output {output.aa_alignment} \
            --internal-nodes
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
        lat_longs = files.lat_longs,
        auspice_config = files.auspice_config
    output:
        auspice_json = "auspice/flu_avian_{subtype}_{segment}.json"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts}\
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json}
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"

SUBTYPES = config.get('subtypes', ["h5nx", "h5n1", "h7n9", "h9n2"])
SEGMENTS = config.get('segments', ["pb2", "pb1", "pa", "ha","np", "na", "mp", "ns"])
TIME =     config.get('time',     ["all-time","2y"])
TARGET_SEQUENCES_PER_TREE = config.get('n_seqs', 3000)
S3_SRC = config.get('s3_src', "s3://nextstrain-data-private/files/workflows/avian-flu")

# The config option `same_strains_per_segment=True'` (e.g. supplied to snakemake via --config command line argument)
# will change the behaviour of the workflow to use the same strains for each segment. This is achieved via these steps:
# (1) Filter the HA segment as normal plus filter to those strains with 8 segments
# (2) Filter the other segments by simply force-including the same strains as (1)
SAME_STRAINS = bool(config.get('same_strains_per_segment', False))


path_to_fauna = '../fauna'


def all_targets():
    return [
        *expand("auspice/avian-flu_{subtype}_{segment}_{time}.json", subtype=[s for s in SUBTYPES if s in ["h5nx","h5n1"]], segment=SEGMENTS,time=TIME),
        *expand("auspice/avian-flu_{subtype}_{segment}_{time}.json", subtype=[s for s in SUBTYPES if s in ['h7n9', 'h9n2']], segment=SEGMENTS,time=[t for t in TIME if t=='all-time'])
    ]

rule all:
    input:
        auspice_json = all_targets()

rule test_target:
    """
    For testing purposes such as CI workflows.
    """
    input: "auspice/avian-flu_h5n1_ha_all-time.json"

rule files:
    params:
        dropped_strains = "config/dropped_strains_{subtype}.txt",
        include_strains = "config/include_strains_{subtype}_{time}.txt",
        reference = "config/reference_{subtype}_{segment}.gb",
        colors = "config/colors_{subtype}.tsv",
        lat_longs = "config/lat_longs_{subtype}.tsv",
        auspice_config = "config/auspice_config_{subtype}.json",
        clades_file = "clade-labeling/{subtype}-clades.tsv",
        description = "config/description.md"

files = rules.files.params

def subtypes_by_subtype_wildcard(wildcards):
    db = {
        'h5nx': ['h5n1', 'h5n2', 'h5n3', 'h5n4', 'h5n5', 'h5n6', 'h5n7', 'h5n8', 'h5n9'],
        'h5n1': ['h5n1'],
        'h7n9': ['h7n9'],
        'h9n2': ['h9n2'],
    }
    return(db[wildcards.subtype])

def metadata_by_wildcards(wildcards):
    if wildcards.subtype in ("h5n1", "h5nx"):
        return "results/metadata-with-clade_{subtype}.tsv"
    else:
        return "data/metadata_{subtype}.tsv"

def group_by(w):
    gb = {
        'h5nx': {'all-time': 'subtype country year', '2y': 'subtype region month host'},
        'h5n1': {'all-time': 'region country year', '2y': 'subtype region month host'},
        'h7n9': {'all-time': 'division year'},
        'h9n2': {'all-time': 'country year'}
        }
    return gb[w.subtype][w.time]

def min_length(w):
    len_dict = {"pb2": 2100, "pb1": 2100, "pa": 2000, "ha":1600, "np":1400, "na":1270, "mp":900, "ns":800}
    length = len_dict[w.segment]
    return(length)

def min_date(w):
    date = {
        'h5nx': {'all-time': '1996', '2y': '2Y'},
        'h5n1': {'all-time': '1996', '2y': '2Y'},
        'h7n9': {'all-time': '2013'},
        'h9n2': {'all-time': '1966'}
        }
    return date[w.subtype][w.time]

def traits_columns(w):
    traits = {'h5nx':'region','h5n1': 'region country', 'h7n9': 'country division', 'h9n2': 'region country'}
    return traits[w.subtype]

def clock_rate(w):
    clock_rates_h5nx = {
        'pb2': '--clock-rate 0.00287',
        'pb1': '--clock-rate 0.00267',
        'pa': '--clock-rate 0.00238',
        'ha': '--clock-rate 0.0048',
        'np': '--clock-rate 0.0022',
        'na': '--clock-rate 0.0028',
        'mp': '--clock-rate 0.0017',
        'ns': '--clock-rate 0.0017'
        }

    clock_rates_h5n1 = {
        'pb2': '--clock-rate 0.00287',
        'pb1': '--clock-rate 0.00264',
        'pa': '--clock-rate 0.00248',
        'ha': '--clock-rate 0.00455',
        'np': '--clock-rate 0.00252',
        'na': '--clock-rate 0.00349',
        'mp': '--clock-rate 0.00191',
        'ns': '--clock-rate 0.00249'
        }

    clock_rate = {
        'h5nx': {'all-time':'', '2y': clock_rates_h5nx[w.segment]},
        'h5n1': {'all-time':'', '2y': clock_rates_h5n1[w.segment]},
        'h7n9': {'all-time':''},
        'h9n2': {'all-time':''}
        }

    return clock_rate[w.subtype][w.time]


def clock_rate_std_dev(w):
    clock_rate_std_dev = {
        'h5nx': {'all-time': '', '2y': '--clock-std-dev 0.00211'},
        'h5n1': {'all-time': '', '2y': '--clock-std-dev 0.00211'},
        'h7n9': {'all-time': ''},
        'h9n2': {'all-time': ''}
        }

    return clock_rate_std_dev[w.subtype][w.time]


rule download_sequences:
    output:
        sequences = "data/{segment}/sequences.fasta.zst",
    params:
        s3_src=S3_SRC,
    shell:
        """
        aws s3 cp {params.s3_src:q}/{wildcards.segment}/sequences.fasta.zst {output.sequences}
        """

rule download_metadata:
    output:
        metadata = "data/metadata.tsv.zst",
    params:
        s3_src=S3_SRC,
    shell:
        """
        aws s3 cp {params.s3_src:q}/metadata.tsv.zst {output.metadata}
        """

rule filter_sequences_by_subtype:
    input:
        sequences="data/{segment}/sequences.fasta.zst",
        metadata="data/metadata.tsv.zst",
    output:
        sequences = "data/sequences_{subtype}_{segment}.fasta",
    params:
        subtypes=subtypes_by_subtype_wildcard,
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --query "subtype in {params.subtypes!r}" \
            --output-sequences {output.sequences}
        """

rule filter_metadata_by_subtype:
    input:
        metadata="data/metadata.tsv.zst",
    output:
        metadata = "data/metadata_{subtype}.tsv",
    params:
        subtypes=subtypes_by_subtype_wildcard,
    shell:
        """
        augur filter \
            --metadata {input.metadata} \
            --query "subtype in {params.subtypes!r}" \
            --output-metadata {output.metadata}
        """

rule add_h5_clade:
    message: "Adding in a column for h5 clade numbering"
    input:
        metadata = "data/metadata_{subtype}.tsv",
        clades_file = files.clades_file
    output:
        metadata= "results/metadata-with-clade_{subtype}.tsv"
    shell:
        """
        python clade-labeling/add-clades.py \
            --metadata {input.metadata} \
            --output {output.metadata} \
            --clades {input.clades_file}
        """

def _filter_params(wildcards, input, output, threads, resources):
    """
    Generate the arguments to `augur filter`. When we are running independent analyses
    (i.e. not using the same strains for each segment), then we generate a full set of
    filter parameters here.
    When we are using the same sequences for each segment, then for HA we use a full
    filter call and for the rest of the segments we filter to the strains chosen for HA
    """
    # For non-HA segments when we are using the SAME_STRAINS for all segments
    # we have a simple filtering approach: match what we're using for HA!
    # We also include the "force-include" list; 99% of the time these strains will already be in
    # the HA strain list (as they were force-included there too) but there may be rare occasions
    # where we force-include a strain which does not have a HA sequence.
    if input.strains:
        # some basic error checking to guard against inadvertent changes in the future
        if not (SAME_STRAINS and wildcards.segment!='ha'):
            raise Exception("A strains input should only be present for SAME_STRAINS + HA!")
        return f"--exclude-all --include {input.strains} {input.include}"

    # If SAME_STRAINS (and due to the above conditional we have the HA segment at this point)
    # then we want to restrict to strains present in all 8 segments. Note that force-included
    # strains may not have all segments, but that's preferable to filtering them out.
    restrict_n_segments = f"n_segments!={len(SEGMENTS)}" if SAME_STRAINS else ''

    # formulate our typical filtering parameters
    cmd  = f" --group-by {group_by(wildcards)}"
    cmd += f" --subsample-max-sequences {TARGET_SEQUENCES_PER_TREE}"
    cmd += f" --min-date {min_date(wildcards)}"
    cmd += f" --include {input.include}"
    cmd += f" --exclude-where host=laboratoryderived host=ferret host=unknown host=other country=? region=? gisaid_clade=3C.2 {restrict_n_segments}"
    cmd += f" --min-length {min_length(wildcards)}"
    cmd += f" --non-nucleotide"
    return cmd

rule filter:
    input:
        sequences = "data/sequences_{subtype}_{segment}.fasta",
        metadata = metadata_by_wildcards,
        exclude = files.dropped_strains,
        include = files.include_strains,
        strains = lambda w: f"results/filtered_{w.subtype}_ha_{w.time}.txt" if (SAME_STRAINS and w.segment!='ha') else [],
    output:
        sequences = "results/filtered_{subtype}_{segment}_{time}.fasta",
        strains = "results/filtered_{subtype}_{segment}_{time}.txt",
    params:
        args = _filter_params,
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --output-sequences {output.sequences} \
            --output-strains {output.strains} \
            {params.args}
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
        alignment = "results/aligned_{subtype}_{segment}_{time}.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --remove-reference \
            --nthreads 1
        """


rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree-raw_{subtype}_{segment}_{time}.nwk"
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
        metadata = metadata_by_wildcards,
    output:
        tree = "results/tree_{subtype}_{segment}_{time}.nwk",
        node_data = "results/branch-lengths_{subtype}_{segment}_{time}.json"
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_filter_iqd = 4,
        clock = clock_rate,
        clock_std_dev = clock_rate_std_dev
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
            {params.clock} \
            {params.clock_std_dev} \
            --clock-filter-iqd {params.clock_filter_iqd}
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.align.output
    output:
        node_data = "results/nt-muts_{subtype}_{segment}_{time}.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}\
            --keep-ambiguous
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = files.reference
    output:
        node_data = "results/aa-muts_{subtype}_{segment}_{time}.json"
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
        metadata = metadata_by_wildcards,
    output:
        node_data = "results/traits_{subtype}_{segment}_{time}.json",
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

rule cleavage_site:
    message: "determining sequences that harbor furin cleavage sites"
    input:
        alignment = "results/aligned_{subtype}_ha_{time}.fasta"
    output:
        cleavage_site_annotations = "results/cleavage-site_{subtype}_ha_{time}.json",
        cleavage_site_sequences = "results/cleavage-site-sequences_{subtype}_ha_{time}.json"
    shell:
        """
        python scripts/annotate-ha-cleavage-site.py \
            --alignment {input.alignment} \
            --furin_site_motif {output.cleavage_site_annotations} \
            --cleavage_site_sequence {output.cleavage_site_sequences}
        """

rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = metadata_by_wildcards,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        fcs = rules.cleavage_site.output.cleavage_site_annotations,
        cleavage_site_sequences = rules.cleavage_site.output.cleavage_site_sequences,
        colors = files.colors,
        lat_longs = files.lat_longs,
        auspice_config = files.auspice_config,
        description = files.description
    output:
        auspice_json = "auspice/avian-flu_{subtype}_{segment}_{time}.json"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} {input.fcs} {input.cleavage_site_sequences}\
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --description {input.description} \
            --include-root-sequence-inline \
            --output {output.auspice_json}
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"

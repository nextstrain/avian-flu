include: "rules/common.smk"
include: "rules/cattle-flu.smk"

SUBTYPES = config.get('subtypes', ["h5nx", "h5n1", "h7n9", "h9n2"])
SEGMENTS = config.get('segments', ["pb2", "pb1", "pa", "ha","np", "na", "mp", "ns"])
TIME =     config.get('time',     ["all-time","2y"])
TARGET_SEQUENCES_PER_TREE = config.get('n_seqs', 3000)

# The config option `same_strains_per_segment=True'` (e.g. supplied to snakemake via --config command line argument)
# will change the behaviour of the workflow to use the same strains for each segment. This is achieved via these steps:
# (1) Filter the HA segment as normal plus filter to those strains with 8 segments
# (2) Filter the other segments by simply force-including the same strains as (1)
SAME_STRAINS = bool(config.get('same_strains_per_segment', False))

# constrain the wildcards to not include `_` which we use to separate "parts" of filenames (where a part may be a wildcard itself)
wildcard_constraints:
    subtype = "[^_]+",
    segment = "[^_]+",
    time = "[^_]+",

def all_targets():
    return [
        *expand("auspice/avian-flu_{subtype}_{segment}_{time}.json", subtype=[s for s in SUBTYPES if s in ["h5nx","h5n1"]], segment=SEGMENTS,time=TIME),
        *expand("auspice/avian-flu_{subtype}_{segment}_{time}.json", subtype=[s for s in SUBTYPES if s in ['h7n9', 'h9n2']], segment=SEGMENTS,time=[t for t in TIME if t=='all-time'])
    ]

rule all:
    input:
        auspice_json = all_targets()

# This must be after the `all` rule above since it depends on its inputs
include: "rules/deploy.smk"

rule test_target:
    """
    For testing purposes such as CI workflows.
    """
    input: "auspice/avian-flu_h5n1_ha_all-time.json"

rule files:
    params:
        dropped_strains = lambda w: "config/dropped_strains_{subtype}.txt" if not w.subtype=='h5n1-cattle-outbreak' else "config/dropped_strains_h5n1.txt",
        include_strains = lambda w: "config/include_strains_{subtype}_{time}.txt" if not w.subtype=='h5n1-cattle-outbreak' else "config/include_strains_{subtype}.txt",
        reference = lambda w: "config/reference_{subtype}_{segment}.gb" if not w.subtype=='h5n1-cattle-outbreak' else "config/reference_h5n1_{segment}.gb",
        colors = lambda w: "config/colors_{subtype}.tsv" if not w.subtype=='h5n1-cattle-outbreak' else "config/colors_h5n1.tsv",
        # TODO - Augur 24.4.0 includes extensive lat-longs by default - can we drop the following avian-flu specific ones?
        lat_longs = lambda w: "config/lat_longs_{subtype}.tsv" if not w.subtype=='h5n1-cattle-outbreak' else "config/lat_longs_h5n1.tsv",
        auspice_config = "config/auspice_config_{subtype}.json",
        clades_file = lambda w: "clade-labeling/{subtype}-clades.tsv" if not w.subtype=='h5n1-cattle-outbreak' else "clade-labeling/h5n1-clades.tsv",
        description = lambda w: "config/description.md" if not w.subtype=='h5n1-cattle-outbreak' else "config/description_{subtype}.md",

files = rules.files.params

def metadata_by_wildcards(wildcards):
    # H5 builds have extra clade-level metadata added to the metadata TSV.
    # We may move this to a node-data JSON which would simplify the snakemake logic
    # a bit -- see <https://github.com/nextstrain/avian-flu/issues/25>
    if wildcards.subtype in ("h5n1", "h5nx", "h5n1-cattle-outbreak"):
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
    if w.subtype == 'h5n1-cattle-outbreak':
        return ''
    return gb[w.subtype][w.time]

def min_length(w):
    len_dict = {"pb2": 2100, "pb1": 2100, "pa": 2000, "ha":1600, "np":1400, "na":1270, "mp":900, "ns":800}
    length = len_dict[w.segment]
    return(length)

def min_date(w):
    if w.subtype == 'h5n1-cattle-outbreak':
        return "2024"
    date = {
        'h5nx': {'all-time': '1996', '2y': '2Y'},
        'h5n1': {'all-time': '1996', '2y': '2Y'},
        'h7n9': {'all-time': '2013'},
        'h9n2': {'all-time': '1966'}
        }
    return date[w.subtype][w.time]

def traits_columns(w):
    traits = {'h5nx':'region','h5n1': 'region country', 'h7n9': 'country division', 'h9n2': 'region country'}
    traits['h5n1-cattle-outbreak'] = traits['h5n1']
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
        'h9n2': {'all-time':''},
        'h5n1-cattle-outbreak': {'default': clock_rates_h5n1[w.segment]}
        }

    return clock_rate[w.subtype][w.time]


def clock_rate_std_dev(w):
    clock_rate_std_dev = {
        'h5nx': {'all-time': '', '2y': '--clock-std-dev 0.00211'},
        'h5n1': {'all-time': '', '2y': '--clock-std-dev 0.00211'},
        'h7n9': {'all-time': ''},
        'h9n2': {'all-time': ''},
        'h5n1-cattle-outbreak': {'default': '--clock-std-dev 0.00211'}
        }

    return clock_rate_std_dev[w.subtype][w.time]

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

    ## By default we filter out unknown country/region, but skip this for the cattle-flu outbreak just
    ## in case there's something interesting with limited metadata
    require_location = ' country=? region=? ' if wildcards.subtype!='h5n1-cattle-outbreak' else ''

    grouping = f" --group-by {group_by(wildcards)}" if group_by(wildcards) else ""

    # formulate our typical filtering parameters
    cmd  = grouping
    cmd += f" --subsample-max-sequences {TARGET_SEQUENCES_PER_TREE}"
    cmd += f" --min-date {min_date(wildcards)}"
    cmd += f" --include {input.include}"
    cmd += f" --exclude-where host=laboratoryderived host=ferret host=unknown host=other host=host {require_location} gisaid_clade=3C.2 {restrict_n_segments}"
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
    threads:
        4
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --remove-reference \
            --nthreads {threads}
        """


rule tree:
    message: "Building tree"
    input:
        alignment = rules.align.output.alignment
    output:
        tree = "results/tree-raw_{subtype}_{segment}_{time}.nwk"
    params:
        method = "iqtree"
    threads:
        4
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method {params.method} \
            --nthreads {threads}
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

def refined_tree(w):
    """
    Return the refined tree to be used for export, traits, ancestry reconstruction etc
    The cattle-flu build introduces an additional step beyond `augur refine`, which is
    why this function exists.
    """
    if w.subtype=='h5n1-cattle-outbreak':
        return "results/tree_{subtype}_{segment}_{time}_outbreak-clade.nwk"
    return "results/tree_{subtype}_{segment}_{time}.nwk"

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = refined_tree,
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
        tree = refined_tree,
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
        tree = refined_tree,
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

def export_node_data_files(wildcards):
    nd = [
        rules.refine.output.node_data,
        rules.traits.output.node_data,
        rules.ancestral.output.node_data,
        rules.translate.output.node_data,
        rules.cleavage_site.output.cleavage_site_annotations,
        rules.cleavage_site.output.cleavage_site_sequences,
    ]
    if wildcards.subtype=="h5n1-cattle-outbreak":
        nd.append("results/tree_{subtype}_{segment}_{time}_outbreak-clade.json")
    return nd


def additional_export_config(wildcards):
    args = ""

    if wildcards.subtype == "h5n1-cattle-outbreak":
        # The auspice-config is for the whole genome analysis, so override the title
        segment = wildcards.segment.upper()
        args += f"--title 'Ongoing influenza A/H5N1 cattle outbreak in North America ({segment} segment)'"

    return args

rule export:
    """
    Export the files into results/ and then use a subsequent rule to move these to the
    auspice/ directory
    """
    input:
        tree = refined_tree,
        metadata = metadata_by_wildcards,
        node_data = export_node_data_files,
        colors = files.colors,
        lat_longs = files.lat_longs,
        auspice_config = files.auspice_config,
        description = files.description
    output:
        auspice_json = "results/avian-flu_{subtype}_{segment}_{time}.json"
    params:
        additional_config = additional_export_config
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --description {input.description} \
            --include-root-sequence-inline \
            {params.additional_config} \
            --output {output.auspice_json}
        """

def auspice_name_to_wildcard_name(wildcards):
    """
    Used to link Auspice JSONs filenames to their intermediate filename which includes all wildcards.
    Examples:
    1. subtype + segment + time in their filename / URL,
        e.g. "avian-flu_h5n1_ha_2y.json" (nextstrain.org/avian-flu/h5n1/ha/2y)
        maps to subtype=h5n1, segment=ha, time=2y
    2. subtype + segment in their filename / URL,
        e.g. "avian-flu_h5n1-cattle-outbreak_ha.json" (nextstrain.org/avian-flu/h5n1-cattle-outbreak/ha)
        maps to subtype=h5n1-cattle-outbreak, segment=ha, time=default
    """
    parts = wildcards.parts.split("_")
    if len(parts)==3:
        [subtype, segment, time] = parts
        assert segment!='genome', "Genome builds are not available for this build"
        return f"results/avian-flu_{subtype}_{segment}_{time}.json"
    if len(parts)==2:
        [subtype, segment] = parts
        assert subtype=='h5n1-cattle-outbreak', "Only h5n1 builds produce an Auspice dataset without a time component in the filename"
        return f"results/avian-flu_{subtype}_{segment}_default.json"
    raise Exception("Auspice JSON filename requested with an unexpected number of (underscore-separated) parts")


rule rename_auspice_datasets:
    """
    This allows us to create files in auspice/ which mirror the intended URL structure rather than
    the wildcard structure we use in the workflow.
    """
    input:
        json = auspice_name_to_wildcard_name
    output:
        json = "auspice/avian-flu_{parts}.json"
    wildcard_constraints:
        timepart = ".*"
    shell:
        """
        cp {input.json} {output.json}
        """


rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"

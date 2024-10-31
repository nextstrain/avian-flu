# constrain the wildcards to not include `_` which we use to separate "parts" of filenames (where a part may be a wildcard itself)
wildcard_constraints:
    subtype = "[^_/]+",
    segment = "[^_/]+",
    time = "[^_/]+",

# defined before extra rules `include`d as they reference this constant
SEGMENTS = ["pb2", "pb1", "pa", "ha","np", "na", "mp", "ns"]

# The config option `same_strains_per_segment=True'` (e.g. supplied to snakemake via --config command line argument)
# will change the behaviour of the workflow to use the same strains for each segment. This is achieved via these steps:
# (1) Filter the HA segment as normal plus filter to those strains with 8 segments
# (2) Filter the other segments by simply force-including the same strains as (1)
SAME_STRAINS = bool(config.get('same_strains_per_segment', False))

NEXTSTRAIN_PUBLIC_BUCKET = "s3://nextstrain-data/"
S3_SRC = config.get('s3_src', {})
LOCAL_INGEST = config.get('local_ingest', None)

def sanity_check_config():
    assert LOCAL_INGEST or S3_SRC, "The config must define either 's3_src' or 'local_ingest'"
    # NOTE: we could relax the following exclusivity of S3_SRC and LOCAL_INGEST
    # if we want to use `--config local_ingest=gisaid` overrides.
    assert not (S3_SRC and LOCAL_INGEST), "The config defined both 'local_ingest' and 's3_src', which are mutually exclusive"
    if S3_SRC:
        assert isinstance(S3_SRC, dict) and all([k in S3_SRC for k in ("name", "sequences", "metadata")]), \
            "Config 's3_src' must be a dict with 'name', 'sequences' and 'metadata' keys"

sanity_check_config()

def collect_builds():
    """
    iteratively create workflow targets from config.builds for the `all` rule
    you can over-ride this by specifying targets (filenames) on the command line
    """
    targets = []
    for subtype,times in config.get('builds', {}).items():
        for segment in config.get('segments', []):
            if len(times):
                for time in times:
                    targets.append(f"auspice/avian-flu_{subtype}_{segment}_{time}.json")
            else:
                targets.append(f"auspice/avian-flu_{subtype}_{segment}.json")
    return targets

rule all:
    input:
        auspice_json = collect_builds()

# This must be after the `all` rule above since it depends on its inputs
include: "rules/deploy.smk"

rule test_target:
    """
    For testing purposes such as CI workflows.
    """
    input: "auspice/avian-flu_h5n1_ha_all-time.json"

rule files:
    params:
        dropped_strains = config['dropped_strains'],
        include_strains = config['include_strains'],
        reference = config['reference'],
        colors = config['colors'],
        # TODO - Augur 24.4.0 includes extensive lat-longs by default - can we drop the following avian-flu specific ones?
        lat_longs = config['lat_longs'],
        auspice_config = config['auspice_config'],
        clades_file = config['clades_file'],
        description = config['description'],

files = rules.files.params


def subtypes_by_subtype_wildcard(wildcards):
    db = {
        'h5nx': ['h5n1', 'h5n2', 'h5n3', 'h5n4', 'h5n5', 'h5n6', 'h5n7', 'h5n8', 'h5n9'],
        'h5n1': ['h5n1'],
        'h7n9': ['h7n9'],
        'h9n2': ['h9n2'],
    }
    db['h5n1-cattle-outbreak'] = [*db['h5nx']]
    assert wildcards.subtype in db, (f"Subtype {wildcards.subtype!r} is not defined in the snakemake function "
        "`subtypes_by_subtype_wildcard` -- is there a typo in the subtype you are targetting?")
    return(db[wildcards.subtype])

rule download_sequences:
    output:
        sequences = f"data/{S3_SRC.get('name', None)}/sequences_{{segment}}.fasta",
    params:
        address=lambda w: S3_SRC.get('sequences', None).format(segment=w.segment),
        no_sign_request=lambda w: "--no-sign-request" if S3_SRC.get('sequences', "").startswith(NEXTSTRAIN_PUBLIC_BUCKET) else ""
    shell:
        """
        aws s3 cp {params.no_sign_request:q} {params.address:q} - | zstd -d > {output.sequences}
        """

rule download_metadata:
    output:
        metadata = f"data/{S3_SRC.get('name', None)}/metadata.tsv",
    params:
        address=S3_SRC.get('metadata', None),
        no_sign_request=lambda w: "--no-sign-request" if S3_SRC.get('metadata', "").startswith(NEXTSTRAIN_PUBLIC_BUCKET) else ""
    shell:
        """
        aws s3 cp {params.no_sign_request:q} {params.address:q} - | zstd -d > {output.metadata}
        """


def input_metadata(wildcards):
    if S3_SRC:
        return f"data/{S3_SRC['name']}/metadata.tsv",
    elif LOCAL_INGEST:
        return f"ingest/{LOCAL_INGEST}/results/metadata.tsv",
    raise Exception() # already caught by `sanity_check_config` above, this is just being cautious

def input_sequences(wildcards):
    if S3_SRC:
        return f"data/{S3_SRC['name']}/sequences_{wildcards.segment}.fasta",
    elif LOCAL_INGEST:
        return f"ingest/{LOCAL_INGEST}/results/sequences_{wildcards.segment}.fasta"
    raise Exception() # already caught by `sanity_check_config` above, this is just being cautious


rule filter_sequences_by_subtype:
    input:
        sequences = input_sequences,
        metadata = input_metadata,
    output:
        sequences = "results/{subtype}/{segment}/sequences.fasta",
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
        metadata = input_metadata,
    output:
        metadata = "results/{subtype}/metadata.tsv",
    params:
        subtypes=subtypes_by_subtype_wildcard,
    shell:
        """
        augur filter \
            --metadata {input.metadata} \
            --query "subtype in {params.subtypes!r}" \
            --output-metadata {output.metadata}
        """

def metadata_by_wildcards(wildcards):
    # H5 builds have extra clade-level metadata added to the metadata TSV.
    # We may move this to a node-data JSON which would simplify the snakemake logic
    # a bit -- see <https://github.com/nextstrain/avian-flu/issues/25>
    if wildcards.subtype in ("h5n1", "h5nx", "h5n1-cattle-outbreak"):
        return "results/{subtype}/metadata-with-clade.tsv"
    else:
        return "results/{subtype}/metadata.tsv",


def get_config(rule_name, rule_key, wildcards, segment=None, fallback="FALLBACK"):
    assert rule_name in config, f"Config missing top-level {rule_name} key"
    assert rule_key in config[rule_name], f"Config missing entry for {rule_name}.{rule_key}"
    try:
        return config[rule_name][rule_key][wildcards.subtype][wildcards.time]
    except KeyError:
        assert fallback in config[rule_name][rule_key], f"config.{rule_name!r}.{rule_key!r} either needs " \
            f"an entry for {wildcards.subtype!r}.{wildcards.time!r} added or a (default) {fallback!r} key."
        return config[rule_name][rule_key][fallback]

def refine_clock_rates(w):
    info = get_config('refine', 'clock_rates', w)

    if w.segment == 'genome':
        # calculate the genome rate via a weighted average of the segment rates
        assert 'genome' not in info, ("This snakemake pipeline is currently set up to calculate the genome clock rate "
            "based on the segment rates, however you have provided a genome rate in the config.")
        try:
            segment_lengths= get_config('refine', 'segment_lengths', w)
        except AssertionError as e:
            # py11 brings e.add_note() which is nicer
            e.args = (*e.args, "NOTE: For segment=genome we require the segment_lengths to be defined as we use them to calculate the clock rate")
            raise
        mean = sum([info[seg][0]*length for seg,length in segment_lengths.items()])/sum(segment_lengths.values())
        stdev = mean/2
        return f"--clock-rate {mean} --clock-std-dev {stdev}"

    assert w.segment in info, 'config.refine.clock_rates: data must be provided for each segment. Use "" for inference.'
    if info[w.segment] == "":
        return ""

    assert isinstance(info[w.segment], list), "The clock rates for {w.subtype!r} {w.time!r} {w.segment!r} must be a list of (rate, std-dev)"
    assert len(info[w.segment])==2, "The clock rates for {w.subtype!r} {w.time!r} {w.segment!r} must be a list of (rate, std-dev)"
    return f"--clock-rate {info[w.segment][0]} --clock-std-dev {info[w.segment][1]}"

def refine_clock_filter(w):
    filter = get_config('refine', 'clock_filter_iqd', w)
    return f"--clock-filter-iqd {filter}" if filter else ""


rule add_h5_clade:
    message: "Adding in a column for h5 clade numbering"
    input:
        metadata = "results/{subtype}/metadata.tsv",
        clades_file = files.clades_file
    output:
        metadata= "results/{subtype}/metadata-with-clade.tsv"
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

    assert wildcards.segment!='genome', "Don't use this function for genome builds!"

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

    exclude_where = get_config('filter', 'exclude_where', wildcards)
    # If SAME_STRAINS (and due to the above conditional we have the HA segment at this point)
    # then we want to restrict to strains present in all 8 segments. Note that force-included
    # strains may not have all segments, but that's preferable to filtering them out.
    if SAME_STRAINS:
        exclude_where += f" n_segments!={len(SEGMENTS)}"


    cmd = ""

    group_by_value = get_config('filter', 'group_by', wildcards)
    cmd += f" --group-by {group_by_value}" if group_by_value else ""

    cmd += f" --subsample-max-sequences {config['target_sequences_per_tree']}"
    cmd += f" --min-date {get_config('filter', 'min_date', wildcards)}"
    cmd += f" --include {input.include}"
    cmd += f" --exclude-where {exclude_where}"
    cmd += f" --min-length {get_config('filter', 'min_length', wildcards)[wildcards.segment]}"
    cmd += f" --non-nucleotide"
    return cmd

rule filter:
    input:
        sequences = "results/{subtype}/{segment}/sequences.fasta",
        metadata = metadata_by_wildcards,
        exclude = files.dropped_strains,
        include = files.include_strains,
        strains = lambda w: f"results/{w.subtype}/ha/{w.time}/filtered.txt" if (SAME_STRAINS and w.segment!='ha') else [],
    output:
        sequences = "results/{subtype}/{segment}/{time}/filtered.fasta",
        strains = "results/{subtype}/{segment}/{time}/filtered.txt",
        metadata = "results/{subtype}/{segment}/{time}/metadata.tsv",
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
            --output-metadata {output.metadata} \
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
        alignment = "results/{subtype}/{segment}/{time}/aligned.fasta"
    wildcard_constraints:
        # for genome builds we don't use this rule; see `rule join_segments`
        segment = "|".join(seg for seg in SEGMENTS)
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
        # Note that the alignment input may come from `rule align` or `rule join_segments`
        alignment = rules.align.output.alignment
    output:
        tree = "results/{subtype}/{segment}/{time}/tree-raw.nwk"
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


def refine_root(wildcards):
    root = get_config('refine', 'genome_root', wildcards) \
        if wildcards.segment=='genome' \
        else get_config('refine', 'root', wildcards)
    return f"--root {root}" if root else ""

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
        metadata = rules.filter.output.metadata,
    output:
        tree = "results/{subtype}/{segment}/{time}/tree.nwk",
        node_data = "results/{subtype}/{segment}/{time}/branch-lengths.json"
    params:
        coalescent = config['refine']['coalescent'],
        date_inference = config['refine']['date_inference'],
        clock_rates = refine_clock_rates,
        clock_filter = refine_clock_filter,
        root = refine_root,
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            {params.root} \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference} \
            {params.clock_rates} \
            {params.clock_filter}
        """


def refined_tree(w):
    """
    Return the refined tree to be used for export, traits, ancestry reconstruction etc
    The cattle-flu build introduces an additional step beyond `augur refine`, which is
    why this function exists.
    """
    if w.subtype=='h5n1-cattle-outbreak' and w.segment!='genome':
        return "results/{subtype}/{segment}/{time}/tree_outbreak-clade.nwk"
    return "results/{subtype}/{segment}/{time}/tree.nwk",

def ancestral_root_seq(wildcards):
    root_seq = get_config('ancestral', 'genome_root_seq', wildcards) \
        if wildcards.segment=='genome' \
        else get_config('ancestral', 'root_seq', wildcards)
    return f"--root-sequence {root_seq}" if root_seq else ""

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = refined_tree,
        alignment = rules.align.output
    output:
        node_data = "results/{subtype}/{segment}/{time}/nt-muts.json"
    params:
        inference = config['ancestral']['inference'],
        root_seq = ancestral_root_seq,
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference} \
            {params.root_seq} \
            --keep-ambiguous
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = refined_tree,
        node_data = rules.ancestral.output.node_data,
        reference = lambda w: config['genome_reference'] if w.segment=='genome' else files.reference
    output:
        node_data = "results/{subtype}/{segment}/{time}/aa-muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data}
        """

def traits_params(wildcards):
    columns = get_config('traits', 'genome_columns', wildcards) \
        if wildcards.segment=='genome' \
        else get_config('traits', 'columns', wildcards)

    bias = get_config('traits', 'genome_sampling_bias_correction', wildcards) \
        if wildcards.segment=='genome' \
        else get_config('traits', 'sampling_bias_correction', wildcards)

    confidence = get_config('traits', 'confidence', wildcards)

    args = f"--columns {columns}"
    if bias:
        args += f" --sampling-bias-correction {bias}"
    if confidence:
        args += f" --confidence"
    return args


rule traits:
    input:
        tree = refined_tree,
        metadata = "results/{subtype}/{segment}/{time}/metadata.tsv",
    output:
        node_data = "results/{subtype}/{segment}/{time}/traits.json",
    params:
        info = traits_params,
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            {params.info}
        """


rule cleavage_site:
    """
    Cleavage sites are only computed for HA and are unchanged depending on the build. As such, they are not
    parameterised in the config file
    """
    input:
        alignment = "results/{subtype}/ha/{time}/aligned.fasta"
    output:
        cleavage_site_annotations = "results/{subtype}/ha/{time}/cleavage-site.json",
        cleavage_site_sequences = "results/{subtype}/ha/{time}/cleavage-site-sequences.json"
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

    if wildcards.subtype=="h5n1-cattle-outbreak" and wildcards.segment!='genome':
        nd.append(rules.prune_tree.output.node_data)
    return nd


def additional_export_config(wildcards):
    args = ""

    title_overrides = get_config('export', 'genome_title', wildcards) \
        if wildcards.segment=='genome' \
        else get_config('export', 'title', wildcards)
    if title_overrides:
        args += ("--title '" +
            title_overrides.format(segment=wildcards.segment.upper(), subtype=wildcards.subtype.upper(), time=wildcards.time) +
            "'")

    return args

rule auspice_config:
    """
    Depending on the build we may make wildcard-dependent modifications to the auspice config JSON.
    If we implement config overlays in augur this rule will become unnecessary.
    """
    input:
        auspice_config = files.auspice_config,
    output:
        auspice_config = "results/{subtype}/{segment}/{time}/auspice-config.json",
    run:
        import json
        with open(input.auspice_config) as fh:
            auspice_config = json.load(fh)
        if wildcards.subtype == "h5n1-cattle-outbreak":
            if wildcards.segment == "genome":
                auspice_config['display_defaults']['distance_measure'] = "num_date"
                division_idx = next((i for i,c in enumerate(auspice_config['colorings']) if c['key']=='division'), None)
                assert division_idx!=None, "Auspice config did not have a division coloring!"
                auspice_config['colorings'].insert(division_idx+1, {
                    "key": "division_metadata",
                    "title": auspice_config['colorings'][division_idx]["title"] + " (metadata)",
                    "type": "categorical",
                })
                auspice_config['colorings'][division_idx]["title"] += " (inferred)"
            else:
                auspice_config['display_defaults']['distance_measure'] = "div"
        with open(output.auspice_config, 'w') as fh:
            json.dump(auspice_config, fh, indent=2)


rule export:
    """
    Export the files into results/ and then use a subsequent rule to move these to the
    auspice/ directory
    """
    input:
        tree = refined_tree,
        metadata = rules.filter.output.metadata,
        node_data = export_node_data_files,
        colors = files.colors,
        lat_longs = files.lat_longs,
        auspice_config = rules.auspice_config.output.auspice_config,
        description = files.description
    output:
        auspice_json = "results/{subtype}/{segment}/{time}/auspice-dataset.json"
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
        return f"results/{subtype}/{segment}/{time}/auspice-dataset.json"
    if len(parts)==2:
        [subtype, segment] = parts
        assert subtype=='h5n1-cattle-outbreak', "Only h5n1 builds produce an Auspice dataset without a time component in the filename"
        return f"results/{subtype}/{segment}/default/auspice-dataset.json"
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


for rule_file in config.get('custom_rules', []):
    include: rule_file

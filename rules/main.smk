# constrain the wildcards to not include `_` which we use to separate "parts" of filenames (where a part may be a wildcard itself)
wildcard_constraints:
    subtype = "[^_/]+",
    segment = "[^_/]+",
    time = "[^_/]+",

# defined before extra rules `include`d as they reference this constant
SEGMENTS = ["pb2", "pb1", "pa", "ha","np", "na", "mp", "ns"]
#SUBTYPES = ["h5n1", "h5nx", "h7n9", "h9n2"]

# The config option `same_strains_per_segment=True'` (e.g. supplied to snakemake via --config command line argument)
# will change the behaviour of the workflow to use the same strains for each segment. This is achieved via these steps:
# (1) Filter the HA segment as normal plus filter to those strains with 8 segments
# (2) Filter the other segments by simply force-including the same strains as (1)
SAME_STRAINS = bool(config.get('same_strains_per_segment', False))

NEXTSTRAIN_PUBLIC_BUCKET = "s3://nextstrain-data/"

rule all:
    input:
        auspice_json = expand_target_patterns()


# This must be after the `all` rule above since it depends on its inputs
include: "deploy.smk"

rule test_target:
    """
    For testing purposes such as CI workflows.
    """
    input: "auspice/avian-flu_h5n1_ha_all-time.json"

class InvalidConfigError(Exception):
    pass

# This uses the `InvalidConfigError` defined above
include: "merge_inputs.smk"

rule filter_sequences_by_subtype:
    input:
        sequences = input_sequences,
        metadata = input_metadata,
    output:
        sequences = "results/{subtype}/{segment}/sequences.fasta",
    params:
        # We don't have all wildcards set here (too early!) so we need to manually specify them for `resolve_config_value`
        # (Note that we do have w.segment set, but we deliberately don't use it as the query must not vary by segment)
        subtypes = lambda w: resolve_config_value('subtype_query')({'subtype': w.subtype, 'segment': '*', 'time': '*'})
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --query {params.subtypes!r} \
            --output-sequences {output.sequences}
        """

rule filter_metadata_by_subtype:
    input:
        metadata = input_metadata,
    output:
        metadata = "results/{subtype}/metadata.tsv",
    params:
        # We don't have all wildcards set here (too early!) so we need to manually specify them for `resolve_config_value`
        subtypes = lambda w: resolve_config_value('subtype_query')({'subtype': w.subtype, 'segment': '*', 'time': '*'})
    shell:
        """
        augur filter \
            --metadata {input.metadata} \
            --query {params.subtypes!r} \
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

def refine_clock_rates(w):
    if w.segment == 'genome':
        # TODO - can we shift this logic into `rules/genome.smk` ?
        # calculate the genome rate via a weighted average of the segment rates

        # Sanity check: check that the config hasn't (mistakenly?) set a genome rate
        # Note: this relies on the genome builds not setting a catch-all "*/*/*" config
        try:
            info = resolve_config_value('refine', 'clock_rates')(w)
        except InvalidConfigError as e:
            pass # this is behaviour we expect - i.e. no clock rates should be provided for genome
        else:
            raise InvalidConfigError("config['refine']['clock_rates'] should not contain a key which can be used for 'genome' builds as this is calculated by the workflow")

        segment_lengths = [
            resolve_config_value('refine', 'segment_lengths')(custom_wildcards)
            for custom_wildcards in
            [{"subtype": w.subtype, "time": w.time, "segment": segment} for segment in SEGMENTS]
        ]
        clock_rates = [
            resolve_config_value('refine', 'clock_rates')(custom_wildcards)
            for custom_wildcards in
            [{"subtype": w.subtype, "time": w.time, "segment": segment} for segment in SEGMENTS]
        ]
        for rates in clock_rates:
            assert isinstance(rates, list) and len(rates)==2, "The clock rates for each segment must be a list of (rate, std-dev), not {rates!r}"

        mean = sum([length*rate[0] for length,rate in zip(segment_lengths, clock_rates)]) / sum(segment_lengths)
        stdev = mean/2
        return f"--clock-rate {mean} --clock-std-dev {stdev}"

    info = resolve_config_value('refine', 'clock_rates')(w)
    if info == "":
        return ""

    assert isinstance(info, list) and len(info)==2, "The clock rates for {w.subtype!r}/{w.time!r}/{w.segment!r} must be a list of (rate, std-dev), not {info!r}"
    return f"--clock-rate {info[0]} --clock-std-dev {info[1]}"

def refine_clock_filter(w):
    filter = resolve_config_value('refine', 'clock_filter_iqd')(w)
    return f"--clock-filter-iqd {filter}" if filter else ""


rule add_h5_clade:
    message: "Adding in a column for h5 clade numbering"
    input:
        metadata = "results/{subtype}/metadata.tsv",
        clades_file = resolve_config_fields_path('clades_file')
    output:
        metadata= "results/{subtype}/metadata-with-clade.tsv"
    params:
        script = script("clade-labeling/add-clades.py")
    shell:
        r"""
        python {params.script} \
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

    exclude_where = resolve_config_value('filter', 'exclude_where')(wildcards)
    # If SAME_STRAINS (and due to the above conditional we have the HA segment at this point)
    # then we want to restrict to strains present in all 8 segments. Note that force-included
    # strains may not have all segments, but that's preferable to filtering them out.
    if SAME_STRAINS:
        exclude_where += f" n_segments!={len(SEGMENTS)}"


    cmd = ""

    group_by_value = resolve_config_value('filter', 'group_by')(wildcards)
    cmd += f" --group-by {group_by_value}" if group_by_value else ""

    cmd += f" --subsample-max-sequences {resolve_config_value('filter', 'target_sequences_per_tree')(wildcards)}"
    cmd += f" --min-date {resolve_config_value('filter', 'min_date')(wildcards)}"
    cmd += f" --include {input.include}"
    cmd += f" --exclude-where {exclude_where}"
    cmd += f" --min-length {resolve_config_value('filter', 'min_length')(wildcards)}"
    cmd += f" --non-nucleotide"
    return cmd

rule filter:
    input:
        sequences = "results/{subtype}/{segment}/sequences.fasta",
        metadata = metadata_by_wildcards,
        exclude = resolve_config_fields_path('dropped_strains'),
        include = resolve_config_fields_path('include_strains'),
        strains = lambda w: f"results/{w.subtype}/ha/{w.time}/filtered.txt" if (SAME_STRAINS and w.segment!='ha') else [],
    output:
        sequences = "results/{subtype}/{segment}/{time}/filtered.fasta",
        strains = "results/{subtype}/{segment}/{time}/filtered.txt",
        metadata = "results/{subtype}/{segment}/{time}/metadata.tsv",
    params:
        args = _filter_params,
    wildcard_constraints:
        # The genome build has a different approach to filtering (see genome.smk)
        segment="(?!genome)[^_/]+"
    shell:
        r"""
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
        reference = resolve_config_fields_path('reference'),
    output:
        alignment = "results/{subtype}/{segment}/{time}/aligned.fasta"
    wildcard_constraints:
        # for genome builds we don't use this rule; see `rule join_segments`
        segment = "|".join(seg for seg in SEGMENTS)
    threads:
        4
    shell:
        r"""
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
    """
    Returns the config-defined rooting parameter (for `augur refine`) together with the argument,
    i.e. "--root <parameter>". If none is defined we return the empty string.
    """
    root = resolve_config_value('refine', 'root')(wildcards)
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
        coalescent = resolve_config_value('refine', 'coalescent'),
        date_inference = resolve_config_value('refine', 'date_inference'),
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
    The genome biulds introduces an additional step beyond `augur refine`, which is
    why this function exists.
    """
    if w.subtype=='h5n1-cattle-outbreak' and w.segment!='genome':
        return "results/{subtype}/{segment}/{time}/tree_outbreak-clade.nwk"
    return "results/{subtype}/{segment}/{time}/tree.nwk",

def ancestral_root_seq(wildcards):
    # The root-seq(uence) is a file, not a name
    value = resolve_config_fields_path('ancestral', 'root_seq')(wildcards)
    if not value: # falsey values result in an empty string path, i.e. we skip the --root-sequence argument
        return ""
    return f"--root-sequence {value}"

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = refined_tree,
        alignment = rules.align.output
    output:
        node_data = "results/{subtype}/{segment}/{time}/nt-muts.json"
    params:
        inference = resolve_config_value('ancestral', 'inference'),
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
        reference = resolve_config_fields_path('reference')
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
    args = ""

    if columns:=resolve_config_value('traits', 'columns')(wildcards):
        args += f"--columns {columns}"
    else:
        raise InvalidConfigError(f"You must define columns for 'augur traits' to run on via config['traits']['columns']")

    if bias:=resolve_config_value('traits', 'sampling_bias_correction')(wildcards):
        args += f" --sampling-bias-correction {bias}"

    if confidence:=resolve_config_value('traits', 'confidence')(wildcards):
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
    params:
        script = script("annotate-ha-cleavage-site.py")
    shell:
        r"""
        python {params.script} \
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


def additional_args_export(wildcards):
    args = []

    if title:=resolve_config_value('export', 'title')(wildcards):
        # Config defined title strings may include wildcards to be expanded
        title = title.format(segment=wildcards.segment.upper(), subtype=wildcards.subtype.upper(), time=wildcards.time)
        args += ["--title", title]

    return args

rule auspice_config:
    """
    Depending on the build we may make wildcard-dependent modifications to the auspice config JSON.
    If we implement config overlays in augur this rule will become unnecessary.
    """
    input:
        auspice_config = resolve_config_fields_path('auspice_config'),
    output:
        auspice_config = "results/{subtype}/{segment}/{time}/auspice-config.json",
    run:
        import json
        with open(input.auspice_config) as fh:
            auspice_config = json.load(fh)
        if wildcards.subtype in ["h5n1-cattle-outbreak", "h5n1-d1.1"]:
            if wildcards.segment == "genome":
                auspice_config['display_defaults']['distance_measure'] = "num_date" if wildcards.subtype == "h5n1-cattle-outbreak" else "div"
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

        # If we have a coloring for 'genoflu' and we're not a genome build then export the genoflu call for the segment
        genoflu_idx = next((i for i,c in enumerate(auspice_config['colorings']) if c['key']=='genoflu'), None)
        if wildcards.segment!='genome' and genoflu_idx is not None:
            auspice_config['colorings'].insert(genoflu_idx+1, {
                    "key": f"genoflu_{wildcards.segment.upper()}",
                    "title": f"GenoFLU ({wildcards.segment.upper()} segment)",
                    "type": "categorical",
                })
        with open(output.auspice_config, 'w') as fh:
            json.dump(auspice_config, fh, indent=2)

rule colors:
    input:
        # TODO - revisit the metadata input future work. The TSV used here determines the colours used, so if we use the same metadata
        # (and it's a superset) for a set of builds (e.g. genome + 8 segments) then the colours are consistent.
        metadata = lambda w: "results/{subtype}/genome/{time}/metadata.tsv" \
            if w.subtype in ['h5n1-cattle-outbreak', 'h5n1-d1.1'] \
            else rules.filter.output.metadata,
        colors = resolve_config_fields_path('colors', 'hardcoded'),
        ordering = resolve_config_fields_path('colors', 'ordering'),
        schemes = resolve_config_fields_path('colors', 'schemes'),
    output:
        colors = "results/{subtype}/{segment}/{time}/colors.tsv",
    params:
        duplications = lambda w: ["=".join(pair) for pair in resolve_config_value('colors', 'duplications')(w)],
        script = script("assign-colors.py"),
    shell:
        r"""
        cp {input.colors} {output.colors} && \
        python3 {params.script} \
            --metadata {input.metadata} \
            --ordering {input.ordering} \
            --color-schemes {input.schemes} \
            --duplications {params.duplications} \
        >> {output.colors}
        """

rule export:
    """
    Export the files into results/ and then use a subsequent rule to move these to the
    auspice/ directory
    """
    input:
        tree = refined_tree,
        metadata = rules.filter.output.metadata,
        node_data = export_node_data_files,
        colors = "results/{subtype}/{segment}/{time}/colors.tsv",
        lat_longs = resolve_config_fields_path("lat_longs"),
        auspice_config = rules.auspice_config.output.auspice_config,
        description = resolve_config_fields_path("description"),
    output:
        auspice_json = "results/{subtype}/{segment}/{time}/auspice-dataset.json"
    params:
        additional_args = additional_args_export
    shell:
        r"""
        export AUGUR_RECURSION_LIMIT=20000
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --description {input.description} \
            --include-root-sequence-inline \
            {params.additional_args:q} \
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
        assert subtype=='h5n1-cattle-outbreak' or subtype=='h5n1-d1.1', \
            "Only h5n1 builds produce an Auspice dataset without a time component in the filename"
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

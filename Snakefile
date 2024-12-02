# constrain the wildcards to not include `_` which we use to separate "parts" of filenames (where a part may be a wildcard itself)
wildcard_constraints:
    subtype = "[^_/]+",
    segment = "[^_/]+",
    time = "[^_/]+",

# defined before extra rules `include`d as they reference this constant
SEGMENTS = ["pb2", "pb1", "pa", "ha","np", "na", "mp", "ns"]
#SUBTYPES = ["h5n1", "h5nx", "h7n9", "h9n2"]

CURRENT_BASEDIR = workflow.current_basedir # TODO XXX store this value here - can't access within functions because workflow.included_stack is empty

# Load the base config.yaml relative to the entry snakefile (i.e. not this snakefile)
if os.path.exists(os.path.join(workflow.basedir, 'config.yaml')):
    configfile: os.path.join(workflow.basedir, 'config.yaml')

# load a config.yaml file if it exists in the current working directory
if os.path.exists("config.yaml"):
    configfile: "config.yaml"

from pprint import pp; pp(config, stream=sys.stderr) # TODO XXX remove

class InvalidConfigError(Exception):
    pass

def resolve_config_path(path):
    """
    Resolve a relative *path* given in a configuration value. Returns a
    function which takes a single argument *wildcards* and returns the resolved
    path with any '{x}' substrings are replaced by their corresponding wildcards
    filled in

    Search order (first match returned):
    1. Relative to the analysis directory
    2. Relative to the directory the entry snakefile was in. Typically this is
       not the Snakefile you are looking at now but (e.g.) the one in
       avian-flu/gisaid
    3. Relative to where this Snakefile is (i.e. `avian-flu/`)
    """
    if not isinstance(path, str):
        raise InvalidConfigError(f"Config path provided to resolve_config_path must be a string. Provided value: {str(path)}")

    def resolve(wildcards):
        try:
            path_expanded = expand(path, **wildcards)[0]
        except snakemake.exceptions.WildcardError as e:
            # str(e) looks like "No values given for wildcard 'subtypes'."
            raise InvalidConfigError(f"resolve_config_path called with path {path!r} however {str(e)}")

        # check if the path exists relative to the working analysis directory
        if os.path.exists(path_expanded):
            return path_expanded

        # Check if the path exists relative to the subdir where the entry snakefile is
        # (e.g. avian-flu/gisaid). If you want to use further subdirectories (e.g. avian-flu/gisaid/config/x.tsv)
        # you're expected to supply the 'config/x.tsv' as the value in the config YAML
        # NOTE: this means analysis directory overrides have to use that same 'config/x.tsv' structure, but
        # given the different directories avian-flu uses that's acceptable. In other words, if we standardised
        # avian-flu then we could add subdirectories to the search order here
        basepath = os.path.join(workflow.basedir, path_expanded)
        if os.path.exists(basepath):
            return basepath

        # Check if the path exists relative to where _this snakefile is_
        # NOTE: We may move this Snakefile to `rules/base.smk`, or perhaps go further and split
        # these functions out into `rules/common.smk` etc. If we do so then then I think the correct behvaiour for
        # step (3) is to search for paths relative to `avian-flu` (not `avian-flu/rules`)
        if workflow.basedir != CURRENT_BASEDIR:
            basepath = os.path.join(CURRENT_BASEDIR, path_expanded)
            if os.path.exists(basepath):
                return basepath

        raise InvalidConfigError(f"Unable to resolve the config-provided path {path!r}, expanded to {path_expanded!r} after filling in wildcards. "
            f"The following directories were searched:\n"
            f"\t1. {os.path.abspath(os.curdir)} (current working directory)\n"
            f"\t2. {workflow.basedir} (where the entry snakefile is)\n"
            f"\t3. {CURRENT_BASEDIR} (where the main avian-flu snakefile is)\n")
    return resolve

def resolve_config_value(rule_parts, wildcards, sep="/"):
    """
    Resolve a config value defined by the *rule_parts of the config,
    e.g. rule_parts = ['filter', 'min_length'] then we expect a scalar or
    a dictionary to be present at config['filter']['min_length'].

    If a scalar then that value is returned, i.e. it's always the same no matter
    what the wildcards are.

    If a dictionary we search it for the relevant value by finding the closest matching
    key once wildcards have been considered. For instance in a three-tiered wildcard pipeline
    such as this with subtype, segment and time we interpret these as ordered in specificity.
    For instance, if only one wildcard value is specified in the config (the others are '*')
    then matching subtype is more specific than segment. Given example 
    wildcard values of {subtype=h5nx, segment=pb2, time=2y} then we have a search order of:
    - 'h5nx/pb2/2y'   ─ all 3 wildcard values specified
    - 'h5nx/pb2/*'    ┐
    - 'h5nx/*/2y'     ├ 2/3 wildcard values specified
    - '*/pb2/2y'      ┘
    - 'h5nx/*/*'      ┐
    - '*/pb2/*'       ├ 1/3 wildcard values specified
    - '*/*/2y'        ┘
    - '*/*/*'         ─ default / fall-back
    and the first key present in the config is used.
    """
    try:
        config_lookup = config
        for i,rule_key in enumerate(rule_parts): # or use functools.reduce etc
            config_lookup = config_lookup[rule_key]
    except KeyError:
        raise InvalidConfigError('Config missing entire entry for config'+''.join(['["'+rule_parts[j]+'""]' for j in range(0,i+1)]))

    if any([isinstance(config_lookup, t) for t in [float, int, bool, str]]):
        return config_lookup

    if not isinstance(config_lookup, dict):
        raise InvalidConfigError(f"ERROR: config under {'.'.join(rule_parts)} must be a scalar value or a dictionary")

    wild_keys = ['subtype', 'segment', 'time'] # workflow specific
    search_keys = [                            # workflow independent, as long as there are 3 or fewer wildcard categories
        sep.join([wildcards[k] for k in wild_keys]),
        *([sep.join(['*' if i==k else wildcards[key] for k,key in enumerate(wild_keys)])
            for i in range(len(wild_keys)-1, -1, -1)] if len(wild_keys)>=2 else []),
        *([sep.join(['*' if i!=k else wildcards[key] for k,key in enumerate(wild_keys)])
            for i in range(0, len(wild_keys))] if len(wild_keys)==3 else []),
        sep.join(['*']*len(wild_keys))
    ]

    for key in search_keys:
        if key in config_lookup:
            return config_lookup[key]
    msg  =  'Config structure incorrect or incomplete for config'+''.join(['["'+rule_parts[j]+'"]' for j in range(0,i+1)])
    msg += f'\n\tThe dictionary is missing a matching key for the current target of {search_keys[0]!r}, or a fuzzy match (i.e. using "*" placeholders)'
    msg +=  '\n\tP.S. If you want to use a single value across all builds then set a scalar value (number, string, boolean)'
    raise InvalidConfigError(msg)



# The config option `same_strains_per_segment=True'` (e.g. supplied to snakemake via --config command line argument)
# will change the behaviour of the workflow to use the same strains for each segment. This is achieved via these steps:
# (1) Filter the HA segment as normal plus filter to those strains with 8 segments
# (2) Filter the other segments by simply force-including the same strains as (1)
SAME_STRAINS = bool(config.get('same_strains_per_segment', False))

NEXTSTRAIN_PUBLIC_BUCKET = "s3://nextstrain-data/"
S3_SRC = config.get('s3_src', {})
LOCAL_INGEST = config.get('local_ingest', None)

def sanity_check_config():
    if not len(config.keys()):
        print("-"*80 + "\nNo config loaded!", file=sys.stderr)
        print("Avian-flu is indented to be run from the snakefile inside a subdir " 
            "(e.g. gisaid/Snakefile) which will pick up the default configfile for that workflow. " 
            "Alternatively you can pass in the config via `--configfile`", file=sys.stderr)
        print("-"*80, file=sys.stderr)
        raise InvalidConfigError("No config")

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
        #sequences = expand("results/{subtype}/{segment}/sequences.fasta", segment=SEGMENTS, subtype=SUBTYPES),
        #metadata = expand("results/{subtype}/metadata.tsv", segment=SEGMENTS, subtype=SUBTYPES)


# This must be after the `all` rule above since it depends on its inputs
# Note this is relative to the `workflow.current_basedir`
include: "rules/deploy.smk"

rule test_target:
    """
    For testing purposes such as CI workflows.
    """
    input: "auspice/avian-flu_h5n1_ha_all-time.json"

# TODO - I find this indirection more confusing than helpful and I'd rather just
# specify `config['colors']` in the rules which use it (e.g.)
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
    # TODO - shift this to config
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

def refine_clock_rates(w):

    if w.segment == 'genome':
        # TODO XXX - can this be part of cattle-flu.smk?
        # calculate the genome rate via a weighted average of the segment rates

        # Sanity check: check that the config hasn't (mistakenly?) set a genome rate
        try:
            info = resolve_config_value(['refine', 'clock_rates'], w)
        except InvalidConfigError as e:
            pass # this is behaviour we expect - i.e. no clock rates should be provided for genome
        else:
            raise InvalidConfigError("config['refine']['clock_rates'] should not contain a key which can be used for 'genome' builds as this is calculated by the workflow")

        segment_lengths = [
            resolve_config_value(['refine', 'segment_lengths'], wildcards)
            for wildcards in
            [{"subtype": w.subtype, "time": w.time, "segment": segment} for segment in SEGMENTS]
        ]
        clock_rates = [
            resolve_config_value(['refine', 'clock_rates'], wildcards)
            for wildcards in
            [{"subtype": w.subtype, "time": w.time, "segment": segment} for segment in SEGMENTS]
        ]
        for rates in clock_rates:
            assert isinstance(rates, list) and len(rates)==2, "The clock rates for each segment must be a list of (rate, std-dev), not {rates!r}"

        mean = sum([length*rate[0] for length,rate in zip(segment_lengths, clock_rates)]) / sum(segment_lengths)
        stdev = mean/2
        return f"--clock-rate {mean} --clock-std-dev {stdev}"

    info = resolve_config_value(['refine', 'clock_rates'], w)
    if info == "":
        return ""

    assert isinstance(info, list) and len(info)==2, "The clock rates for {w.subtype!r}/{w.time!r}/{w.segment!r} must be a list of (rate, std-dev), not {info!r}"
    return f"--clock-rate {info[0]} --clock-std-dev {info[1]}"

def refine_clock_filter(w):
    filter = resolve_config_value(['refine', 'clock_filter_iqd'], w)
    return f"--clock-filter-iqd {filter}" if filter else ""


rule add_h5_clade:
    message: "Adding in a column for h5 clade numbering"
    input:
        metadata = "results/{subtype}/metadata.tsv",
        clades_file = resolve_config_path(files.clades_file)
    output:
        metadata= "results/{subtype}/metadata-with-clade.tsv"
    params:
        script = os.path.join(workflow.current_basedir, "clade-labeling/add-clades.py")
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

    exclude_where = resolve_config_value(['filter', 'exclude_where'], wildcards)
    # If SAME_STRAINS (and due to the above conditional we have the HA segment at this point)
    # then we want to restrict to strains present in all 8 segments. Note that force-included
    # strains may not have all segments, but that's preferable to filtering them out.
    if SAME_STRAINS:
        exclude_where += f" n_segments!={len(SEGMENTS)}"


    cmd = ""

    group_by_value = resolve_config_value(['filter', 'group_by'], wildcards)
    cmd += f" --group-by {group_by_value}" if group_by_value else ""

    cmd += f" --subsample-max-sequences {config['target_sequences_per_tree']}"
    cmd += f" --min-date {resolve_config_value(['filter', 'min_date'], wildcards)}"
    cmd += f" --include {input.include}"
    cmd += f" --exclude-where {exclude_where}"
    cmd += f" --min-length {resolve_config_value(['filter', 'min_length'], wildcards)}"
    cmd += f" --non-nucleotide"
    return cmd

rule filter:
    input:
        sequences = "results/{subtype}/{segment}/sequences.fasta",
        metadata = metadata_by_wildcards,
        exclude = resolve_config_path(files.dropped_strains),
        include = resolve_config_path(files.include_strains),
        strains = lambda w: f"results/{w.subtype}/ha/{w.time}/filtered.txt" if (SAME_STRAINS and w.segment!='ha') else [],
    output:
        sequences = "results/{subtype}/{segment}/{time}/filtered.fasta",
        strains = "results/{subtype}/{segment}/{time}/filtered.txt",
        metadata = "results/{subtype}/{segment}/{time}/metadata.tsv",
    params:
        args = _filter_params,
    wildcard_constraints:
        # The genome build has a different approach to filtering (see cattle-flu.smk)
        segment="(?!genome)[^_/]+"
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
        reference = resolve_config_path(files.reference)
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
        root = lambda w: f"--root {resolve_config_value(['refine', 'root'], w)}" if resolve_config_value(['refine', 'root'], w) else ''
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
    config_value = resolve_config_value(['ancestral', 'root_seq'], wildcards)
    if not config_value: # falsey values skip the --root-sequence argument
        return ""
    return f"--root-sequence {resolve_config_path(config_value)(wildcards)}"

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
        reference = lambda w: resolve_config_path(config['genome_reference'] if w.segment=='genome' else files.reference)(w)
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
    args = f"--columns {resolve_config_value(['traits', 'columns'], wildcards)}"
    if bias:=resolve_config_value(['traits', 'sampling_bias_correction'], wildcards):
        args += f" --sampling-bias-correction {bias}"
    if confidence:=resolve_config_value(['traits', 'confidence'], wildcards):
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
        r"""
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
        script = os.path.join(workflow.current_basedir, "scripts/annotate-ha-cleavage-site.py")
    shell:
        """
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


def additional_export_config(wildcards):
    args = ""

    title = resolve_config_value(['export', 'title'], wildcards)
    if title:
        # Config defined title strings may include wildcards to be expanded
        title = title.format(segment=wildcards.segment.upper(), subtype=wildcards.subtype.upper(), time=wildcards.time)
        args += f"--title {title!r}"
        
    return args

rule auspice_config:
    """
    Depending on the build we may make wildcard-dependent modifications to the auspice config JSON.
    If we implement config overlays in augur this rule will become unnecessary.
    """
    input:
        auspice_config = resolve_config_path(files.auspice_config),
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

rule colors:
    input:
        colors = resolve_config_path(files.colors),
    output:
        colors = "results/{subtype}/{segment}/{time}/colors.tsv",
    shell:
        """
        cp {input.colors} {output.colors}
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
        lat_longs = resolve_config_path(files.lat_longs),
        auspice_config = rules.auspice_config.output.auspice_config,
        description = resolve_config_path(files.description),
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

# Add in additional rules at/towards the end of the Snakefile so that
# they have access to functions and values declared above
for rule_file in config.get('custom_rules', []):
    # Relative custom rule paths in the config are expected to be relative to the analysis (working) directory
    include: os.path.join(os.getcwd(), rule_file)

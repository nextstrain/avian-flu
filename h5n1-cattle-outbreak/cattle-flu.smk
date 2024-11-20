# This rule file is conditionally included for the h5n1-cattle-outbreak build

rule filter_segments_for_genome:
    # Note for developers: The {genome_seg} wildcard here is not the {segment}
    # wildcard that is used throughtout the pipeline. This rule is invoked when
    # we are running a genome build (so {segment}=genome for most rules) and the
    # `join_segments` rule ultimately requests filtered sequences / alignments
    # for each constituent segment. We call this {genome_seg} just to
    # distinguish it when reading the code. 
    input:
        sequences = "results/{subtype}/{genome_seg}/sequences.fasta",
        metadata = "results/{subtype}/metadata-with-clade.tsv", # TODO: use a function here instead of hardcoding
        include = resolve_config_path(config['include_strains']),
        exclude = resolve_config_path(config['dropped_strains']),
    output:
        sequences = "results/{subtype}/{segment}/{time}/filtered_{genome_seg}.fasta"
    params:
        min_date = "2024-01-01",
        query = 'region == "North America"'
    wildcard_constraints:
        subtype = 'h5n1-cattle-outbreak',
        segment = 'genome',
        time = 'default',
    log: "logs/{subtype}/{segment}/{time}/filtered_{genome_seg}.txt",
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --include {input.include} \
            --exclude {input.exclude} \
            --min-date {params.min_date} \
            --query {params.query:q} \
            --output-log {log} \
            --output-sequences {output.sequences}
        """

rule align_segments_for_genome:
    input:
        sequences = "results/{subtype}/{segment}/{time}/filtered_{genome_seg}.fasta",
        # Use the H5N1 reference sequences for alignment
        reference = lambda w: [
            resolve_config_path(expanded)(w)
            for expanded in
            expand(config['reference'], subtype='h5n1', segment=w.genome_seg)
        ]
    output:
        alignment = "results/{subtype}/{segment}/{time}/aligned_{genome_seg}.fasta"
    wildcard_constraints:
        subtype = 'h5n1-cattle-outbreak',
        segment = 'genome',
        time = 'default',
    threads:
        8
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --remove-reference \
            --fill-gaps \
            --nthreads {threads}
        """

rule join_segments:
    # the output of this rule is the same as 'rule align' but the wildcard constraints
    # allow snakemake to choose the correct rule to run. Note that `wildcards.segment="genome"`
    # here, and for that we need alignments for 8 individual segments, which we refer to as `wildcards.genome_seg`
    input:
        alignment = expand("results/{{subtype}}/{{segment}}/{{time}}/aligned_{genome_seg}.fasta", genome_seg=SEGMENTS) 
    output:
        alignment = "results/{subtype}/{segment}/{time}/aligned.fasta"
    wildcard_constraints:
        subtype = 'h5n1-cattle-outbreak',
        segment = 'genome',
        time = 'default',
    params:
        script = os.path.join(workflow.current_basedir, "../scripts/join-segments.py")
    shell:
        """
        python {params.script} \
            --segments {input.alignment} \
            --output {output.alignment}
        """

rule genome_metadata:
    input:
        sequences = "results/{subtype}/{segment}/{time}/aligned.fasta",
        metadata = metadata_by_wildcards,
    output:
        metadata = temp("results/{subtype}/{segment}/{time}/metadata_intermediate.tsv")
    wildcard_constraints:
        subtype = 'h5n1-cattle-outbreak',
        segment = 'genome',
        time = 'default',
    shell:
        """
        augur filter --metadata {input.metadata} --sequences {input.sequences} --output-metadata {output.metadata}
        """

def metadata_columns_to_duplicate(wildcards):
    column = resolve_config_value(['traits', 'columns'], wildcards)
    assert isinstance(column, str) and ' ' not in column, \
        "cattle-flu.smk::metadata_columns_to_duplicate: Genome workflow only expects there to be a single column to run `augur traits` on."
    return [column, column+"_metadata"]

rule add_metadata_columns_to_show_non_inferred_values:
    """
    Genome builds run `augur traits` for "division" (we assert this below) so we want to add a metadata
    column `division_metadata` which is a duplicate of `division`.

    NOTE: long-term we should be consulting `traits_params()` to work out the columns to duplicate, but
    that function's not visible to this .smk file so would require deeper refactoring.
    """
    input:
        metadata = "results/{subtype}/{segment}/{time}/metadata_intermediate.tsv"
    output:
        metadata = "results/{subtype}/{segment}/{time}/metadata.tsv"
    wildcard_constraints:
        subtype="h5n1-cattle-outbreak",
        segment="genome",
        time="default",
    params:
        columns = lambda w: metadata_columns_to_duplicate(w),
    shell:
        """
        cat {input.metadata} | csvtk mutate -t -f {params.columns[0]} -n {params.columns[1]} > {output.metadata}
        """

ruleorder: add_metadata_columns_to_show_non_inferred_values > filter

rule prune_tree:
    input:
        tree = "results/{subtype}/{segment}/{time}/tree.nwk",
        strains = "auspice/avian-flu_h5n1-cattle-outbreak_genome.json",
    output:
        tree = "results/{subtype}/{segment}/{time}/tree_outbreak-clade.nwk",
        node_data = "results/{subtype}/{segment}/{time}/outbreak-clade-strains-in-genome-tree.json",
    wildcard_constraints:
        subtype="h5n1-cattle-outbreak",
        time="default",
    params:
        script = os.path.join(workflow.current_basedir, "../scripts/restrict-via-common-ancestor.py")
    shell:
        r"""
        python3 {params.script} \
            --tree {input.tree} \
            --strains {input.strains} \
            --output-tree {output.tree} \
            --output-metadata {output.node_data}
        """

rule colors_genome:
    # TODO: add these input files / params to the config YAML. The config YAML must also
    # define the concept of whether this rule should run so this isn't trivial and is
    # thus left as a to-do. Once they are in the YAML we should switch to `resolve_config_path`
    input:
        metadata = "results/{subtype}/genome/{time}/metadata.tsv", # Always use the genome metadata, even for segment builds
        ordering = lambda w: resolve_config_path("config/h5n1-cattle-outbreak/color_ordering.tsv", w),
        schemes = lambda w: resolve_config_path("config/h5n1-cattle-outbreak/color_schemes.tsv", w),
        colors = lambda w: resolve_config_path(files.colors, w)
    output:
        colors = "results/{subtype}/{segment}/{time}/colors.tsv",
    params:
        duplications = lambda w: "=".join(metadata_columns_to_duplicate({**w, 'segment': 'genome'})),
        script = os.path.join(workflow.current_basedir, "../scripts/assign-colors.py")
    wildcard_constraints:
        subtype="h5n1-cattle-outbreak",
        time="default",
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

ruleorder: colors_genome > colors

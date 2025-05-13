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
        include = resolve_config_path('include_strains'),
        exclude = resolve_config_path('dropped_strains'),
    output:
        sequences = "results/{subtype}/{segment}/{time}/filtered_{genome_seg}.fasta"
    wildcard_constraints:
        subtype = 'h5n1-cattle-outbreak|h5n1-d1.1',
        segment = 'genome',
        time = 'default',
    log: "logs/{subtype}/{segment}/{time}/filtered_{genome_seg}.txt",
    shell:
        r"""
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --include {input.include} \
            --exclude {input.exclude} \
            --output-log {log} \
            --output-sequences {output.sequences}
        """

rule align_segments_for_genome:
    input:
        sequences = "results/{subtype}/{segment}/{time}/filtered_{genome_seg}.fasta",
        # Use the H5N1 reference sequences for alignment
        reference = lambda w: resolve_config_path('reference')({'subtype': w.subtype, 'segment': w.genome_seg, 'time': 'default'})
    output:
        alignment = "results/{subtype}/{segment}/{time}/aligned_{genome_seg}.fasta"
    wildcard_constraints:
        subtype = 'h5n1-cattle-outbreak|h5n1-d1.1',
        segment = 'genome',
        time = 'default',
    threads:
        8
    shell:
        r"""
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
        alignment = "results/{subtype}/{segment}/{time}/aligned-unmasked.fasta"
    wildcard_constraints:
        subtype = 'h5n1-cattle-outbreak|h5n1-d1.1',
        segment = 'genome',
        time = 'default',
    params:
        script = script("join-segments.py"),
    shell:
        r"""
        python {params.script} \
            --segments {input.alignment} \
            --output {output.alignment}
        """

rule mask_genome:
    input:
        alignment = "results/{subtype}/{segment}/{time}/aligned-unmasked.fasta"
    output:
        alignment = "results/{subtype}/{segment}/{time}/aligned.fasta",
    params:
        percentage = resolve_config_value('mask', 'min_support'),
        script = script("mask.py"),
    wildcard_constraints:
        subtype = 'h5n1-cattle-outbreak|h5n1-d1.1',
        segment = 'genome',
        time = 'default',
    shell:
        r"""
        python {params.script} \
            --alignment {input.alignment} \
            --percentage {params.percentage} \
            --output {output.alignment}
        """

rule genome_metadata:
    input:
        sequences = "results/{subtype}/{segment}/{time}/aligned.fasta",
        metadata = metadata_by_wildcards,
    output:
        metadata = temp("results/{subtype}/{segment}/{time}/metadata_intermediate.tsv")
    wildcard_constraints:
        subtype = 'h5n1-cattle-outbreak|h5n1-d1.1',
        segment = 'genome',
        time = 'default',
    shell:
        """
        augur filter --metadata {input.metadata} --sequences {input.sequences} --output-metadata {output.metadata}
        """

def _assert_traits_column(wildcards):
    columns = resolve_config_value('traits', 'columns')(wildcards)
    if columns != 'division':
        raise InvalidConfigError("The genome-specific pipeline currently requires the (genome-specific) build to infer " \
            f"'augur traits' for 'division' only, not {columns!r}.")

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
        subtype='h5n1-cattle-outbreak|h5n1-d1.1',
        segment="genome",
        time="default",
    params:
        old_column = "division",
        new_column = "division_metadata",
        assert_traits = _assert_traits_column,
    shell:
        """
        cat {input.metadata} | csvtk mutate -t -f {params.old_column} -n {params.new_column} > {output.metadata}
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
        subtype='h5n1-cattle-outbreak|h5n1-d1.1',
        time="default",
    params:
        script = script("restrict-via-common-ancestor.py")
    shell:
        """
        python3 {params.script} \
            --tree {input.tree} \
            --strains {input.strains} \
            --output-tree {output.tree} \
            --output-metadata {output.node_data}
        """

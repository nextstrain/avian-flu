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
        metadata = "results/{subtype}/metadata-with-clade.tsv",
        include = config['include_strains'],
        exclude = config['dropped_strains'],
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
        reference = lambda w: expand(config['reference'], subtype='h5n1', segment=w.genome_seg)
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
    shell:
        """
        python scripts/join-segments.py \
            --segments {input.alignment} \
            --output {output.alignment}
        """

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
    shell:
        """
        python3 scripts/restrict-via-common-ancestor.py \
            --tree {input.tree} \
            --strains {input.strains} \
            --output-tree {output.tree} \
            --output-metadata {output.node_data}
        """

def assert_expected_config(w):
    try:
        # TODO: once we refactor things we should use `get_config()` here
        # see <https://github.com/nextstrain/avian-flu/pull/100#discussion_r1823047047>
        # but currently this snakefile doesn't have access to that function.
        assert len(config['traits']['genome_columns'])==1 and config['traits']['genome_columns']['FALLBACK']=="division"
    except Exception as err:
        raise Exception("Rule add_metadata_columns_to_show_non_inferred_values expected a certain format for config['traits'] that has since changed") from err

rule add_metadata_columns_to_show_non_inferred_values:
    """
    Genome builds run `augur traits` for "division" (we assert this below) so we want to add a metadata
    column `division_metadata` which is a duplicate of `division`.

    NOTE: long-term we should be consulting `traits_params()` to work out the columns to duplicate, but
    that function's not visible to this .smk file so would require deeper refactoring.
    """
    input:
        metadata = "results/{subtype}/metadata-with-clade.tsv",
    output:
        metadata = "results/{subtype}/{segment}/{time}/metadata-with-clade-and-non-inferred-values.tsv",
    wildcard_constraints:
        subtype="h5n1-cattle-outbreak",
        segment="genome",
        time="default",
    params:
        old_column = "division",
        new_column = "division_metadata",
        assert_traits = assert_expected_config,
    shell:
        """
        cat {input.metadata} | csvtk mutate -t -f {params.old_column} -n {params.new_column} > {output.metadata}
        """

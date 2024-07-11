
rule h5n1_cattle_outbreak:
    """helper rule to avoid typing out all targets on the command line"""
    input:
        expand("auspice/avian-flu_h5n1-cattle-outbreak_{segment}.json", segment=[*SEGMENTS, 'genome'])

rule filter_segments_for_genome:
    # Note for developers: The {genome_seg} wildcard here is not the {segment}
    # wildcard that is used throughtout the pipeline. This rule is invoked when
    # we are running a genome build (so {segment}=genome for most rules) and the
    # `join_segments` rule ultimately requests filtered sequences / alignments
    # for each constituent segment. We call this {genome_seg} just to
    # distinguish it when reading the code. 
    input:
        sequences = "data/sequences_{subtype}_{genome_seg}.fasta",
        metadata = "results/metadata-with-clade_{subtype}.tsv",
        # Following will be part of the config YAML TODO XXX
        exclude = "config/dropped_strains_{subtype}.txt",
        include = "config/include_strains_{subtype}.txt",
    output:
        sequences = "results/filtered-for-genome_{subtype}_{genome_seg}_{time}.fasta",
    params:
        min_date = "2024-01-01",
        query = 'region == "North America"'
    log: "logs/filtered-for-genome_{subtype}_{genome_seg}_{time}.txt",
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
        sequences = "results/filtered-for-genome_{subtype}_{genome_seg}_{time}.fasta",
        reference = "config/reference_h5n1_{genome_seg}.gb",
    output:
        alignment =  "results/aligned-for-genome_{subtype}_{genome_seg}_{time}.fasta",
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
    # allow snakemake to choose the correct rule to run
    input:
        alignment = expand("results/aligned-for-genome_{{subtype}}_{genome_seg}_{{time}}.fasta",
                           genome_seg=["pb2", "pb1", "pa", "ha", "np", "na", "mp", "ns"]) 
    output:
        alignment = "results/aligned_{subtype}_{segment}_{time}.fasta"
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
        tree = "results/tree_{subtype}_{segment}_{time}.nwk",
        strains = "auspice/avian-flu_h5n1-cattle-outbreak_genome.json",
    output:
        tree = "results/tree_{subtype}_{segment}_{time}_outbreak-clade.nwk",
        node_data = "results/tree_{subtype}_{segment}_{time}_outbreak-clade.json",
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


def calc_genome_clock_rate(segment_rates, segment_lengths):
    mean = sum([cr * segment_lengths[seg] for seg,cr in segment_rates.items()])/sum(segment_lengths.values())
    stdev = mean/2
    return (mean, stdev)

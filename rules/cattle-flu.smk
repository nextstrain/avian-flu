
rule download_tree:
    """
    Downloads the tree behind nextstrain.org/avian-flu/h5n1-cattle-outbreak/genome
    so that the segment-level builds can use it to restrict the strains used.
    TODO: if the whole-genome analysis has been run locally we should (optionally) use that tree.
    """
    output:
        tree = "results/tree_{subtype}_genome.json",
    params:
        dataset="https://data.nextstrain.org/avian-flu_h5n1-cattle-outbreak_genome.json"
    wildcard_constraints:
        subtype="h5n1-cattle-outbreak",
        time="default",
    shell:
        """
        curl --compressed {params.dataset} -o {output.tree}
        """


rule prune_tree:
    input:
        tree = "results/tree_{subtype}_{segment}_{time}.nwk",
        strains = "results/tree_{subtype}_genome.json",
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

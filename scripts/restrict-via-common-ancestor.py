"""
Returns a monophyletic clade from the provided tree which encompasses all the (terminal)
strains provided via a list or an Auspice dataset.
"""

import argparse
import json
from Bio import Phylo
from sys import exit, stderr
from typing import cast, Set, Dict, Any
from augur.io import read_strains


def strains_in_auspice_json(fname: str) -> Set[str]:
    with open(fname) as fh:
        j = json.load(fh)
    assert isinstance(j['tree'], dict), "Only supports datasets with a single dataset"
    strains = set()
    nodes = [j['tree']]
    while len(nodes):
        node = nodes.pop(0)
        children = node.get("children", [])
        nodes.extend(children)
        if not children:
            strains.add(node['name'])
    return strains

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--tree', type=str, metavar="NEWICK", required=True, help="Input tree")
    parser.add_argument('--strains', type=str, metavar="TXT|JSON", required = True,
        help="A list of strains or Auspice JSON whose strains are used to define the common ancestor")
    parser.add_argument('--output-tree', type=str, metavar="NEWICK", required = True, help="output tree")
    parser.add_argument('--output-metadata', type=str, metavar="JSON", help="output node-data JSON")
    args = parser.parse_args()


    tree = Phylo.read(args.tree, "newick")
    if args.strains.lower().endswith(".txt"):
        all_defining_strains = cast(Set[str], read_strains(args.strains))
    elif args.strains.lower().endswith(".json"):        
        all_defining_strains = strains_in_auspice_json(args.strains)
    else:
        print("Unknown format of '--strains' - must be .txt or .json", file=stderr)
        exit(2)

    print(f"Defining strains in input strains: {len(all_defining_strains)}")
    defining_strains = [node for node in tree.get_terminals() if node.name in all_defining_strains]
    print(f"Defining strains found in input tree: {len(defining_strains)}")
    # print("DEFINING STRAINS IN TREE", defining_strains)

    if len(defining_strains)<2:
        print("ERROR: Not enough defining strains", file=stderr)
        exit(2)

    ca = tree.common_ancestor(defining_strains)
    print(f"Input tree root node: {tree.root.name!r}, num tips: {len(tree.get_terminals())}")
    print(f"Restricted tree root node: {ca.name!r}, num tips: {len(ca.get_terminals())}")

    # TODO -- (the root of) `ca` may be part of a polytomy when viewed using divergence, and if so we probably
    # want to walk up the tree and return the entire polytomy. This is consistent with the general idea of
    # returning terminal nodes which may be basal to the "true" cattle-flu outbreak

    Phylo.write(ca, args.output_tree, "newick")

    if (args.output_metadata):
        meta: Dict[str, Any] = {"nodes": {}}
        # TODO add "generated_by" section once <https://github.com/nextstrain/augur/issues/1476>
        # is resolved

        for node in ca.get_terminals():
            # booleans don't flow nicely through to Auspice so use strings here
            meta['nodes'][node.name] = {'genome_tree': "Yes" if node.name in all_defining_strains else "No"}

        # document the strains in the genome tree (--strains) but which are missing from the
        # input tree. Ideally there are none!
        strains_in_tree = set([node.name for node in ca.get_terminals()])
        meta['nodes']['extensions'] = {
            'genome_strains_missing_from_segment_tree': [x for x in all_defining_strains if x not in strains_in_tree]
        }

        with open(args.output_metadata, 'w') as fh:
            json.dump(meta, fh)

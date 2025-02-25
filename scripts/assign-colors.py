"""
Generates a colors TSV (on STDOUT) for use in `augur export`.

Colors are generated from the provided color-scheme via an ordering TSV which
defines the expected values for each category. Only those values which are
observed in the (provided) metadata file are assigned colors, and this means
that colors are not stable over time or datasets. Observations which are missing
from the ordering TSV will print a warning on STDERR and be sorted to the end of
the order for color-assignment.
"""


import argparse
from sys import stderr, stdout
from augur.io import read_metadata
from collections import defaultdict
import csv

def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--ordering', type=str, required=True, nargs="+", metavar="TSV",
                        help="input ordering file(s). First column: category name, second column: value. Multiple files are essentially concatenated together.")
    parser.add_argument('--color-schemes', type=str, required=True, metavar="TSV",
                        help="input color schemes file. Each line is a list of colour hexes, in order. Each line should have a different number of entries.")
    parser.add_argument('--metadata', type=str, required=True, metavar="TSV",
                        help="restrict colors to only those found in metadata")
    parser.add_argument('--duplications', type=str, required=False, nargs="*", default=[], metavar="categoryX=categoryY",
                        help="Duplicate the colours generated for categoryX under the name categoryY")
    return parser.parse_args()

def read_ordering(fnames):
    # copied from <https://github.com/nextstrain/ncov/blob/80c06c15f0669c5990e0cc03628deeeedca7b720/scripts/assign-colors.py>
    # with modifications
    assignment = defaultdict(list)
    for fname in fnames:
        with open(fname) as f:
            for line in f.readlines():
                array = line.lstrip().rstrip().split("\t")
                if len(array) == 2 and not line.startswith("#"):
                    name = array[0]
                    trait = array[1]
                    assignment[name].append(trait)
    return assignment

def read_schemes(fname):
    # copied from <https://github.com/nextstrain/ncov/blob/80c06c15f0669c5990e0cc03628deeeedca7b720/scripts/assign-colors.py>
    schemes = {}
    counter = 0
    with open(args.color_schemes) as f:
        for line in f.readlines():
            counter += 1
            array = line.lstrip().rstrip().split("\t")
            schemes[counter] = array
    return schemes

def index_of(search_list, search_value, not_found_value):
    try:
        return search_list.index(search_value)
    except ValueError:
        return not_found_value

if __name__ == '__main__':
    args = parse_args()
    metadata = read_metadata(args.metadata)
    ordering = read_ordering(args.ordering)
    schemes = read_schemes(args.color_schemes)
    tsv_writer = csv.writer(stdout, delimiter='\t', lineterminator='\n')
    duplications = {x.split('=')[0]: x.split('=')[1] for x in args.duplications}

    for category in ordering:
        try:
            observations = list(metadata[category].unique())
        except KeyError:
            print(f"WARNING: Ordering specified for '{category}' but this wasn't found in the metadata")
            continue
        observations.sort(key=lambda x: index_of(ordering[category], x, len(ordering[category])))
        if (missing:=set(observations) - set(ordering[category])):
            print(f"WARNING: For trait {category} we encountered {len(missing)} value(s) not present in the ordering TSV: {missing}", file=stderr)
        try:
            hexes = schemes[len(observations)]
        except KeyError:
            raise Exception(f"Ordering '{category}' had {len(observations)} values but the color-schemes don't go this high!")
        for pair in zip(observations, hexes):
            tsv_writer.writerow([category, *pair])

        if category in duplications:
            for pair in zip(observations, hexes):
                tsv_writer.writerow([duplications[category], *pair])


import argparse
from sys import stderr, stdout
from augur.io import read_metadata
import csv

def parse_args():
    parser = argparse.ArgumentParser(
        description="Assign colours based on ordering and metadata TSVs. Outputs a TSV on STDOUT.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--ordering', type=str, required=True, metavar="TSV",
                        help="input ordering file. First column: category name, second column: value")
    parser.add_argument('--color-schemes', type=str, required=True, metavar="TSV",
                        help="input color schemes file. Each line is a list of colour hexes, in order. Each line should have a different number of entries.")
    parser.add_argument('--metadata', type=str, required=True, metavar="TSV",
                        help="restrict colors to only those found in metadata")
    parser.add_argument('--duplications', type=str, required=False, nargs="+", metavar="categoryX=categoryY",
                        help="Duplicate the colours generated for categoryX under the name categoryY")
    return parser.parse_args()

def read_ordering(fname):
    # copied from <https://github.com/nextstrain/ncov/blob/80c06c15f0669c5990e0cc03628deeeedca7b720/scripts/assign-colors.py>
    assignment = {}
    with open(fname) as f:
        for line in f.readlines():
            array = line.lstrip().rstrip().split("\t")
            if len(array) == 2:
                name = array[0]
                trait = array[1]
                if name not in assignment:
                    assignment[name] = [trait]
                else:
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
            raise Exception(f"Ordering specified for '{category}' but this wasn't found in the metadata")
        observations.sort(key=lambda x: index_of(ordering[category], x, len(ordering[category])))
        if (missing:=set(observations) - set(ordering[category])):
            print(f"WARNING: For trait {category} we encountered {len(missing)} value(s) not present in the ordering TSV: {missing}", file=stderr)
            # -- temporory --
            # reduce the values we provide colorings for to just those present in the orderings TSV
            # (values pruned out here will still be shown in Auspice, but they'll be grey)
            # This addresses bug <https://github.com/nextstrain/avian-flu/issues/131> where
            # we had more division values than colors and thus the pipeline failed
            # Hopefully the underlying cause will be fixed by improvements in how we ingest data, e.g.
            # <https://github.com/nextstrain/augur/issues/1578>
            # and we can remove this restriction.
            observations = [obs for obs in observations if obs in ordering[category]]
            print(f"These will not be assigned colors and thus will appear grey in Auspice", file=stderr)
        try:
            hexes = schemes[len(observations)]
        except KeyError:
            raise Exception(f"Ordering '{category}' had {len(observations)} values but the color-schemes don't go this high!")
        for pair in zip(observations, hexes):
            tsv_writer.writerow([category, *pair])

        if category in duplications:
            for pair in zip(observations, hexes):
                tsv_writer.writerow([duplications[category], *pair])

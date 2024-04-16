"""
Takes in a set of metadata TSVs corresponding to segments (i.e. typically 8 TSVs)
and adds a column to the input `--metadata` TSV with the number of segments
that strain appears in.
"""

import argparse
import csv
from collections import defaultdict

def collect_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--segments', type=str, nargs='+', help='Metadata TSVs for all segments')
    parser.add_argument('--metadata', type=str, help='Metadata TSV which will be amended for output. Must also appear in --segments.')
    parser.add_argument('--output',   type=str, help='metadata file name')
    return parser.parse_args()

def read_metadata(fname, strains_only=False):
    strains = set()
    rows = []
    with open(fname) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            strains.add(row['strain'])
            if not strains_only:
                rows.append(row)
    if strains_only:
        return strains
    return (strains, reader.fieldnames, rows)

def summary(strain_count):
    ## Print some basic stats!
    print("Num strains observed (across all segments): ", len(strain_count.keys()))
    counts = [0]*9 # 1-indexed
    for n in strain_count.values():
        counts[n]+=1
    for i in range(1,9):
        print(f"Num strains observed in {i} segments: ", counts[i])


if __name__=="__main__":
    args = collect_args()
    strain_count = defaultdict(int)
    matching_metadata_file_provided = False

    for fname in args.segments:
        if fname==args.metadata:
            matching_metadata_file_provided = True
            _strains, fieldnames, rows = read_metadata(fname)
        else:
            _strains = read_metadata(fname, strains_only=True)
        for s in _strains:
            strain_count[s]+=1
    if not matching_metadata_file_provided:
        raise Exception("You must provide a --metadata file which exactly matches one of the provided --segments files")
    
    summary(strain_count)

    # append count to data for output
    column = "n_segments"
    fieldnames.append(column)
    for row in rows:
        row[column]=strain_count[row['strain']]
    
    with open(args.output, 'w') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print("Output written to", args.output)

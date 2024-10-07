"""
Given a metadata TSV on STDIN, write out the same metadata TSV but with
an extra column "n_segments" which is the sum of all "segment_X" columns
for all X in --segments
TODO - use proper quoting once it's been decided
"""

import csv
import argparse
from sys import stdin,stdout

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument('--segments', metavar="SEGMENT", nargs="+", required=True, help="Segment names")
    args = parser.parse_args()
    reader = csv.DictReader(stdin, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
    writer = csv.DictWriter( # copy/paste from the (unexported) `augur.io.metadata.write_records_to_tsv` 
        stdout,
        [*list(reader.fieldnames), "n_segments"],
        extrasaction='ignore',
        delimiter='\t',
        lineterminator='\n',
        quoting=csv.QUOTE_NONE,
        quotechar=None,
    )
    writer.writeheader()
    for record in reader:
        record['n_segments'] = sum([int(record.get(f"segment_{segment}", '0'))==1 for segment in args.segments])
        writer.writerow(record)

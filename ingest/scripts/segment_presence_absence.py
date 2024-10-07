"""
Given a metadata TSV on STDIN, write a 2-column TSV to STDOUT.
Columns: "strain" and "segment_{segment}" where the segment values "1"
(i.e. any data in the input metadata indicates the segment is present)
TODO - use proper quoting once it's been decided
"""

import csv
import argparse
from sys import stdin,stdout

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument('--segment', metavar="SEGMENT", required=True, help="Segment name")
    args = parser.parse_args()
    reader = csv.DictReader(stdin, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
    print(f"strain\tsegment_{args.segment}", file=stdout)
    for row in reader:
        print(f"{row['strain']}\t1", file=stdout)
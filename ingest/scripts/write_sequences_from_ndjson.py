
"""
Write FASTA sequences from NDJSON

TKTK - Is there a way to do this using simpler tools patched / piped together?
"""

import json
import argparse
from io import TextIOWrapper

def parse_args():
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument('--input', required=True,  metavar="NDJSON", type=str)
    parser.add_argument('--id-key', required=False, nargs='+', type=str, default="strain")
    parser.add_argument('--output', required=True, nargs='+', metavar="SEGMENT=FIILENAME", type=str)
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    segments: dict[str,TextIOWrapper] = {}
    for el in args.output:
        fields = el.split('=')
        if len(fields)!=2:
            raise Exception(f"Unknown output argument '{el}'")
        segments[fields[0]] = open(fields[1], 'w')

    with open(args.input) as metadata_fh:
        for line in metadata_fh:
            record = json.loads(line)
            strain = record[args.id_key]
            for segment,fh in segments.items():
                sequence = record.get(f"sequence_{segment}", None)
                if not sequence:
                    continue
                print(f">{strain}\n{sequence}", file=fh)

    for fh in segments.values():
        fh.close()

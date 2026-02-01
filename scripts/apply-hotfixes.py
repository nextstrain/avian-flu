#!/usr/bin/env python3

"""
Apply hotfixes to metadata TSV files.

Reads a hotfixes file containing strain-level field corrections and applies
them to the input metadata. Hotfixes that are no longer needed (because the
upstream value already matches) are reported and excluded from the output.
"""


import argparse
from sys import stderr
from augur.io import read_metadata
import pandas as pd


def read_tsv(filepath: str) -> list[str]:
    """Read a TSV file line-by-line, preserving all lines including comments."""
    with open(filepath, 'r') as f:
        return f.read().splitlines()

Hotfixes = dict[str, dict[str, dict[str, str|int]]]


def parse_lines(lines: list[str]) -> Hotfixes:
    """
    Parse hotfix lines into a nested dictionary structure.

    Returns hotfixes[strain_name][field_name] = {'value': field_value, 'idx': idx}
    """
    hotfixes: Hotfixes = {}
    for idx, line in enumerate(lines):
        # Skip comments and empty lines
        if line.startswith('#') or not line.strip():
            continue

        parts = line.split('\t')
        if len(parts) < 3:
            print(f"WARNING: Skipping malformed hotfix line {idx}: {line}", file=stderr)
            continue

        strain_name, field_name, field_value = parts[0], parts[1], parts[2]

        if strain_name not in hotfixes:
            hotfixes[strain_name] = {}

        hotfixes[strain_name][field_name] = {'value': field_value, 'idx': idx}

    return hotfixes



def apply_fixes(metadata: pd.DataFrame, hotfixes: Hotfixes, unnecessary_fixes: list[int]) -> None:
    """
    Apply hotfixes to the metadata DataFrame in place.

    Tracks unnecessary fixes (where the value already matches) in the provided list.
    """
    for strain_name, fields in hotfixes.items():
        if strain_name not in metadata.index:
            print(f"WARNING: Strain '{strain_name}' not found in metadata", file=stderr)
            continue

        for field_name, fix_info in fields.items():
            fix_value = fix_info['value']
            idx = fix_info['idx']
            assert isinstance(idx, int)

            if field_name not in metadata.columns:
                print(f"WARNING: Strain '{strain_name}' field '{field_name}' not found in metadata columns", file=stderr)
                continue

            current_value = metadata.loc[strain_name, field_name]

            # Handle NaN comparisons
            if pd.isna(current_value) and pd.isna(fix_value):
                print(f"Unnecessary fix: {strain_name} {field_name} already NaN", file=stderr)
                unnecessary_fixes.append(idx)
            elif str(current_value) == str(fix_value):
                print(f"Unnecessary fix: {strain_name} {field_name} already '{current_value}'", file=stderr)
                unnecessary_fixes.append(idx)
            else:
                metadata.loc[strain_name, field_name] = fix_value
                print(f"Applied fix: {strain_name} {field_name} changed from '{current_value}' to '{fix_value}'", file=stderr)


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--metadata', type=str, required=True, metavar="TSV", help="input metadata")
    parser.add_argument('--hotfixes', type=str, required=True, metavar="TSV", help="hotfixes TSV")
    parser.add_argument('--output', type=str, required=True, metavar="TSV", help="output metadata")
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    metadata: pd.DataFrame = read_metadata(args.metadata)

    unnecessary_fixes: list[int] = []

    # read the hotfixes TSV file, line-by-line, into a list
    lines = read_tsv(args.hotfixes)
    # then extract the hotfixes out of it
    # hotfixes[strain_name][field_name] = {value: field_value, idx: idx} where idx is a pointer to the appropriate lines[idx]
    hotfixes = parse_lines(lines)

    # and apply the hotfixes to the `metadata`, with lots of logging
    apply_fixes(metadata, hotfixes, unnecessary_fixes)

    # write out the updated metadata 
    metadata.to_csv(args.output, sep='\t')

    # and overwrite the updated hotfix lines, exluding those in unnecessary_fixes
    with open(args.hotfixes, 'w') as fh:
        for idx, line in enumerate(lines):
            if idx not in unnecessary_fixes:
                print(line, file=fh)

"""
Transforms a metadata file with one row per accession into a file with
one row per strain. Where a strain contains sequences for multiple segments this
will group those segments together.

Segment-specific field names (i.e. those not in --common-strain-fields) are modified
to ensure their suffix is "_{segment}". For instance "accession → accession_HA".

Rows are matched on strain name and basic sanity checking is performed when grouping.
Strains with zero or multiple matches for any segment are dropped entirely.
Any disagreement within the "--common-strain-fields" will result in the strain being
dropped, however empty values may be replaced and ambiguous dates may be replaced with
specific ones (where appropriate).

There is no way to resolve conflicts yet, however this should be added in the future.
"""

import argparse
import csv
import yaml
from collections import defaultdict
from sys import stderr, exit
from typing import Any
from typing_extensions import TypedDict, NotRequired

def parse_args():
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument('--common-strain-fields', metavar="NAME", required=True, nargs='+', type=str,
                        help="Fields which must match among sequences of a single strain")
    parser.add_argument('--segments', metavar="NAME", required=True, nargs='+', type=str,
                        help="Segment names")
    parser.add_argument('--metadata', metavar="TSV", required=True, type=str,
                        help="Input metadata file. ID column='accession'")
    parser.add_argument('--resolutions', metavar='YAML', required=False,
                        help="Rules to resolve conflicts when grouping")
    parser.add_argument('--output-metadata', metavar="TSV", required=True, type=str,
                        help="Input metadata file. ID column='strain'")
    return parser.parse_args()

def group_by_strain(filename: str) -> dict[str,list]:
    strains = defaultdict(list)
    with open(filename) as fh:
        reader = csv.DictReader(fh, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
        for row in reader:
            strains[row['strain']].append(row)
    return strains


class Resolutions(TypedDict):
    PICK_SEGMENT: NotRequired[list[dict]]
    PICK_FIELD: NotRequired[list[dict]]

def parse_resolutions_yaml(fname:str) -> Resolutions:
    with open(fname) as fh:
        try:
            resolutions = yaml.safe_load(fh)
        except yaml.YAMLError as e:
            print(e)
            exit(2)
    return resolutions

def log(msg)->None:
    # Currently we just dump to STDERR, but this should be formalised
    print(msg, file=stderr)

def get_segment(strain:str, row:dict[str,Any], segment_names: list[str])->str:
    # segment fields are encoded like 'segment_S' with a value of 1 if present.
    segments_present = [seg for seg in segment_names if row[f'segment_{seg}']=="1"]
    accession = row.get('accession', 'unknown')
    if len(segments_present)==0:
        raise AssertionError(f"Accession '{accession}' (strain '{strain}') mapped to no segments. Skipping this accession.")
    if len(segments_present)>1:
        raise AssertionError(f"Accession '{accession}' (strain '{strain}') mapped to multiple segments: {', '.join(segments_present)}. Skipping this accession.")
    return next(iter(segments_present))


def resolve_segment(resolutions: Resolutions, strain:str, segment:str, accessions: list[str]) -> str|None:
    rules = [el for el in resolutions.get('PICK_SEGMENT', []) if el.get('strain')==strain and el.get('segment')==segment]
    if len(rules)==0:
        return None
    if len(rules)>1:
        log(f"Malformed resolutions YAML - multiple PICK_SEGMENT blocks for strain={strain} segment={segment}")
        exit(2)
    rule = rules[0]
    if rule['accession'] not in accessions:
        log(f"ERROR! A PICK_SEGMENT resolution for strain {strain} for segment {segment} specified an accession which wasn't in the metadata.")
        return None
    return rule['accession']

def resolve_field(resolutions: Resolutions, strain:str, field:str, accessions: list[str]) -> str|None:
    # TODO - this is essentially identical to `resolve_segment`. If they don't diverge, we should consolidate.
    rules = [el for el in resolutions.get('PICK_FIELD', []) if el.get('strain')==strain and el.get('field')==field]
    if len(rules)==0:
        return None
    if len(rules)>1:
        log(f"Malformed resolutions YAML - multiple PICK_FIELD blocks for strain={strain} field={field}")
        exit(2)
    rule = rules[0]
    if rule['accession'] not in accessions:
        log(f"ERROR! A PICK_FIELD resolution for strain {strain} for field {field} specified an accession which wasn't in the metadata.")
        return None
    return rule['accession']


def assign_segments(strain_name:str, rows:list, segment_names:list[str], resolutions: Resolutions)->dict[str,dict]|None:
    """
    Given rows (assigned to a strain) assign each to a segment. Error if
    (1) The same row (sequence) is assigned to multiple sequences
    (2) Multiple rows are assigned to the same segment
    NOTE: In the future we should resolve (2). TODO XXX
    """
    segments = {}
    rows_by_segment = {segment: list() for segment in segment_names}
    for row in rows:
        try:
            segment_name = get_segment(strain_name, row, segment_names)
        except AssertionError as e:
            log(e)
            continue # ignore this row
        rows_by_segment[segment_name].append(row)
    # Drop any segments with more than one matching accession / sequence unless there's a rule to resolve it
    for segment_name, seg_rows in rows_by_segment.items():
        accessions = [r['accession'] for r in seg_rows]
        if len(seg_rows)==1:
            segments[segment_name] = seg_rows[0]
        elif len(seg_rows)>1:
            if (accession := resolve_segment(resolutions, strain_name, segment_name, accessions)):
                segments[segment_name] = next(iter([row for row in seg_rows if row['accession']==accession]))
                log(f"Resolving '{strain_name}' to use accession {accession} for segment {segment_name} ")
                continue
            log(f"Strain '{strain_name}' had multiple accessions for segment {segment_name}. "
                f"Accessions: {', '.join(accessions)}. "
                "Skipping this segment.")
            continue # ignore this segment (other segments for this strain may be OK)

    if len(segments)==0:
        log(f"Strain '{strain_name}' had zero or multiple accessions for all segments. Dropping this entire strain.")
        return None
    return segments

class ValueMatchingError(Exception):
    pass

class HeaderInfo(TypedDict):
    fields: list[str]
    common: list[str]
    segment_specific: list[dict[str,str]]


def pick_from_values(strain_name:str, field_name:str, rows:list, resolutions: Resolutions, allow_empty=True)->str:
    values = set(row[field_name] for row in rows)
    if allow_empty and "" in values and len(values)!=1:
        values.remove("")
    if len(values)>1:
        if field_name=='date':
            try:
                return resolve_mismatch_dates(strain_name, values)
            except ValueMatchingError:
                # continue, and use the error message printing below
                pass

        if (accession:=resolve_field(resolutions, strain_name, field_name, [row['accession'] for row in rows])):
            value = next(iter([row[field_name] for row in rows if row['accession']==accession]))
            log(f"Resolving '{strain_name}' to use {field_name}={value}")
            return value

        # want to print out helpful messages about disagreement, so order by most commonly observed
        obs = defaultdict(list)
        for row in rows:
            obs[row[field_name]].append(row['accession'])
        msg = f"Strain '{strain_name}' Disagreement for '{field_name}', {len(obs)} observed values:"
        for v,acc in sorted(obs.items(), key=lambda item: item[1], reverse=True):
            msg+=f"\n\t{', '.join(acc)}: {v}"
        raise ValueMatchingError(msg)
    return values.pop()

def resolve_mismatch_dates(strain_name:str, values:set[str])->str:
    if "XXXX-XX-XX" in values:
        values.remove("XXXX-XX-XX")
    # TODO XXX - this function is incomplete
    if len(values)==1:
        return values.pop()
    raise ValueMatchingError()


def make_wide(strain: str, rows: list, segment_names: list[str], header_info:HeaderInfo, resolutions: Resolutions) -> dict[str,str]|None:
    segments = assign_segments(strain, rows, segment_names, resolutions)
    if not segments:
        return None

    metadata = {
        "strain": strain,
        "n_segments": str(len(segments.keys())),
    }

    observed_mismatches = False
    for field_name in header_info['common']:
        try:
            metadata[field_name] = pick_from_values(strain, field_name, list(segments.values()), resolutions)
        except ValueMatchingError as e:
            log(e)
            observed_mismatches = True
    if observed_mismatches:
        return None

    for info in header_info['segment_specific']:
        if info['segment'] not in segments:
            continue
        metadata[info['output_field']] = segments[info['segment']][info['input_field']]

    return metadata


def header(segments:list[str], tsv_fields:list[str], common_strain_fields:list[str])->HeaderInfo:

    assert all([f in tsv_fields for f in common_strain_fields]), \
        "Names in '--common-strain-fields' must be present in the TSV"

    segment_specific = []
    for segment in segments:
        for field in [f for f in tsv_fields if f not in common_strain_fields and f!="strain"]:
            # if the field is "X_Y" and Y is a segment then handle separately
            if "_" in field and field.split('_')[-1] in segments:
                if field.split('_')[-1]==segment:
                    segment_specific.append({"segment": segment, "input_field": field, "output_field": field})
                continue
            segment_specific.append({"segment": segment, "input_field": field, "output_field": f"{field}_{segment}"})

    return {
        "fields": ["strain", "n_segments", *common_strain_fields, *[x['output_field'] for x in segment_specific]],
        "common": common_strain_fields,
        "segment_specific": segment_specific
    }

if __name__=="__main__":
    args = parse_args()
    resolutions = parse_resolutions_yaml(args.resolutions) if args.resolutions else {}

    strains = group_by_strain(args.metadata)

    header_info = header(args.segments, list(next(iter(strains.values()))[0].keys()), args.common_strain_fields)

    collapsed = [row
                 for row in [make_wide(strain, rows, args.segments, header_info, resolutions) for strain,rows in strains.items()]
                 if row is not None]
    log("If any errors have been printed above, then those strains will have been dropped.")

    with open(args.output_metadata, 'w', newline='') as output_metadata:
        # DictWriter is what `augur.io.metadata.write_records_to_tsv` uses
        # (we can switch to that if we ensure the field order in the first row is correct)
        tsv_writer = csv.DictWriter(
            output_metadata,
            header_info['fields'],
            extrasaction='ignore',
            delimiter='\t',
            lineterminator='\n',
            quoting=csv.QUOTE_NONE,
            quotechar=None,
        )
        tsv_writer.writeheader()
        for record in collapsed:
            tsv_writer.writerow(record)

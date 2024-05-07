"""
Curate the metadata that originated from Andersen Lab's avian-influenza repo
<https://github.com/andersen-lab/avian-influenza>.

Parses NDJSON record from stdin and outputs new record to stdout.
"""
import copy
import json
from datetime import datetime
from enum import Enum
from sys import stdin, stdout, stderr


NEXTSTRAIN_RECORD = {
    'strain': '?',
    'virus': 'avian_flu',
    'isolate_id': '?',
    'date': '?',
    'region': 'North America',
    'country': 'USA',
    'division': '?',
    'location': '?',
    'host': '?',
    'domestic_status': '?',
    'subtype': 'h5n1',
    'originating_lab': '?',
    'submitting_lab': '?',
    'authors': '?',
    'PMID': '?',
    'gisaid_clade': '?',
    'h5_clade': '?',
}


def create_new_record(anderson_record: dict) -> dict:
    """
    Create a new NEXTSTRAIN_RECORD with additional data from the provided
    `andersen_record`.
    """
    new_record = copy.deepcopy(NEXTSTRAIN_RECORD)
    new_record['isolate_id'] = anderson_record['Run']
    new_record['division'] = anderson_record['US State']
    new_record['location'] = anderson_record['US State']

    center_name = parse_center_name(anderson_record['Center Name'])
    new_record['originating_lab'] = center_name
    new_record['submitting_lab'] = center_name

    new_record['host'] = parse_host_group(anderson_record['Host'])
    new_record['date'] = parse_date(anderson_record['Date'])
    new_record['strain'] = f'A/{anderson_record["Host"]}/{new_record["country"]}/{anderson_record["isolate"]}/{parse_year(new_record["date"])}'
    return new_record


def parse_center_name(center_name: str) -> str:
    if center_name == 'USDA-NVSL':
        return center_name.replace('-', ' ')

    return center_name


def parse_host_group(host: str) -> str:
    """
    Bin `host` into the HOST_GROUPS
    """
    # Replace with enum.StrEnum starting with Python 3.11
    class HOST_GROUPS(str, Enum):
        AVIAN = 'Avian'
        CATTLE = 'Cattle'
        NONHUMAN_MAMMAL = 'Nonhuman Mammal'

    known_hosts = {
        'Blackbird': HOST_GROUPS.AVIAN,
        'Cat': HOST_GROUPS.NONHUMAN_MAMMAL,
        'Cattle': HOST_GROUPS.CATTLE,
        'CAGO': HOST_GROUPS.AVIAN,
        'Chicken': HOST_GROUPS.AVIAN,
        'Grackle': HOST_GROUPS.AVIAN,
        'Goose': HOST_GROUPS.AVIAN,
        'PEFA': HOST_GROUPS.AVIAN,
        'Skunk': HOST_GROUPS.NONHUMAN_MAMMAL,
        'Raccoon': HOST_GROUPS.NONHUMAN_MAMMAL,
    }

    host_group = known_hosts.get(host)
    if host_group is None:
        print(f"WARNING: unable to group unknown host {host!r}", file=stderr)
        return host

    return host_group


def parse_date(date_string: str) -> str:
    """
    If date_string is empty, 'NA', or includes a `?`, then returns the
    hardcoded `2024-XX-XX` date.

    Otherwise return the original date_string.
    """
    default_date = '2024-XX-XX'

    if not date_string or date_string == 'NA' or '?' in date_string:
        return default_date

    return date_string


def parse_year(date_string: str) -> str:
    """
    Parse the year from the provided `date_string`
    """
    date_formats = ['%Y-%m-%d', '%Y-XX-XX']
    for date_format in date_formats:
        try:
            parsed_date = datetime.strptime(date_string, date_format)
            return str(parsed_date.year)
        except ValueError:
            continue

    raise ValueError(f"Could not parse year from date string {date_string!r}")


if __name__ == '__main__':

    for record in stdin:
        anderson_record = json.loads(record)
        new_record = create_new_record(anderson_record)
        json.dump(new_record, stdout, allow_nan=False, indent=None, separators=',:')
        print()
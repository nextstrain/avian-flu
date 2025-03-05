"""
Takes a (modified) GenoFLU results TSV on STDIN and writes a TSV to STDOUT
The input TSV is expected to have three fields:
* strain
* genoflu
* details
The output TSV exports 10 fields:
* strain
* genoflu
* genoflu_<SEGMENT> (8 fields)
"""


from augur.io.metadata import read_table_to_dict
from sys import stdin, stdout
from csv import DictWriter

if __name__ == "__main__":
    SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
    HEADER = ['strain', 'genoflu', *[f"genoflu_{s}" for s in SEGMENTS]]
    tsv_writer = DictWriter(stdout,HEADER,extrasaction='ignore',delimiter='\t',lineterminator='\n')
    tsv_writer.writeheader()
    for record in read_table_to_dict(stdin.buffer, ["\t"]):
        if record['details']:
            for segment_name, lineage in (parts.split(':') for parts in record['details'].split(", ")):
                record[f"genoflu_{segment_name}"] = lineage
        # There are a variety of situations where a genotype is not assigned
        # (1) segments have all been assigned, but there's no matching constellation
        if record['genoflu'] == "Not assigned: No Matching Genotypes":
            record['genoflu'] = "Unseen constellation"
        # (2) Not all 8 segments have assignments (some segments were too divergent)
        elif record['genoflu'].startswith("Not assigned: Only"):
            record['genoflu'] = "Not assigned (too divergent)"
        # (3) TODO Ideally the results would also include samples with <8 segments sequenced
        tsv_writer.writerow(record)

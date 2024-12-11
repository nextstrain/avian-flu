
from Bio import SeqIO
from collections import defaultdict
import argparse
import json


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--segments', type = str, required = True, nargs='+', help = "per-segment alignments")
    parser.add_argument('--output', type = str, required = True, help = "output whole genome alignment")
    parser.add_argument('--output-node-data', type = str, required = False, help = "output metadata in node-data JSON format")
    args = parser.parse_args()

    records = {}
    strain_counts: dict[str,int] = defaultdict(int)
    segment_lengths = defaultdict(int)
    atgc = set(['A','T','G','C'])

    assert len(args.segments)==8, "Expected 8 segments!"

    for fname in args.segments:
        records[fname] = {record.name:str(record.seq).upper() for record in SeqIO.parse(fname, 'fasta')}
        for key in records[fname]:
            strain_counts[key]+=1
        assert len({len(seq) for seq in records[fname].values()})==1, f"Different sequence lengths observed in {fname}"
        segment_lengths[fname] = len(next(iter(records[fname].values())))
        print(f"{fname}: parsed {len(records[fname].keys())} sequences, each {segment_lengths[fname]} nt")

    ## how many strains are missing segments?
    num_segments: dict[str, list[str]] = {str(i): [] for i in range(1,8+1)}
    for strain,n in strain_counts.items():
        num_segments[str(n)].append(strain)
        

    def sequence(segment, name):
        if seq:=records[segment].get(name, None):
            return seq
        # `augur ancestral` is run with --keep-ambiguous but without --keep-overhangs
        # So use Ns to represent missing segments rather than gaps
        # https://docs.nextstrain.org/en/latest/guides/bioinformatics/missing-sequence-data.html
        return "N" * segment_lengths[segment]
    
    def atgc_perc(seq):
        atgc = set(['A','T','G','C'])
        len([nt for nt in seq if nt in atgc])

    node_data: dict[str,dict] = {'nodes': {}}

    with open(args.output, 'w') as fh:
        print("writing genome to ", args.output)
        for name,count in strain_counts.items():
            if count<7:
                print(f"Excluding {name} as it only appears in {count} segments")
                continue
            genome = "".join([sequence(seg, name) for seg in args.segments])
            node_data['nodes'][name] = {
                "ATGC_perc": int( len([nt for nt in genome if nt in atgc])/len(genome) * 100),
                "num_segments": strain_counts[name]
            }
            print(f">{name}\n{genome}", file=fh)

    if args.output_node_data:
        with open(args.output_node_data, 'w') as fh:
            json.dump(node_data, fh)
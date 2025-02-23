
import argparse
from augur.io.sequences import read_sequences, write_sequences
from Bio.Seq import MutableSeq

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--alignment', type=str, required=True, help='genome alignment')
    parser.add_argument('--percentage', type=str, required=True, help='positions with less coverage than this will be masked')
    parser.add_argument('--output', type=str, required=True, help='masked genome alignment')

    args = parser.parse_args()

    # Store everything in memory
    alignment = list(read_sequences(args.alignment))
    genome_size = len(alignment[0].seq)
    counts = [0 for _ in range(0, genome_size)] # zero-based
    valid_bases = set(list("ATGCatcg"))
    n_genomes = len(alignment)

    for sample in alignment:
        for idx, base in enumerate(sample.seq):
            if base in valid_bases:
                counts[idx] += 1

    mask_bool = [c/n_genomes*100 < float(args.percentage) for c in counts]
    mask_sites = [i for i,b in enumerate(mask_bool) if b==True]

    print("Masking sites (zero-based):", mask_sites)
    print("Total number of sites to mask:", len(mask_sites))

    for sample in alignment:
        sample.seq = MutableSeq(sample.seq)
        for idx in mask_sites:
            sample.seq[idx] = 'N'

    write_sequences(alignment, args.output)
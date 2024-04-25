
import argparse
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SimpleLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from augur.utils import load_features

## GenBank coordinate systems           Python feature coordinate system: 0-based [start, end)
## start: 1-based, inclusive            0-based, inclusive (genbank start number -1)
## end: 1-based, inclusive              0-based, exclusive (same number as genbank end)
## length: end - start + 1              end - start


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--genbank', type = str, required = True, nargs='+', help = "per-segment genbank files (in order)")
    parser.add_argument('--output', type = str, required = True, help = "output whole genome alignment")
    args = parser.parse_args()

    shift_bp = 0
    record = None
    seqs = []

    for file in args.genbank:
        features = load_features(file) # enforces 'nuc' to be present
        gb = SeqIO.read(file, 'genbank') # read a second time to get info not provided by augur's `load_features`
        seqs.append(gb.seq)
        if shift_bp == 0:
            print(f"Using {file} unmodified ({features['nuc'].location.end}nuc) ({', '.join(features.keys())})\n")
            record = SeqRecord(
                Seq(''),
                id=gb.id,
                name=gb.name,
                features=[features['nuc'], *[feat for name, feat in features.items() if name != 'nuc']], # nuc 1st
                description=f"Combined references for {len(args.genbank)} segments",
                annotations=gb.annotations,
            )
            shift_bp += features['nuc'].location.end
            continue
        
        # Segment 2 (or 3, 4, 5...)
        print(f"Reading {file} ({features['nuc'].location.end}nuc) and shifting by {shift_bp}nuc")
        segment_len = features['nuc'].location.end
        # Update nuc feature
        nuc = record.features[0]
        nuc.location = FeatureLocation(nuc.location.start, nuc.location.end + segment_len, nuc.location.strand)
        # Add shifted segment features to record
        for name,feat in features.items():
            if name=='nuc': continue

            if isinstance(feat.location, SimpleLocation):
                assert feat.location.start < feat.location.end
                assert feat.location.strand == 1 # not yet considering -ve strand CDSs
                _previous_coords = f"[{feat.location.start+1}, {feat.location.end}]" # genbank coordinate system
                feat.location = FeatureLocation(feat.location.start+shift_bp, feat.location.end+shift_bp, feat.location.strand)
                record.features.append(feat)
                print(f"\t{name} shifting from {_previous_coords} to [{feat.location.start+1}, {feat.location.end}] (genbank coordinates)")
            elif isinstance(feat.location, CompoundLocation):
                shifted_parts = []
                _previous_coords = []
                _shifted_coords = []
                for part in feat.location.parts:
                    assert isinstance(part, SimpleLocation)
                    assert part.start < part.end
                    assert part.strand == 1 # not yet considering -ve strand CDSs
                    _previous_coords.append(f"[{part.start+1}, {part.end}]") # genbank coordinate system
                    _shifted_part = FeatureLocation(part.start+shift_bp, part.end+shift_bp, part.strand)
                    _shifted_coords.append(f"[{_shifted_part.start+1}, {_shifted_part.end}]") # genbank coordinate system
                    shifted_parts.append(_shifted_part)
                feat.location = CompoundLocation(shifted_parts)
                record.features.append(feat)
                print(f"\t{name} shifting from join({', '.join(_previous_coords)}) to join({', '.join(_shifted_coords)}) (genbank coordinates)")
            else:
                raise Exception("Unknown location type")

        shift_bp += features['nuc'].location.end
        print()

    record.seq = record.seq.join(seqs)

    print(f"Total sequence length: {len(record.seq)}")
    SeqIO.write(record, args.output, 'genbank')
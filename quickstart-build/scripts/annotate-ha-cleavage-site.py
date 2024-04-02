"""
This script will read in the HA alignment file, translate the sequence to amino acids, 
find the beginning of HA2, and pull out the 4 amino acid sites immediately preceding HA2. 
If HA2 is preceded immediately by amino acids R-X-K/R-R, then it is annotated
as having a furin cleavage motif. Otherwise, it is annotated as wild type. This produces 2
json files, one annotating a binary has a furin site or wild type site, the other coding 
the actual sequence at the cleavage sites. These can both be used in augur export so 
that the annotation and sequence show up in auspice as a color by.  
"""

import Bio
from Bio import SeqIO
import json

import argparse
parser = argparse.ArgumentParser()

parser.add_argument('--alignment', type=str, help='alignment file output by rule augur align')
parser.add_argument('--furin_site_motif', type=str, help='name of output json file that annotates tips as having a furin cleavage site or a wt cleavage site')
parser.add_argument('--cleavage_site_sequence', type=str, help='name of output json file that annotates the cleavage site sequence for each tip')


args = parser.parse_args()
alignment = args.alignment
furin_site_motif_json = args.furin_site_motif
cleavage_site_sequence_json = args.cleavage_site_sequence


def translate_nucleotide_to_aa(input_nucleotide_str):
    # convert to sequence object 
    sequence = Bio.Seq.Seq(input_nucleotide_str)
    # translate
    aa_sequence = str(sequence.translate())
    return(aa_sequence)


def return_ha2_start_position(aa_sequence):
    
    ha2_begin = "GLFG"
    
    # if .find does not match part of the string, it will return a value of -1 
    start_pos_ha2_aa = aa_sequence.find(ha2_begin)
    start_pos_ha2_nt = start_pos_ha2_aa*3
    return(start_pos_ha2_nt)


def output_furin_site_aa_sequence(ha2_nt_start, nt_sequence):

    """In the following function, we are collecting the 12 nucleotides that precede the start 
    of HA2 and translating them to yield the amino acid sequence for the 4 amino acids that 
    precede HA2. The reason we are doing it this way is to control for alignment 
    inconsistences from one build to the next, which sometimes insert gaps in different places
    in this region, resulting in inconsistent translations."""

    # find the 12 nucleotides that precede the start of HA2 that are not Ns 
    cleavage_site_nts = []
    current_site = ha2_nt_start
    # starting at the current base, go backwards until we accumulate 12 nts that aren't Ns
    while len(cleavage_site_nts) < 12: 
        nt_site = nt_sequence[current_site-1]
        if nt_site != "N":
            cleavage_site_nts.append(nt_site)
        current_site -= 1
            
    # reverse the list, translate it, and join into a string
    furin_site_nts = "".join(cleavage_site_nts[::-1])
    furin_site_sequence = Bio.Seq.Seq(furin_site_nts)
    furin_site = str(furin_site_sequence.translate())
    
    return(furin_site)


def output_furin_cleavage_site_jsons(alignment, output_json1, output_json2):
    
    output_dict_furin = {"nodes":{}}
    output_dict_seq = {"nodes":{}}
    
    with open(output_json1, "w") as outfile: 
        outfile.write("")
    with open(output_json2, "w") as outfile: 
        outfile.write("")


    for seq in SeqIO.parse(alignment, "fasta"):

        strain_name = seq.description

        # convert gaps to Ns to avoid translation errors
        nt_sequence = str(seq.seq).upper().replace("-","N")

        # translate to amino acids and find the start of HA2 in aa and nt coordinates
        aa_sequence = translate_nucleotide_to_aa(nt_sequence)

        start_pos_ha2_nt = return_ha2_start_position(aa_sequence)
        
        # if .find returns a null, it outputs a value of -1. This means that there is not adequate sequence data at the cleavage site, the start_pos will be a negative number. If it is positive, then the start of HA2 was present and we can annotate
        if start_pos_ha2_nt > 0: 

            # collect nts and translate the furin cleavage site sequence
            furin_site = output_furin_site_aa_sequence(start_pos_ha2_nt, nt_sequence)
        
            # if those 4 preceding amino acids have the pattern R-X-R/K-R, then it is cleavable 
            # by furin. Here, X is any amino acid (but not a gap), and the 3rd position can be 
            # K or R
            if furin_site[0] == "R" and furin_site[3] == "R" and (furin_site[2] == "R" or furin_site[2]=="K") and furin_site[1]!="X":
                furin_site_annotation = "present"
            else:
                furin_site_annotation = "absent"
        
        # if there was not sequence data at the cleavage site, annotate as Ns and missing data
        else:
            furin_site = "NNNN"
            furin_site_annotation = "missing data"
        
        output_dict_furin["nodes"][strain_name] = {"furin_cleavage_motif":furin_site_annotation}
        output_dict_seq["nodes"][strain_name] = {"cleavage_site_sequence":furin_site.replace("X","-")}
        
        
    f = open(output_json1, "w")
    json.dump(output_dict_furin, f)
    f.close()

    f = open(output_json2, "w")
    json.dump(output_dict_seq, f)
    f.close()


output_furin_cleavage_site_jsons(alignment, furin_site_motif_json, cleavage_site_sequence_json)



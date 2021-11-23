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
        sequence = str(seq.seq).upper().replace("-","N")

        # convert back to sequence object
        sequence = Bio.Seq.Seq(sequence) 

        # translate and find beginning of ha2
        aa = str(sequence.translate())
        ha2_begin = "GLFG"

        start_pos_ha2 = aa.find(ha2_begin)
        
        # define the furin site as the 4 positions prior to the start of HA2
        furin_site = aa[start_pos_ha2-4:start_pos_ha2]
		
		# if those 4 preceding amino acids have the pattern R-X-R/K-R, then it is cleavable 
		# by furin. Here, X is any amino acid (but not a gap), and the 3rd position can be 
		# K or R
        if furin_site[0] == "R" and furin_site[3] == "R" and (furin_site[2] == "R" or furin_site[2]=="K") and furin_site[1]!="X":
            furin_site_annotation = "present"
        else:
            furin_site_annotation = "absent"

        
        output_dict_furin["nodes"][strain_name] = {"furin_cleavage_motif":furin_site_annotation}
        output_dict_seq["nodes"][strain_name] = {"cleavage_site_sequence":furin_site.replace("X","-")}
        
    f = open(output_json1, "w")
    json.dump(output_dict_furin, f)
    f.close()

    f = open(output_json2, "w")
    json.dump(output_dict_seq, f)
    f.close()


output_furin_cleavage_site_jsons(alignment, furin_site_motif_json, cleavage_site_sequence_json)



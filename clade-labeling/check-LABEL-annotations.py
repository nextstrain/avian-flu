"""This will go through the output from LABEL and check clade assignments against pre-labeled clades.
Pre-labeled clade information is only available for sequences from the Influenza Research
Database. Still, this will check and print out any instances in which the clades do not 
match. It will also generate the {subtype}-clades.tsv file which is used in the Snakemake
pipeline. In general, it seems like clade labels match quite well."""


import argparse
parser = argparse.ArgumentParser()

parser.add_argument('--label_output', type=str, help='path to the output file from label that ends in "_final.tab"')
parser.add_argument('--output', type=str, help='output file name for clades file')

args = parser.parse_args()
label_output = args.label_output
output_filename = args.output


with open(output_filename, "w") as outfile: 
    outfile.write("name\tclade\n")

print("\n\nPrinting strains for which LABEL annotation does not match what was annotated in IRD")
print("strain name\tIRD clade\tLABEL clade")

with open(label_output, "r") as infile: 
    for line in infile: 
        if "VIRUS STRAIN" not in line:
            strain_name = line.split("\t")[0].strip().replace(" ","")
            clade_assignment = line.split("\t")[1].replace("_","-").replace(" ","").strip().replace("UNRECOGNIZABLE","?")
        
            with open(output_filename, "a") as outfile: 
                outfile.write(strain_name + "\t" + clade_assignment + "\n")

"""This script reads in the clades file and adds an H5 clade column to the metadata tsv. 
This creates a new metadata file that is then used in downstream snakemake fules for the 
H5N1 or H5NX builds. The clades file is a tab-delimited file that contains the strain
name of each H5N1 or H5NX strain with its assigned clade. Clade assignments for now 
are performed with LABEL, a CDC developed software for assigning H5 clades. This is used in 
the avian-flu pipeline in rule add_h5_clade. 

usage: add-clades.py --clades clade-file.tsv --metadata input-metadata.tsv --output output-file-name"""


import argparse
parser = argparse.ArgumentParser()

parser.add_argument('--clades', type=str, help='a tab-delimited text file containing H5 HA clade assignments by LABEL')
parser.add_argument('--metadata', type=str, help='metadata file output from augur parse')
parser.add_argument('--output', type=str, help='output file name for metadata with annotated clades')


args = parser.parse_args()
metadata_infile = args.metadata
clades_file = args.clades
metadata_outfile = args.output


def read_in_clades_file(clades_file):
	
	clade_assignments = {}
	with open(clades_file,"r") as infile: 
		for line in infile: 
			strain_name = line.split("\t")[0]
			clade = line.strip().split("\t")[1]
		
			clade_assignments[strain_name] = clade
			
	return(clade_assignments)
		

		
def annotate_metadata_file(metadata_infile, metadata_outfile, clade_assignments):
	known = 0
	unknown_clades = 0
	
	with open(metadata_outfile, "w") as outfile: 
		outfile.write("")
	
	with open(metadata_infile, "r") as infile:
		linecount = 0
		for line in infile:
			linecount += 1
			if linecount == 1:
				new_line = line.strip() + "\th5_label_clade\n"
			else:
				strain_name = line.split("\t")[0]
				if strain_name in clade_assignments:
					clade = clade_assignments[strain_name]
					known += 1
				else:
					clade = "?"
					unknown_clades += 1
					print("unknown clade for ", strain_name)
				new_line = line.strip() + "\t" + clade + "\n"
		
			with open(metadata_outfile, "a") as outfile:
				outfile.write(new_line)
	print("was able to assign clades to", str(known), "strains")
	print("unable to find clades for", str(unknown_clades), "strains")

				

# run 
clade_assignments = read_in_clades_file(clades_file)
annotate_metadata_file(metadata_infile, metadata_outfile, clade_assignments)

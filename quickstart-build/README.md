# Avian flu quickstart build

This is a simplified version of the avian-flu builds that are hosted on Nextstrain. The public avian influenza builds have become somewhat complicated, and include features that are likely unnecessary for most users. For example, we include calls to download data from our databases, and use a separate tool for clade annotation. We have received some requests from users for a simpler build that will immediately work, given a set of example data. This is that build. The intention for this simplified build is to provide a ready to use pipeline for users to use as a baseline to develop their own custom builds, without having to weed through excessive Snakemake documentation or delete parts of the pipeline that are useful only to the Nextstrain team. This build includes the following alterations:

1. **A simplified and heavily commented Snakefile** 
This Snakefile retains the wildcard structure in the full avian-flu build, but removes features that are unnecessary for most users. I've also added in comments throughout the Snakefile that break down the parts of the Snakefile, annotates where inputs and outputs are specified, and explains a few spots that can be edited to customize builds. The main feature that has been removed from this build is clade annotation. Many users will start their builds will pre-annotated sequences, and because the clade annotation significantly complicates the pipeline, it is removed here. 

2. **An input dataset of H5N1 sequences from GenBank**
The pipeline as currently included in this Snakefile will perform the pipeline beginning with the fasta files in the `example_data` folder. To run with your own custom data, simply put your fasta files into this folder instead.  

3. **Complete output folders**
These builds are small and based completely off data from Genbank, allowing us to include all of the intermediate files that are generated as part of the pipeline. It is my hope that this allows users to compare steps in the pipeline where they may be encountering errors. 


## To make your own, custom build
The [Nextstrain docs](https://docs.nextstrain.org/en/latest/index.html) are a fantastic resource for getting started with the Nextstrain pipeline, and include some [great tutorials](https://docs.nextstrain.org/en/latest/install.html) to get you started. Once you have those tools installed, do the following: 

1. clone the `avian-flu` repo: 

`git clone https://github.com/nextstrain/avian-flu.git`

2. Navigate to the quickstart: 

`cd avian-flu/quickstart-build/`

3. Open the `Snakefile` in a text editor, and read through the format and comments. 

4. Test that the build works with the example data:

`conda activate nextstrain`

`snakemake -p --cores 1`

5. Replace the example data with your own data, and rerun: 

`snakemake -p --cores 1`


#### A note on metadata customization
This pipeline is extremely flexible, and can be used to build trees based on any set of metadata you have available for your input data. In the example data, we've included several columns of metadata, but these are not all necessary for you to run a custom build, and you can easily add more. For example, if your sequences have been annotated with clades, add that clade information into the sequence headers and into the `augur parse` rule, and it will be included in the build. Strain name is a required field, and collection date will be necessary for inferring time trees. If you lack information on geography, you can still build a tree, but the map functionality will not work, and you will need to change the subsampling rules. To alter colorby options, edit/add/remove fields in `config/auspice_config_{subtype}.json`.

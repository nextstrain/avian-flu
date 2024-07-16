We gratefully acknowledge the authors, originating and submitting laboratories of the genetic sequences and metadata for sharing their work. Please note that although data generators have generously shared data in an open fashion, that does not mean there should be free license to publish on this data. Data generators should be cited where possible and collaborations should be sought in some circumstances. Please try to avoid scooping someone else's work. Reach out if uncertain.

Genomic data from the ongoing outbreak of H5N1 in cattle in the US was shared by the [National Veterinary Services Laboratories (NVSL)](https://www.aphis.usda.gov/labs/about-nvsl) of the [Animal and Plant Health Inspection Service (APHIS)](https://www.aphis.usda.gov/) of the U.S. Department of Agriculture (USDA) in an open fashion to NCBI GenBank (consensus genomes and complete metadata) and to the SRA (raw reads with redacted metadata) in [BioProject PRJNA1102327](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1102327). Other groups have contributed sequence data here, but the majority of viral genomes have been shared by the USDA. The Andersen Lab has assembled raw reads from this SRA BioProject and publicly shared consensus genomes to [GitHub](https://github.com/andersen-lab/avian-influenza). We thank the USDA for genomic data sharing and the Andersen Lab for sharing assembled consensus genomes.

In this analysis, we've curated data from NCBI GenBank and merged this data with SRA data via the Andersen Lab GitHub repository. Curated sequence data is available as:
 - [data.nextstrain.org/files/workflows/avian-flu/h5n1/ha/sequences.fasta.zst](https://data.nextstrain.org/files/workflows/avian-flu/h5n1/ha/sequences.fasta.zst)
 - [data.nextstrain.org/files/workflows/avian-flu/h5n1/na/sequences.fasta.zst](https://data.nextstrain.org/files/workflows/avian-flu/h5n1/na/sequences.fasta.zst)
 - etc...

and metadata is available as:
 - [data.nextstrain.org/files/workflows/avian-flu/h5n1/metadata.tsv.zst](https://data.nextstrain.org/files/workflows/avian-flu/h5n1/metadata.tsv.zst)


These data are [updated daily](https://github.com/nextstrain/avian-flu/actions/workflows/ingest-ncbi.yaml) pulling from GenBank and GitHub. Data source as GenBank vs SRA-via-Andersen-Lab is included in this metadata and is available as a [coloring to this page](?c=data_source). Importantly, SRA-derived genomes only have "2024-XX-XX" as collection date and "USA" as collection location. In this analysis, we've inferred collection date and collection location for these samples along with confidence in date and location.

In addition to this cattle outbreak specific view, we have broader views of H5N1 evolution available as:
 - [nextstrain.org/avian-flu/h5n1/ha/2y](https://nextstrain.org/avian-flu/h5n1/ha/2y)
 - [nextstrain.org/avian-flu/h5n1/na/2y](https://nextstrain.org/avian-flu/h5n1/na/2y)
 - etc...

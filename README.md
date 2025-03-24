# Long read data standardization SOP for downstream analysis
## Mike Nalls, Jon Moller, Rylee Genner, Melissa Meredith, and Kimberley Billingsley
This repository summarizes the CARD Applied Neurogenomics Group's long read sequencing data standardization SOP for downstream analyses like QTLs.
# File formats and examples
## Metadata
## Omics
## Genetic data
## Genetic maps
# Preprocessing scripts and usage
## Metadata
## Omics
## Genetic data
Genetic data and maps are generated together in two steps. First, input SV or SNV VCF variant files from Napu are normalized to biallelic variants and then converted to a CSV file with the necessary fields using bcftools query. Sample data is included in the repository as example_sv_input.vcf.gz, example_snv_input.vcf.gz, example_sv_preprocess.csv, and example_snv_preprocess.csv. Based on preliminary testing and available NIH HPC resources, we recommend parallelizing genetic map and data generation by chromosome. Tests on the HBCC cohort SNV VCF indicated chromosome 1 genetic map/data processing took 15 minutes on a 32GB RAM/32 CPU allocation.
```
Usage: scripts/vcf_preprocess_for_genetic_data_map.sh -v variant_type -c chromosome -i input.vcf(.gz) -o output.csv
	-v Variant type - structural variants (SV) or single nucleotide variants (SNV)
	-c Chromosome to subset (optional)
	-i Input VCF file
	-o Output CSV file for genetic data/map generation with helper Python script
```
```
usage: make_genetic_data_and_map.py [-h] -s {SV,SNV} -i INPUT_FILE -o OUTPUT_PREFIX [-c CHROMOSOME]

Convert preprocessed VCF input in CSV data table to genetic map and data matrix CSV for downstream analysis (e.g., QTLs).

optional arguments:
  -h, --help            show this help message and exit
  -s {SV,SNV}, --file_type {SV,SNV}
                        Specify the file type: SVs for Structural Variants or SNVs for Single Nucleotide Variants.
  -i INPUT_FILE, --input_file INPUT_FILE
                        Path to the input preprocessed CSV file for parsing (from vcf_preprocess_for_genetic_data_map.sh).
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        Prefix for the output genetic data and genetic map files.
  -c CHROMOSOME, --chromosome CHROMOSOME
                        Chromosome(s) to include (optional).
```
## Genetic maps
Genetic maps are made together with genetic data above. Sample output genetic maps are included in the repository as example_sv_genetic_map.csv and example_snv_genetic_map.csv.

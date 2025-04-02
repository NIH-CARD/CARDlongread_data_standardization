# Long read data standardization SOP for downstream analysis
## Mike Nalls, Jon Moller, Rylee Genner, Melissa Meredith, and Kimberley Billingsley
This repository summarizes the CARD Applied Neurogenomics Group's long read sequencing data standardization SOP for downstream analyses like QTLs. It includes both examples of file formats for each downstream analysis input data type and preprocessing scripts to convert raw data files (e.g., VCFs, methylation BEDs) into tables suitable for downstream analysis.
# File formats and examples
## Metadata
Metadata files contain individual level data and potentially protected health information (PHI) so be cautious with sharing. In general these are sample level metrics not stratified by haplotype. 

Mandatory columns are as follows:
SAMPLE = sample ID in alphanumeric string format.
MALE = binary indicator confirming XY sex chromosomes
STATUS = disease of interest, controls for the population are always designated by the string value “control”.
CHEMISTRY = an alphanumeric code for the chemistry used such as “R9” or “R10”.
Any additional columns are considered covariates generally.

An example of one such file is below:
SAMPLE,MALE,STATUS,CHEMISTRY,PC1,PC2,PC3,PC4,PC5
indiv1,1,control ,R9,0.001,-0.002,0.005.0.006,-0.1
indiv2,1,PD ,R9,0.01,-0.003,0.006.0.0068,-0.4,
indiv3,0,AD ,R10,0.3,-0.004,0.0065.0.08,0.002

## Omics
## Genetic data
## Genetic maps
# Preprocessing scripts and usage
## Metadata
## Omics
We have so far developed a script to convert methylation BED files from Napu into methylation data and map files as described above. This script takes methylation BED files for each haplotype as input and outputs methylation data and map files as CSV output. It is written to process CpG island (CGI), gene body (GB), and promoter (PROM) regions. It performs a full outer join on each haplotype (union of methylation regions) and fills in missing values in the respective haplotype as NA. Sample methylation data and map files are provided as example_methylation_data.csv and example_methylation_map.csv.
```
usage: make_methylation_data_and_map.py [-h] -t {CGI,GB,PROM} -h1 HAPLOTYPE_1 -h2 HAPLOTYPE_2 -p REGION_PREFIX -o OUTPUT_PREFIX

Convert methylation BED files to methylation map and data matrix CSVs for downstream analysis (e.g., QTLs).

optional arguments:
  -h, --help            show this help message and exit
  -t {CGI,GB,PROM}, --region_type {CGI,GB,PROM}
                        Specifying region type for parsing and joining. Currently supported are CpG islands (CGI), gene bodies (GB), and promoters (PROM).
  -h1 HAPLOTYPE_1, --haplotype_1 HAPLOTYPE_1
                        Methylation BED file for the first haplotype (H1).
  -h2 HAPLOTYPE_2, --haplotype_2 HAPLOTYPE_2
                        Methylation BED file for the second haplotype (H2).
  -p REGION_PREFIX, --region_prefix REGION_PREFIX
                        Prefix for the output methylation region names (e.g., CGIs, gene bodies, promoters).
  -o OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                        Prefix for the output methylation data and genetic map files.
```
## Genetic data
Genetic data and maps are generated together in two steps. First, input SV or SNV VCF variant files from Napu are normalized to biallelic variants, subset by chromosome, and then converted to a CSV file with the necessary fields using bcftools query (`vcf_preprocess_for_genetic_data_map.sh`). Next, preprocessed CSV files are converted to per sample, per haplotype genetic data matrix and per variant genetic map files using `make_genetic_data_and_map.py`. Sample preprocessed data is included in the repository as example_sv_preprocess.csv and example_snv_preprocess.csv. Based on preliminary testing and available NIH HPC resources, we recommend parallelizing genetic map and data generation by chromosome. Tests on the HBCC cohort SNV VCF indicated chromosome 1 genetic map/data processing took 15 minutes on a 32GB RAM/32 CPU allocation. Example output genetic data are also provided as example_sv_data.csv and example_snv_data.csv.
```
Usage: vcf_preprocess_for_genetic_data_map.sh -v variant_type -c chromosome -i input.vcf(.gz) -o output.csv
	-v Variant type - structural variants (SV) or single nucleotide variants (SNV)
	-c Chromosome to subset (optional)
	-i Input VCF file
	-o Output CSV file for genetic data/map generation with helper Python script
```
Below is an example of preprocessed output from the above script:
```
ID,CHROM,POS,END,REF,ALT,HBCC_81935_FTX,HBCC_81931_FTX,HBCC_81929_FTX,HBCC_81934_FTX
chr1_10611_C_G,chr1,10611,10611,C,G,0/1,./.,./.,0/1
chr1_10622_T_G,chr1,10622,10622,T,G,./.,1/0,1/0,./.
chr1_10623_T_C,chr1,10623,10623,T,C,./.,1/1,1/1,0/1
chr1_10627_A_AG,chr1,10627,10627,A,AG,0/0,1/0,1/0,./.
chr1_10815_T_TC,chr1,10815,10815,T,TC,./.,./.,./.,1/0
chr1_10927_A_G,chr1,10927,10927,A,G,./.,1/0,1/0,./.
chr1_10936_G_C,chr1,10936,10936,G,C,./.,./.,./.,./.
chr1_10968_A_AG,chr1,10968,10968,A,AG,./.,./.,1/0,0/0
chr1_11002_A_C,chr1,11002,11002,A,C,./.,1/0,1/0,./.
```
```
usage: make_genetic_data_and_map.py [-h] -s {SV,SNV} -i INPUT_FILE -o OUTPUT_PREFIX [-c CHROMOSOME]

Convert preprocessed VCF input in CSV data table (from vcf_preprocess_for_genetic_data_map.sh) to genetic map and data matrix CSVs for downstream analysis (e.g., QTLs).

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
Genetic maps are made together with genetic data above. Sample output genetic maps corresponding to above genetic data are included in the repository as example_sv_genetic_map.csv and example_snv_genetic_map.csv.

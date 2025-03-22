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
```
usage: make_genetic_data.py [-h] -s {SV,SNV} -i INPUT_FILE -o OUTPUT_FILE [-c CHROMOSOME] [-t THREADS]

Convert input VCF file to genetic data matrix for downstream analysis (e.g., QTLs).

optional arguments:
  -h, --help            show this help message and exit
  -s {SV,SNV}, --file_type {SV,SNV}
                        Specify the file type: SVs for Structural Variants or SNVs for Single Nucleotide Variants.
  -i INPUT_FILE, --input_file INPUT_FILE
                        Path to the input VCF file for parsing.
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Path to the output .csv file where the final converted data will be saved.
  -c CHROMOSOME, --chromosome CHROMOSOME
                        Chromosome(s) to include (optional)
  -t THREADS, --threads THREADS
                        Threads for VCF decompression (optional)
```
## Genetic maps
```
usage: make_genetic_map.py [-h] -s {SV,SNV} -i INPUT_FILE -o OUTPUT_FILE [-c CHROMOSOME] [-t THREADS]

Convert input VCF file to genetic map for downstream analysis (e.g., QTLs).

optional arguments:
  -h, --help            show this help message and exit
  -s {SV,SNV}, --file_type {SV,SNV}
                        Specify the file type: SVs for Structural Variants or SNVs for Single Nucleotide Variants.
  -i INPUT_FILE, --input_file INPUT_FILE
                        Path to the input VCF file for parsing.
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Path to the output .csv file where the final converted data will be saved.
  -c CHROMOSOME, --chromosome CHROMOSOME
                        Chromosome(s) to include (optional)
  -t THREADS, --threads THREADS
                        Threads for VCF decompression (optional)
```

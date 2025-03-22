#!/usr/bin/python 

import argparse
import pysam
from pysam import VariantFile
import pandas as pd
import numpy as np
import os
import csv

"""
    convert SV/SNV variant VCF file into genetic map CSV file following standards set in March 2025.
    
    Usage: 
    python make_genetic_map.py -s [SV | SNV] -i input.vcf.gz -o output.csv

    output: 
    A csv file with following columns: variant name (NAME), CHROM, START, STOP, REF, ALT
    variant name (SNV) = SNV_CHROM_START_REF_ALT
    variant name (SV) = SV_CHROM_SVTYPE_START_STOP

    Rylee Genner 03/2025
    Jon Moller edits 03/2025
"""

def write_output_header(output_file):
    """Print opening header of CSV"""
    with open(output_file, "w") as outfile:
        outfile.write("NAME,CHROM,START,STOP,REF,ALT\n")

# routine to process SV records appropriately		
def sv_rec_process(variant_rec,alt_index):
	# if SV is insertion, stop and start positions are the same
    if (variant_rec.info["SVTYPE"]=="INS"):
       stop=variant_rec.pos   
	# if SV is deletion, stop position is start position plus absolute value of SV length minus 1
    elif (variant_rec.info["SVTYPE"]=="DEL"):
       stop=variant_rec.pos+abs(variant_rec.info["SVLEN"])-1
    else:
       stop=variant_rec.pos
	# sv name is SV_CHROM_SVTYPE_START_STOP   
    sv_name="SV_" + variant_rec.chrom + "_" + variant_rec.info["SVTYPE"] + "_" + str(variant_rec.pos) + "_" + str(stop)
    sv_chrom=variant_rec.chrom
    sv_start=variant_rec.pos
    sv_stop=stop
    sv_ref=variant_rec.ref
    # just use single ALT allele (rec.alts[alt_index]) with option to specify which ALT allele
    sv_alt=variant_rec.alts[alt_index]
    # output row is NAME,CHROM,START,STOP,REF,ALT
    sv_rec_row=(sv_name,sv_chrom,sv_start,sv_stop,sv_ref,sv_alt)
    return(sv_rec_row)

# routine to process SNV records appropriately
def snv_rec_process(variant_rec,alt_index):
    # sv name is SNV_CHROM_START_REF_ALT  
    snv_name="SNV_" + variant_rec.chrom + "_" + str(variant_rec.pos) + "_" + variant_rec.ref + "_" + variant_rec.alts[alt_index]
    snv_chrom=variant_rec.chrom
    snv_start=variant_rec.pos
    # just use single ALT allele (rec.alts[alt_index]) with option to specify which ALT allele
    snv_alt=variant_rec.alts[alt_index]
    snv_ref=variant_rec.ref
    # snv stop depends upon whether insertion, deletion, or substitution
    # insertion: stop and start are both record position
    if (len(snv_alt) > len(snv_ref)):
        snv_stop=variant_rec.pos
    # deletion: start is record position plus alternate allele length minus 1; stop is record position plus reference allele length minus 1
    elif (len(snv_alt) < len(snv_ref)):
        snv_start=variant_rec.pos+len(snv_alt)-1
        snv_stop=variant_rec.pos+len(snv_ref)-1
    # substitution: start and stop are both record position
    elif (len(snv_alt) == len(snv_ref)):
        snv_stop=variant_rec.pos
    # output row is NAME,CHROM,START,STOP,REF,ALT
    snv_rec_row=(snv_name,snv_chrom,snv_start,snv_stop,snv_ref,snv_alt)
    return(snv_rec_row)
    
# routine to process whole vcf
def vcf_process_loop(file_type,chromosome,variant_file,output_file):
    # print(file_type)
    # procedure for structural variants
    if (file_type == "SV"):
        for rec in variant_file.fetch(chromosome):
            # get number of alternate alleles per record (per position)
            num_sv_alts=len(rec.alts)
            # print(num_sv_alts)
            # loop through alternate alleles by index (range 0 to number of alternate alleles)
            for i in range(0,num_sv_alts):
                # rec row for SVs
                rec_row=sv_rec_process(rec,i)
                # output rec_row to csv
                with open(output_file, "a") as outfile:
                    csv_writer = csv.writer(outfile)
                    csv_writer.writerow(rec_row)
    # procedure for single nucleotide variants
    elif (file_type == "SNV"):
        for rec in variant_file.fetch(chromosome):
            # get number of alternate alleles per record (per position)
            num_snv_alts=len(rec.alts)
            # print(num_snv_alts)
            # loop through alternate alleles by index (0 to number of alternate alleles)
            for i in range(0,num_snv_alts):
                # rec row for SNVs
                rec_row=snv_rec_process(rec,i)
                # output rec_row to csv
                with open(output_file, "a") as outfile:
                    csv_writer = csv.writer(outfile)
                    csv_writer.writerow(rec_row)

parser = argparse.ArgumentParser(description="Convert input VCF file to genetic map for downstream analysis (e.g., QTLs).")

# Add argument for file type (SV or SNV) first
parser.add_argument(
    "-s", "--file_type",
    required=True,
    choices=["SV", "SNV"],
    help="Specify the file type: SVs for Structural Variants or SNVs for Single Nucleotide Variants."
)

# Add argument for input file
parser.add_argument(
    "-i", "--input_file",
    required=True,
    help="Path to the input VCF file for parsing."
)

# Add argument for final output tsv file
parser.add_argument(
    "-o", "--output_file",
    required=True,
    help="Path to the output .csv file where the final converted data will be saved."
)

# Add argument for chromosome 
parser.add_argument(
    "-c", "--chromosome",
    required=False,
    help="Chromosome(s) to include (optional)"
)

# Add argument for threads
parser.add_argument(
    "-t", "--threads",
    required=False,
    default=1,
    type=int,
    help="Threads for VCF decompression (optional)"
)

# Parse the arguments
args = parser.parse_args()

# import VCF file
input_variant_file=VariantFile(args.input_file,threads=args.threads)

# make output header
write_output_header(args.output_file)

# loop through VCF records to make map
vcf_process_loop(args.file_type,args.chromosome,input_variant_file,args.output_file)

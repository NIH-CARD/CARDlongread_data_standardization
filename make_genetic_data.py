#!/usr/bin/python 

import argparse
import pysam
from pysam import VariantFile
import pandas as pd
import numpy as np
import os
import csv

"""
    convert SV/SNV variant VCF file into genetic data (variant presence/absence/missingness per sample, per haplotype) CSV file following standards set in March 2025.
    
    Usage: 
    python make_genetic_data.py -s [SVs | SNVs] -i input.vcf.gz -o output.csv

    output: 
    A csv file with following columns: sample (SAMPLE), haplotype (HAPLOTYPE), and each variant in input file
    variant name (SNV) = SNV_CHROM_START_REF_ALT
    variant name (SV) = SV_CHROM_SVTYPE_START_STOP
    Each row contains sample (e.g., HBCC_81925_FTX), haplotype (H1 or H2), and presence (1), absence (0), or missingness (NA) of given variant
 
    Rylee Genner 03/2025
    Jon Moller edits 03/2025
"""

# initialize haplotype matrix with sample and haplotype names
def initialize_matrix(variant_file):
    # prepare first two columns of output matrix, SAMPLE and HAPLOTYPE
    # sample list
    sample_list=list(variant_file.header.samples)
    # sample number
    number_of_samples=len(sample_list)
    # create output matrix
    output_matrix_df=pd.DataFrame(dtype='Int64')
    # initialize SAMPLE column with sample list repeated once
    output_matrix_df['SAMPLE']=sample_list+sample_list
    # initialize HAPLOTYPE column with H1 repeated for each sample and H2 repeated for each sample
    haplotype_list=['H1'] * number_of_samples + ['H2'] * number_of_samples
    output_matrix_df['HAPLOTYPE']=haplotype_list
    return(output_matrix_df)
    
# routine to process whole vcf into output genetic data matrix
def make_genetic_haplotype_matrix(file_type,chromosome,variant_file,output_matrix_df):
    # procedure for structural variants
    if (file_type == "SV"):
        for rec in variant_file.fetch(chromosome):
            # exclude multiallelic sites
            if (len(rec.alts) == 1):
                # name SV the same way I would in the genetic map
                # first get stop position
                # if SV is insertion, stop and start positions are the same
                if (rec.info["SVTYPE"]=="INS"):
                   stop=rec.pos   
            	# if SV is deletion, stop position is start position plus absolute value of SV length minus 1
                elif (rec.info["SVTYPE"]=="DEL"):
                   stop=rec.pos+abs(rec.info["SVLEN"])-1
                else:
                   stop=rec.pos
                # now concatenate elements into full SV name (SV_CHROM_TYPE_START_STOP)
                sv_name="SV_" + rec.chrom + "_" + rec.info["SVTYPE"] + "_" + str(rec.pos) + "_" + str(stop)
                # get all sample haplotype H1 genotypes of variant
                gts_h1 = [s['GT'][0] for s in rec.samples.values()]
                # get all sample haplotype H2 genotypes of variant
                gts_h2 = [s['GT'][1] for s in rec.samples.values()]
                # make genotype column (h1 first, then h2)
                gt_column = gts_h1 + gts_h2
                # append column to output matrix with sv_name as column name
                output_matrix_df[sv_name]=gt_column
                # create genotype column dataframe (h1 first, then h2)
                # gt_column_df=pd.DataFrame(gts_h1 + gts_h2,columns=[sv_name])
                # output_matrix_df_concat=pd.concat([output_matrix_df,gt_column_df])
                # force to be integer
                output_matrix_df_concat[sv_name].astype('Int64')
    # procedure for single nucleotide variants
    elif (file_type == "SNV"):
        for rec in variant_file.fetch(chromosome):
            # exclude multiallelic sites
            if (len(rec.alts) == 1):
                # name SNV the same way I would in the genetic map
                snv_name="SNV_" + rec.chrom + "_" + str(rec.pos) + "_" + rec.ref + "_" + rec.alts[0]
                # get haplotype H1 genotypes of variant
                gts_h1 = [s['GT'][0] for s in rec.samples.values()]
                # get haplotype H2 genotypes of variant
                gts_h2 = [s['GT'][1] for s in rec.samples.values()]
                # make genotype column (h1 first, then h2)
                gt_column = gts_h1 + gts_h2
                # append column to output matrix with snv_name as column name
                output_matrix_df[snv_name]=gt_column
                # create genotype column dataframe (h1 first, then h2)
                # gt_column_df=pd.DataFrame(gts_h1 + gts_h2,columns=[snv_name])
                # output_matrix_df=pd.concat([output_matrix_df,gt_column_df])
                # force to be integer
                output_matrix_df[snv_name].astype('Int64')
    return(output_matrix_df)

parser = argparse.ArgumentParser(description="Convert input VCF file to genetic data matrix for downstream analysis (e.g., QTLs).")

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

# initialize matrix
genetic_data_matrix=initialize_matrix(input_variant_file)

# loop through VCF records to make map
make_genetic_haplotype_matrix(args.file_type,args.chromosome,input_variant_file,genetic_data_matrix)

# export to csv
# represent absent values as NA
genetic_data_matrix.to_csv(args.output_file,na_rep='NA',index=False)
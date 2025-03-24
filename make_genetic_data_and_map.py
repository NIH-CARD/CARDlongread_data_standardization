#!/usr/bin/python 

# make_genetic_data_and_map.py
# script to convert bcftools query output from 

import argparse
import pandas as pd
import numpy as np
import os
import re

# routine to name variants
def name_genetic_variants(variant_type,input_data_frame):
    # make variant names depending upon variant type
    # if structural variants, SV_CHROM_SVTYPE_START_STOP
    if (variant_type == "SV"):
        variant_name="SV_"+input_data_frame['CHROM']+"_"+input_data_frame['SVTYPE']+"_"+input_data_frame['POS'].astype(str)+"_"+input_data_frame['END'].astype(str)
    # if single nucleotide variants, SNV_CHROM_START_STOP_REF_ALT
    elif (variant_type == "SNV"):
        variant_name="SNV_"+input_data_frame['CHROM']+"_"+input_data_frame['POS'].astype(str)+"_"+input_data_frame['REF']+"_"+input_data_frame['ALT']
    return(variant_name)
    
# routine to process whole imported data frame into genetic map
def make_genetic_map(variant_type,input_data_frame):
    # initialize output genetic map data frame
    genetic_map_df=pd.DataFrame(dtype='Int64',columns=['NAME','CHROM','START','STOP','REF','ALT'])
    # name variants depending upon variant type
    genetic_map_df['NAME']=name_genetic_variants(variant_type,input_data_frame)
    # select columns of interest
    genetic_map_df['CHROM']=input_data_frame['CHROM']
    genetic_map_df['START']=input_data_frame['POS']
    genetic_map_df['STOP']=input_data_frame['END']
    genetic_map_df['REF']=input_data_frame['REF']
    genetic_map_df['ALT']=input_data_frame['ALT']
    # return genetic map data frame
    return(genetic_map_df)
    
# routine to process whole imported data frame into output genetic data matrix
def make_genetic_haplotype_matrix(variant_type,input_data_frame):
    # get sample names - exclude first six columns (ID, CHROM, POS, END, REF, ALT) - elements 0 through 5
    if (variant_type=='SNV'):
        sample_names=input_data_frame.columns.tolist()[6:]
    elif (variant_type=='SV'):
        sample_names=input_data_frame.columns.tolist()[8:]
    # prepare column names for haplotype matrix to be transposed and appended with sample name and haplotype info
    haplotype_sample_names=sum(([ [i+"_H1", i+"_H2"] for i in sample_names ]),[])
    # prepare initial untransposed haplotype matrix
    transposed_genetic_haplotype_matrix_df = pd.DataFrame(dtype='Int64',columns=haplotype_sample_names)
    # split genotypes per sample
    for i in sample_names:
        # split by EITHER "/" or "|" (e.g., 0/1 or 0|1)
        transposed_genetic_haplotype_matrix_df[[i+"_H1",i+"_H2"]]=input_data_frame[i].str.split(r"\/|\|",expand=True)
    # transpose data frame into final output
    genetic_haplotype_matrix_df=transposed_genetic_haplotype_matrix_df.transpose()
    # remove original data frame to save memory
    del(transposed_genetic_haplotype_matrix_df)
    # name variants depending upon variant type AFTER transpose; make column names match these variant names
    genetic_haplotype_matrix_df.columns=name_genetic_variants(variant_type,input_data_frame)
    # add sample and haplotype columns
    haplotype_list = ['H1','H2'] * len(sample_names)
    genetic_haplotype_matrix_df.insert(loc=0,column='HAPLOTYPE',value=haplotype_list)
    # sample list - two of each sample name given two haplotypes
    doubled_sample_name_list=sum(([i] * 2 for i in sample_names),[])
    genetic_haplotype_matrix_df.insert(loc=0,column='SAMPLE',value=doubled_sample_name_list)
    # replace periods with NAs if necessary
    genetic_haplotype_matrix_df.replace('.', np.nan, inplace=True)
    # return genetic haplotype matrix data frame
    return(genetic_haplotype_matrix_df)
    
parser = argparse.ArgumentParser(description="Convert preprocessed VCF input in CSV data table (from vcf_preprocess_for_genetic_data_map.sh) to genetic map and data matrix CSVs for downstream analysis (e.g., QTLs).")

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
    help="Path to the input preprocessed CSV file for parsing (from vcf_preprocess_for_genetic_data_map.sh)."
)

# Add argument for final output tsv file
parser.add_argument(
    "-o", "--output_prefix",
    required=True,
    help="Prefix for the output genetic data and genetic map files."
)

# Add argument for chromosome 
parser.add_argument(
    "-c", "--chromosome",
    required=False,
    help="Chromosome(s) to include (optional)."
)

# Parse the arguments
args = parser.parse_args() 
# import VCF preprocessed into CSV with bcftools query   
vcf_preprocess_df=pd.read_csv(args.input_file)
# if chromosome specified
if (args.chromosome is not None):
    vcf_preprocess_df=vcf_preprocess_df[vcf_preprocess_df['CHROM']==args.chromosome]
# make genetic map
output_genetic_map_df=make_genetic_map(args.file_type,vcf_preprocess_df)
# make genetic data
output_genetic_data_df=make_genetic_haplotype_matrix(args.file_type,vcf_preprocess_df)
# export genetic map as CSV
# 1e6 lines at a time
output_genetic_map_df.to_csv(args.output_prefix+"_genetic_map.csv",index=False,header=True,na_rep="NA",chunksize=1000000)
# export genetic data as CSV
# 1e3 lines at a time
output_genetic_data_df.to_csv(args.output_prefix+"_genetic_data.csv",index=False,header=True,na_rep="NA",chunksize=1000)
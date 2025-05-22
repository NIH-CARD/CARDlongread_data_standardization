#!/usr/bin/python 

# make_methylation_data_and_map.py
# script to convert methylation bed per haplotype output (for CpG islands, gene bodies, or promoters) into methylation data and map for QTLs
import argparse
import pandas as pd
import numpy as np
import os
import re

# routine to filter initial methylation bed by missing information rate
def filter_regions_by_missing_info_rate(region_type,input_data_frame,missing_info_rate):
    if (region_type=="CGI"):
        # calculate missing information rate for each row
        df_missing_info_rates=input_data_frame.iloc[:,10:].isna().sum(axis=1)/input_data_frame.iloc[:,10:].shape[1]
        # subset data frame based on missing information rate
        filtered_data_frame=input_data_frame[df_missing_info_rates < missing_info_rate]
    else:
        # calculate missing information rate for each row
        df_missing_info_rates=input_data_frame.iloc[:,6:].isna().sum(axis=1)/input_data_frame.iloc[:,6:].shape[1]
        # subset data frame based on missing information rate
        filtered_data_frame=input_data_frame[df_missing_info_rates < missing_info_rate]
    # return filtered bed lacking regions with missing methylation information for proportion of samples above missing information rate
    return(filtered_data_frame)

# routine to rename methylation regions for each haplotype uniquely
def rename_methylation_regions(input_data_frame,region_prefix):
    # create copy of data frame so original is left unchanged
    renamed_data_frame=input_data_frame.copy()
    # make new region name in following format: "PREFIX_CHROM_START_STOP"
    new_region_name=region_prefix + "_" + input_data_frame['chrom'] + "_" + input_data_frame['start'].astype(str) + "_" + input_data_frame['end'].astype(str)
    # replace name column with standardized region name
    renamed_data_frame['name']=new_region_name
    return(renamed_data_frame)

# routine to relabel and outer join per haplotype bed files 
def join_relabeled_haplotype_beds(region_type,hap1_data_frame,hap2_data_frame):
    # create copies of data frames so original is left unchanged
    hap1_data_frame_corrected = hap1_data_frame.copy()
    hap2_data_frame_corrected = hap2_data_frame.copy()
    # correct names of columns (e.g, remove avgMod_ prefix from column names)
    hap1_data_frame_corrected.columns = hap1_data_frame_corrected.columns.str.removeprefix("avgMod_") 
    hap2_data_frame_corrected.columns = hap2_data_frame_corrected.columns.str.removeprefix("avgMod_")
    # correct names of columns in alternate case (e.g., remove _GRCh38_1_modFraction suffix from column names)
    hap1_data_frame_corrected.columns = hap1_data_frame_corrected.columns.str.removesuffix("_GRCh38_1_modFraction") 
    hap2_data_frame_corrected.columns = hap2_data_frame_corrected.columns.str.removesuffix("_GRCh38_1_modFraction")
    # change chromStart and chromEnd to start and end if necessary
    if ('chromStart' in hap1_data_frame_corrected.columns) and ('chromEnd' in hap1_data_frame_corrected.columns):
        hap1_data_frame_corrected.rename(columns={'chromStart': 'start', 'chromEnd' : 'end'}, inplace=True)
    # repeat process for hap2_data_frame_corrected
    if ('chromStart' in hap2_data_frame_corrected.columns) and ('chromEnd' in hap2_data_frame_corrected.columns):
        hap2_data_frame_corrected.rename(columns={'chromStart': 'start', 'chromEnd' : 'end'}, inplace=True)
    # fix #chrom column to chrom
    if ('#chrom' in hap1_data_frame_corrected.columns):
        hap1_data_frame_corrected.rename(columns={'#chrom': 'chrom'}, inplace=True)
    # repeat process for hap2_data_frame_corrected
    if ('#chrom' in hap2_data_frame_corrected.columns):
        hap2_data_frame_corrected.rename(columns={'#chrom': 'chrom'}, inplace=True)
    # labeling, column drop, and joining depending on region
    if (region_type=="CGI"):
        # label sample columns as H1 for hap1 data frame - must remake whole column names list
        hap1_data_frame_corrected.columns=list(hap1_data_frame_corrected.columns[0:10])+list(hap1_data_frame_corrected.columns[10:] + "_H1")
        # label sample columns as H2 for hap2 data frame - must remake whole column names list
        hap2_data_frame_corrected.columns=list(hap2_data_frame_corrected.columns[0:10])+list(hap2_data_frame_corrected.columns[10:] + "_H2")
        # use FULL outer join of each hap data frame to fill in missing methylation regions as NA in respective hap
        # join on region name corrected above in rename_methylation_regions (data_frame['name'])
        joined_haps_data_frame=pd.merge(hap1_data_frame_corrected,hap2_data_frame_corrected,on=['chrom','start','end','name','length','cpgNum','gcNum','perCpg','perGc','obsExp'],how='outer')
    elif (region_type=="GB"):
        # label sample columns as H1 for hap1 data frame - must remake whole column names list
        hap1_data_frame_corrected.columns=list(hap1_data_frame_corrected.columns[0:6])+list(hap1_data_frame_corrected.columns[6:] + "_H1")
        # drop strand column - inconsistent between haplotypes in NABEC test...
        hap1_data_frame_corrected=hap1_data_frame_corrected.drop('strand',axis=1)
        # label sample columns as H2 for hap2 data frame - must remake whole column names list
        hap2_data_frame_corrected.columns=list(hap2_data_frame_corrected.columns[0:6])+list(hap2_data_frame_corrected.columns[6:] + "_H2")
        # drop strand column - inconsistent between haplotype in NABEC test...
        hap2_data_frame_corrected=hap2_data_frame_corrected.drop('strand',axis=1)
        # use FULL outer join of each hap data frame to fill in missing methylation regions as NA in respective hap
        # join on region name corrected above in rename_methylation_regions (data_frame['name'])
        joined_haps_data_frame=pd.merge(hap1_data_frame_corrected,hap2_data_frame_corrected,on=['chrom','start','end','name','dot'],how='outer')
    elif (region_type=="PROM"):
        # label sample columns as H1 for hap1 data frame - must remake whole column names list
        hap1_data_frame_corrected.columns=list(hap1_data_frame_corrected.columns[0:6])+list(hap1_data_frame_corrected.columns[6:] + "_H1")
        # drop strand column - inconsistent between haplotypes in HBCC test...
        hap1_data_frame_corrected=hap1_data_frame_corrected.drop('strand',axis=1)
        # label sample columns as H2 for hap2 data frame - must remake whole column names list
        hap2_data_frame_corrected.columns=list(hap2_data_frame_corrected.columns[0:6])+list(hap2_data_frame_corrected.columns[6:] + "_H2")
        # drop strand column - inconsistent between haplotype in HBCC test...
        hap2_data_frame_corrected=hap2_data_frame_corrected.drop('strand',axis=1)
        # use FULL outer join of each hap data frame to fill in missing methylation regions as NA in respective hap
        # join on region name corrected above in rename_methylation_regions (data_frame['name'])
        joined_haps_data_frame=pd.merge(hap1_data_frame_corrected,hap2_data_frame_corrected,on=['chrom','start','end','promoter_name','number'],how='outer')
    # return joined haps data frame
    return(joined_haps_data_frame)
    
# routine to prepare methylation map data frame from input methylation BED file(s)
def make_methylation_map(region_type,joined_haps_data_frame):
    # extract map from joined haps data frame
    methylation_map_df=pd.DataFrame(columns=['NAME','CHROM','START','STOP'])
    if (region_type=="CGI" or region_type=="GB"):
        methylation_map_df['NAME']=joined_haps_data_frame['name']
    elif (region_type=="PROM"):
        methylation_map_df['NAME']=joined_haps_data_frame['promoter_name']
    methylation_map_df['CHROM']=joined_haps_data_frame['chrom']
    methylation_map_df['START']=joined_haps_data_frame['start']
    methylation_map_df['STOP']=joined_haps_data_frame['end']
    return(methylation_map_df)

# routine to prepare methylation data matrix data frame from input methylation BED file(s)
def make_methylation_data(region_type,joined_haps_data_frame):
    # initialize data frame
    transposed_methylation_data_df=pd.DataFrame()
    # get column names for final matrix (region names from haplotype joined, sample relabeled, region renamed data frame)
    if (region_type=="CGI" or region_type=="GB"):
        methylation_data_col_names=joined_haps_data_frame['name'].tolist()
    elif (region_type=="PROM"):
        methylation_data_col_names=joined_haps_data_frame['promoter_name'].tolist()
    # add per sample per haplotype data to initial data frame
    if (region_type=="CGI"):
        transposed_methylation_data_df=joined_haps_data_frame.iloc[:,10:]
    elif (region_type=="GB"):
        transposed_methylation_data_df=joined_haps_data_frame.iloc[:,6:]
    elif (region_type=="PROM"):
        transposed_methylation_data_df=joined_haps_data_frame.iloc[:,6:]
    # transpose initial data frame
    methylation_data_df=transposed_methylation_data_df.transpose()
    # set column names of transposed data frame
    methylation_data_df.columns=methylation_data_col_names
    # turn index to column
    methylation_data_df['INDEX']=methylation_data_df.index
    # add haplotype column
    # add sample column
    methylation_data_df[['SAMPLE','HAPLOTYPE']] = methylation_data_df['INDEX'].str.rsplit("_",n=1,expand=True)
    # drop index column
    methylation_data_df=methylation_data_df.drop('INDEX',axis=1)
    # move sample and haplotype columns to first two columns
    methylation_data_df.insert(0, 'HAPLOTYPE', methylation_data_df.pop('HAPLOTYPE'))
    methylation_data_df.insert(0, 'SAMPLE', methylation_data_df.pop('SAMPLE'))
    # sort by BOTH sample and haplotype values (e.g., HBCC_81925_FTX,H1 then HBCC_81925_FTX,H2, then next sample)
    methylation_data_df=methylation_data_df.sort_values(by=["SAMPLE","HAPLOTYPE"])
    # output matrix
    return(methylation_data_df)
    
parser = argparse.ArgumentParser(description="Convert methylation BED files to methylation map and data matrix CSVs for downstream analysis (e.g., QTLs).")

# Define region type
parser.add_argument(
    "-t", "--region_type",
    required=True,
    choices=["CGI", "GB", "PROM"],
    help="Specifying region type for parsing and joining. Currently supported are CpG islands (CGI), gene bodies (GB), and promoters (PROM)."
)

# Add argument for haplotype 1 BED file
parser.add_argument(
    "-h1", "--haplotype_1",
    required=True,
    help="Methylation BED file for the first haplotype (H1)."
)

# Add argument for haplotype 2 BED file
parser.add_argument(
    "-h2", "--haplotype_2",
    required=True,
    help="Methylation BED file for the second haplotype (H2)."
)

# Add argument for region prefix
parser.add_argument(
    "-p", "--region_prefix",
    required=True,
    help="Prefix for the output methylation region names (e.g., CGIs, gene bodies, promoters)."
)

# Add argument for final output csv file
parser.add_argument(
    "-o", "--output_prefix",
    required=True,
    help="Prefix for the output methylation data and genetic map files."
)

# Add argument for missingness filter
parser.add_argument(
    "-m", "--missing_info_rate",
    required=False,
    default=0.05,
    type=float,
    help="Filter out regions lacking methylation information for higher than this proportion of samples (default 0.05 or 5%%)."
)

# Parse the arguments
args = parser.parse_args() 
# import haplotype 1 methylation bed file
# include period as NA value
hap1_meth_bed=pd.read_csv(args.haplotype_1,sep="\t",na_values=['.'])
# import haplotype 2 methylation bed file
# include period as NA value
hap2_meth_bed=pd.read_csv(args.haplotype_2,sep="\t",na_values=['.'])
# filter out regions with missing information rate over parameter
# default rate is 5% missing or less
hap1_meth_bed=filter_regions_by_missing_info_rate(args.region_type,hap1_meth_bed,args.missing_info_rate)
hap2_meth_bed=filter_regions_by_missing_info_rate(args.region_type,hap2_meth_bed,args.missing_info_rate)
# rename methylation regions
if (args.region_type=="CGI" or args.region_type=="GB"):
    # rename haplotype 1 methylation regions
    hap1_meth_bed_renamed=rename_methylation_regions(hap1_meth_bed,args.region_prefix)
    # rename haplotype 2 methylation regions
    hap2_meth_bed_renamed=rename_methylation_regions(hap2_meth_bed,args.region_prefix)
elif (args.region_type=="PROM"):
    # no need to rename either, so just assign new data frame with copy
    hap1_meth_bed_renamed=hap1_meth_bed.copy()
    hap2_meth_bed_renamed=hap2_meth_bed.copy()
# outer join haplotype files
joined_relabeled_hap_df=join_relabeled_haplotype_beds(args.region_type,hap1_meth_bed_renamed,hap2_meth_bed_renamed)
# make methylation map
output_methylation_map_df=make_methylation_map(args.region_type,joined_relabeled_hap_df)
# make methylation data matrix
output_methylation_data_df=make_methylation_data(args.region_type,joined_relabeled_hap_df)
# export methylation map as csv
output_methylation_map_df.to_csv(args.output_prefix+"_methylation_map.csv",index=False,header=True,na_rep="NA",chunksize=1000000)
# export methylation data matrix as csv
output_methylation_data_df.to_csv(args.output_prefix+"_methylation_data.csv",index=False,header=True,na_rep="NA",chunksize=1000)

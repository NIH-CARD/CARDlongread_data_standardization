# script to filter QTL output by p-value using BH FDR correction

import pandas as pd
import numpy as np
import time
import psutil
from statsmodels.stats import multitest
import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Perform multiple hypothesis testing false-discovery rate (FDR) correction for allele-specific methylation QTL output from all chromosomes.")
    parser.add_argument("--input", required=True, help="Path to input QTL file for all chromosomes; 'gene,window_size,outcome,predictor,beta,std_err,r2,p_value,N,predictor_freq' as header.")
    parser.add_argument("--output", required=True, help="Path to the output QTL file with additional two fields - rejected and pvalue-corrected.")
    parser.add_argument("--rejected", action=argparse.BooleanOptionalAction, default=False, required=False, help="Set this option to only print results where the null hypothesis was rejected (i.e., statistically significant associations).")
    return parser.parse_args()
    
def main():
    # Parse input and output arguments.
    args = parse_args()
    
    # Start timing and memory tracking
    start_time = time.time()
    process = psutil.Process()
    
    # Load input file
    input_QTL_df = pd.read_csv(args.input)
    
    # perform Benjamini-Hochberg multiple hypothesis testing correction
    # BH correction is default. Also unsorted p-values is default
    # fdr correction output is two arrays - first rejection (True or False for null hypothesis) and second corrected p-value
    input_QTL_df['rejected']=multitest.fdrcorrection(input_QTL_df.p_value)[0]
    input_QTL_df['pvalue-corrected']=multitest.fdrcorrection(input_QTL_df.p_value)[1]
    
    if args.rejected is True:
        # only print results where null hypothesis is rejected
        input_QTL_df=input_QTL_df[input_QTL_df['rejected']==True]
    
    # save output file; no indexes
    input_QTL_df.to_csv(args.output,index=False)
    
    # Print runtime and max RAM usage
    end_time = time.time()
    print(f"Execution Time: {end_time - start_time:.2f} seconds")
    print(f"Max RAM Usage: {process.memory_info().rss / (1024 ** 2):.2f} MB")
    
if __name__ == "__main__":
    main()
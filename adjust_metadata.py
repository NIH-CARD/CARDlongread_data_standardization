#!/usr/bin/python 

# adjust_metadata.py
# script to adjust Rylee/Melissa's metadata to meet that needed for Mike's script

import argparse
import pysam
from pysam import VariantFile
import pandas as pd
import numpy as np
import os
import csv

parser = argparse.ArgumentParser(description="Adjust Rylee/Melissa's metadata CSVs to meet Mike's requirements.")

# Add argument for input file
parser.add_argument(
    "-i", "--input_file",
    required=True,
    help="Path to input CSV metadata file."
)

# Add argument for final output csv file
parser.add_argument(
    "-o", "--output_file",
    required=True,
    help="Path to corrected CSV metadata output file."
)

# Parse the arguments
args = parser.parse_args() 
# import metadata file
input_metadata=pd.read_csv(args.input_file)
# adjust metadata
# first make SAMPLE and HAPLOTYPE columns
input_metadata[['SAMPLE','HAPLOTYPE']]=input_metadata['sample'].str.rsplit("_",n=1,expand=True)
# drop sample and HAPLOTYPE columns
input_metadata = input_metadata.drop(['sample','HAPLOTYPE'],axis=1)
# drop duplicate rows
input_metadata = input_metadata.drop_duplicates()
# rearrange SAMPLE column
input_metadata.insert(0, 'SAMPLE', input_metadata.pop('SAMPLE'))

# export fixed csv
input_metadata.to_csv(args.output_file,index=False,header=True)
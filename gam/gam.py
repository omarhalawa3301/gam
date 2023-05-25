# Importing modules
import gzip
import vcf
import pandas as pd
import argparse as ap

# Importing file processing functions
from gam_utils import *

#   __  __      __      __      _____         __      __  ______   _____     _____   _____    ____    _   _ 
#  |  \/  |     \ \    / /     |  __ \        \ \    / / |  ____| |  __ \   / ____| |_   _|  / __ \  | \ | |
#  | \  / |      \ \  / /      | |__) |        \ \  / /  | |__    | |__) | | (___     | |   | |  | | |  \| |
#  | |\/| |       \ \/ /       |  ___/          \ \/ /   |  __|   |  _  /   \___ \    | |   | |  | | | . ` |
#  | |  | |  _     \  /     _  | |       _       \  /    | |____  | | \ \   ____) |  _| |_  | |__| | | |\  |
#  |_|  |_| (_)     \/     (_) |_|      (_)       \/     |______| |_|  \_\ |_____/  |_____|  \____/  |_| \_|

#   _   _    ____      __  __       _             __     __  ______   _______ 
#  | \ | |  / __ \    |  \/  |     | |            \ \   / / |  ____| |__   __|
#  |  \| | | |  | |   | \  / |     | |             \ \_/ /  | |__       | |   
#  | . ` | | |  | |   | |\/| |     | |              \   /   |  __|      | |   
#  | |\  | | |__| |   | |  | |  _  | |____   _       | |    | |____     | |   
#  |_| \_|  \____/    |_|  |_| (_) |______| (_)      |_|    |______|    |_|   . . . (sadge)


"""
    Name:          Pei Ting Chua Chai, Omar Halawa
    Email:         pchai@ucsd.edu, ohalawa@ucsd.edu
    File name:     gam.py
    Project:       GAM (GWAS with Advanced Machine-learning) - CSE 185 Project
    Repository:    https://github.com/omarhalawa3301/gam
    Description:   GAM main python script to: 
                    - Process genotype (vcf.gz) & phenotype (.phen) file inputs
                    - Process optional algorithm specification arguments
                    - Call functions from gam_utils.py to process files as DFs
                    - Perform GWAS (linear-regression based) by:
                      1. Creating linear-regression models of each SNP
                      2. 
                    - Predict feature dataset and compare to "true" target file
                   Outputs an .assoc.linear file.
                   Also outputs a Manhattan plot visualization.
                   Designed to allow for further file type implementation.
                   Created for CLI implementation.
                   
    References:    https://www.cog-genomics.org/plink/1.9/assoc#linear
"""

# Adding arguments to script for genotype & phenotype file inputs,
# .assoc.linear file output path, algorithm specifics parameters, & debugging
parser = ap.ArgumentParser(description='GAM')

# Adding file input arguments (required)
# Genotype file input (.vcf.gz):
# TODO, allow for ".vcf" file input too?
parser.add_argument("-g", "--genotype", help="genotype data filename"
                    + " Valid file format(s): .vcf.gz", required=True)
# Phenotype file input (.phen):
parser.add_argument("-p", "--phenotype", help="phenotype data filename"
                    + " Valid file format(s): .phen", required=True)

# Optional arguments for result file's output directory and  
# Assigning result file's name (.assoc.linear) (optional, has default value):
parser.add_argument("-o", "--out", help="path to output results file",
                    nargs="?", const=1)

# Option of using simple linear regression 
parser.add_argument("-l", "--linear", help="do GWAS using linear regression", nargs="?", required=False)

# Option of using an ensemble of linear regression models
parser.add_argument("-le", "--linear-ensemble", help="do GWAS using an ensemble of linear models", nargs="?", required=False)

# Option of using boosted decision trees (sklearn.ensemble.Adaboost)
parser.add_argument("-bdt", "--boosted", help="do GWAS using boosted decision trees", nargs="?", required=False)


# Parsing arguments for future calls within script to utilize
args = parser.parse_args()


geno_file = file_valid(args.genotype, "GENOTYPE")
pheno_file = file_valid(args.phenotype, "PHENOTYPE")

with gzip.open(input_file, 'rb') as f_in





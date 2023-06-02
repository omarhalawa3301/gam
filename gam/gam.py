# Importing modules
import pandas as pd
import argparse as ap
import time
from scipy import stats

# Importing file processing functions
# from .gam_utils import *
from gam_utils import *

# Importing genotype and phenotype string markers (more to implement for later)
# from .Marker import *
from Marker import *

# Importing version
# from gam import __version__

start = time.time()


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
                      2. Calculating p-value score for each SNP
                      3. Analyze and visualize resulting data  
                   Outputs an .assoc.linear file.
                   Also outputs a Manhattan plot visualization.
                   Designed to allow for further file type implementation.
                   Created for CLI implementation.
                   
    References:    https://www.cog-genomics.org/plink/1.9/assoc#linear
                   https://tiny.cc/gam-float
"""

# Custom function for limiting argument's values between 0 and 1, inclusive of both
# See second reference for code source
def restricted_float(x):
    try:
        x = float(x)
    except ValueError:
        raise ap.ArgumentTypeError("%r not a floating-point literal" % (x,))

    if x <= 0.0 or x >= 1.0:
        raise ap.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x


def main():

    # Adding arguments to script for genotype & phenotype file inputs,
    # .assoc.linear file output path, algorithm specifics parameters, & debugging
    parser = ap.ArgumentParser(description='GAM')

    # Adding file input arguments (required)
    # Genotype file input (.vcf.gz):
    # TODO, allow for ".vcf" file input too?
    parser.add_argument("-g", "--genotype", help="genotype data filename "
                    + "Valid file format(s): .vcf.gz", required=True, type=str)
    # Phenotype file input (.phen):
    parser.add_argument("-p", "--phenotype", help="phenotype data filename "
                    + "Valid file format(s): .phen", required=True, type=str)


    # Optional arguments for result file's output directory and  
    # Assigning result file's name (.assoc.linear) (optional, has default value):
    parser.add_argument("-o", "--out", help="path to output results file", 
                        nargs="?", default=".", type=str, required=False)

    # Assigning parameters to the three possible modes of running GWAS
    modes = parser.add_mutually_exclusive_group(required=True)

    # Option of using simple linear regression 
    modes.add_argument("-l", "--linear", help="do GWAS using linear regression", action='store_true', required=False)

    # Option of using an ensemble of linear regression models
    modes.add_argument("--le", "--linear-ensemble", help="do GWAS using an ensemble of linear models", action='store_true', required=False)

    # Option of using boosted decision trees (sklearn.ensemble.Adaboost)
    modes.add_argument("--bdt", "--boosted", help="do GWAS using boosted decision trees", action='store_true', required=False)


    # Option of specifying MAF (minor allele frequency) by which to filter out SNPs that have a minor allele frequency less than specified
    parser.add_argument("-m", "--maf", help="assign minor allele frequency for filtration", nargs="?", type=restricted_float, required=False)

    # TODO: Covariate implementation

    # Parsing arguments for future calls within script to utilize
    args = parser.parse_args()

    pheno_ext = file_valid(args.phenotype, Marker.PT)
    geno_ext = file_valid(args.genotype, Marker.GT)

    # Obtain pandas data frames from process function (number of samples only matters for genotype file processing)
    # Also obtaining number of samples and snps for further data processing in pipeline 
    pheno_df = process(args.phenotype, pheno_ext, Marker.NA)
    num_samples = len(pheno_df.index)
    geno_df = process(args.genotype, geno_ext, num_samples)
    num_snps = len(geno_df.index)

    # Initializing basename for outputting visuals and 
    basename = args.genotype.replace(geno_ext,"")

    # GWAS mode check:
    # Performing "normal" linear regression
    if (args.linear):

        # Genotype df looks as follows:
        #              SAMPLE1 SAMPLE2 SAMPLE3 . . .
        #        SNP1    0        0       1
        #        SNP2    2        2       1
        #        SNP3    1        1       0
        #          .
        #          .
        #          .  

        # Phenotype df looks as follows:
        #   
        #        SAMPLE1   1.642
        #        SAMPLE2   -0.132
        #           .
        #           .
        #           .

        # values to retain from linear regression for plotting
        chrs, snp_ids, bps, alts, test_labels, nmiss, beta, t_stat, pvals = ([] for i in range(Marker.ASSOC_COLS))

        # Iterating through each biallelic SNP that we have taken in
        for snp in range(0, num_snps):

            # Genotype values for a specific SNP among samples (0, 1, or 2)
            X = geno_df.iloc[snp].values

            # Phenotype values for a specific SNP among samples 
            y = pheno_df["Phenotype"].values

            # SCIPY linear regression model for SNP
            slope, intercept, r_value, p_value, std_err = stats.linregress(X, y)

            # Filtering data needed for assoc.linear output file
            chrs.append(str(geno_df.index[snp][0]))
            bps.append(str(geno_df.index[snp][1]))
            snp_ids.append(str(geno_df.index[snp][2]))
            alts.append(str(geno_df.index[snp][3]))
            test_labels.append("ADD")
            nmiss.append(str(num_samples))
            beta.append(str(slope))
            t_stat.append(str(slope/std_err))
            pvals.append(str(p_value))
            
        # Bonferroni Correction: adjust the given p-value threshold by number of tests (SNPs)
        sig_thresh = Marker.P_STD / len(geno_df.index)
        print("GWAS Bonferroni-Corrected P-Value: ", sig_thresh)

        # Obtaining output assoc.linear filename and creating it
        assoc_filename = basename+".assoc.linear"
        assoc_result(chrs, snp_ids, bps, alts, test_labels, nmiss, beta, t_stat, pvals, assoc_filename)
    
    # # Else , performing
    # elif (args.linear_ensemble):
    #     # TODO: Write the process for linear ensemble GWAS
    #     return
    # elif (args.boosted):
    #     # TODO: Write the process for boosted decision trees GWAS
    #     return

    # Plotting visualizations (Manhattan and qqplot)
    # Uses genotype file's basename for plot graphs
    plot(assoc_filename, basename)
        

if __name__ == "__main__":
    main()

end = time.time()
print("Total runtime:", end - start, "seconds")

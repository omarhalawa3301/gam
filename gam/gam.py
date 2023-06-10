# Importing modules
import pandas as pd
import argparse as ap
import numpy as np
import time
import os
from art import *
from sklearn.linear_model import LinearRegression
from scipy import stats
import statsmodels.api as sm

# Importing file processing functions
from .gam_utils import *
# Importing genotype and phenotype string markers (more to implement for later)
from .Marker import *
# Importing version
from gam import __version__

"""
    Name:          Pei Ting Chua Chai, Omar Halawa
    Email:         pchai@ucsd.edu, ohalawa@ucsd.edu
    File name:     LinearRegression.py
    Project:       GAM (GWAS with Advanced Machine-learning) - CSE 185 Project
    Repository:    https://github.com/omarhalawa3301/gam
    Description:   Class file that extends LinearRegression from Scikit-learn to
                   calculate t-statistics and p-values for model coefficients (betas) 
                   
    References:    Stolen (based) from https://gist.github.com/brentp/5355925 
"""


start = time.time()

tprint("GAM", "isometric1")
catchphrase = text2art("If you didn't want who I GWAS, you don't deserve who I GAM...", font="subscript2") + " " + art("do you even lift bro") + "!!!"
print(catchphrase, "\n")

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
                   https://stackoverflow.com/a/42677750/21989720
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

    # Optional arguments for result files' output directory and names, MAF filtration, and confounding factor control via covariates
    # Assigning result file's name (.assoc.linear) (optional, has default value):
    parser.add_argument("-d", "--dir", help="directory to output result files", 
                        nargs="?", default=".", type=str, required=False)
    parser.add_argument("-o", "--out", help="basename of result files", 
                        nargs="?", type=str, required=False)
    parser.add_argument("--maf", help="minor allele frequency, filter out SNPs with MAF less than assigned threshold", 
                        nargs="?", default=0.0, type=float, required=False)
    parser.add_argument("--covar", help="takes in .eigenvec file of PCA run to account for confounding factors via covariates", 
                        nargs="?", default=None, type=str, required=False)

    # Assigning parameters to the three possible modes of running GWAS
    modes = parser.add_mutually_exclusive_group(required=True)

    # Option of using simple linear regression 
    modes.add_argument("-l", "--linear", help="do GWAS using linear regression", action='store_true', required=False)

    # Option of using an ensemble of linear regression models
    modes.add_argument("--le", "--linear-ensemble", help="do GWAS using an ensemble of linear models", action='store_true', required=False)

    # Option of using boosted decision trees (sklearn.ensemble.Adaboost)
    modes.add_argument("--bdt", "--boosted", help="do GWAS using boosted decision trees", action='store_true', required=False)


    # # Option of specifying MAF (minor allele frequency) by which to filter out SNPs that have a minor allele frequency less than specified
    # parser.add_argument("-m", "--maf", help="assign minor allele frequency for filtration", nargs="?", type=restricted_float, required=False)    

    # Parsing arguments for future calls within script to utilize
    args = parser.parse_args()

    pheno_ext = file_valid(args.phenotype, Marker.PT)
    geno_ext = file_valid(args.genotype, Marker.GT)

    # Exiting after printing error message(s) if invalid inputs
    if (pheno_ext == None or geno_ext == None):
        exit()

    # Obtain pandas data frames from process function (number of samples only matters for genotype file processing)
    # Also obtaining number of samples and snps for further data processing in pipeline 
    pheno_df = process(args.phenotype, pheno_ext, Marker.NA, args.maf)
    num_samples = len(pheno_df.index)
    geno_df = process(args.genotype, geno_ext, num_samples, args.maf)
    num_snps = len(geno_df.index)

    if (args.covar != None):
        covar_df = covar_process(args.covar)

    # Initializing basename for outputting visuals and assoc.linear result file
    if (args.out != None):
        basename = args.out
    else:
        basename = args.genotype.replace(geno_ext,"")

    # Input formatting check as to not cause issues with --dir parameter
    if (basename.startswith("./")):
        basename = basename.replace("./","")

    # Obtaining output assoc.linear filename and creating it
    if (args.dir != "." and args.dir.endswith("/")):
        assoc_filename = args.dir + basename.split("/")[-1] + ".assoc.linear"
        plot_filename = args.dir + basename.split("/")[-1] + ".png"
    elif (args.dir != "."):
        assoc_filename = args.dir + "/" +basename.split("/")[-1] + ".assoc.linear"
        plot_filename = args.dir + "/" +basename.split("/")[-1] + ".png"
    else:    
        assoc_filename = basename + ".assoc.linear"
        plot_filename = basename + ".png"

    # Making output directory if it does not exist
    if (not exists(args.dir)):
        os.makedirs(args.dir)

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
            y = pheno_df["Phenotype"]

            # No covariate processing
            if (args.covar == None):
                # SCIPY linear regression model for SNP
                slope, intercept, r_value, p_value, std_err = stats.linregress(X, y.values)
                t = slope/std_err

            # Covariate processing
            else:

                # Creating dataframe of all variables (genotype column and PC columns)
                new_df = pd.concat([geno_df.iloc[snp], covar_df], axis=1)
                
                # Linear regression on the new dataframe
                model = LinearRegression()
                model.fit(new_df,y)

                # Get the model's intercept and coefficients in one numpy array
                params = np.append(model.intercept_,model.coef_)

                # Predict on the data we fitted the model to
                predictions = model.predict(new_df)

                # Code from last link in references
                newX = pd.DataFrame({"Constant":np.ones(len(new_df))}).join(pd.DataFrame(new_df.reset_index(drop=True)))
                MSE = (sum((y-predictions)**2))/(len(newX)-len(newX.columns))

                var_b = MSE*(np.linalg.inv(np.dot(newX.T,newX)).diagonal())
                sd_b = np.sqrt(var_b)
                ts_b = params/ sd_b
                
                p_values =[2*(1-stats.t.cdf(np.abs(i),(len(newX)-len(newX.columns)-1))) for i in ts_b]

                myDF3 = pd.DataFrame()
                myDF3["Coefficients"],myDF3["Standard Errors"],myDF3["t values"],myDF3["Probabilities"] = [params,sd_b,ts_b,p_values]
                p_value = (myDF3["Probabilities"][1])
                t = (myDF3["t values"][1])
                slope = model.coef_[0]

            # Filtering data needed for assoc.linear output file
            chrs.append(str(geno_df.index[snp][0]))
            bps.append(str(geno_df.index[snp][1]))
            snp_ids.append(str(geno_df.index[snp][2]))
            alts.append(str(geno_df.index[snp][3]))
            test_labels.append("ADD")
            nmiss.append(str(num_samples))
            beta.append(str(slope))
            t_stat.append(str(t))
            pvals.append(str(p_value))
            
        # Bonferroni Correction: adjust the given p-value threshold by number of tests (SNPs)
        sig_thresh = Marker.P_STD / len(geno_df.index)

        assoc_result(chrs, snp_ids, bps, alts, test_labels, nmiss, beta, t_stat, pvals, assoc_filename)
    
    # # Else , performing linear ensemble GWAS
    elif (args.linear_ensemble):

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

        assoc_result(chrs, snp_ids, bps, alts, test_labels, nmiss, beta, t_stat, pvals, assoc_filename)

    # elif (args.boosted):
    #     # TODO: Write the process for boosted decision trees GWAS
    #     return

    # Plotting visualizations (Manhattan and qqplot)
    # Uses genotype file's basename for plot graphs
    plot(assoc_filename, plot_filename)
    
    end = time.time()
    print("\n\nTotal runtime:", end - start, "seconds")

if __name__ == "__main__":
    main()

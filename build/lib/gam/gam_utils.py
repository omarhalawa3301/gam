# Importing modules
from os.path import exists
import pandas as pd
import gzip
from cyvcf2 import VCF
import matplotlib.pyplot as plt
from qqman import qqman

# Importing constant marker and extension strings
from .Marker import *
from .Extension import *
# from Marker import *
# from Extension import *


"""
    Name:          Pei Ting Chua Chai, Omar Halawa
    Email:         pchai@ucsd.edu, ohalawa@ucsd.edu
    File name:     gam_utils.py
    Project:       GAM (GWAS with Advanced Machine-learning) - CSE 185 Project
    Repository:    https://github.com/omarhalawa3301/gam
    Description:   GAM function-containing python script for main gam.py to use

    Functions:     file_valid: checks genotype & phenotype input file validity
                   
    References:    https://www.cog-genomics.org/plink/1.9/assoc#linear
"""

def file_valid(name, marker):
    """ Function that checks for filename validity through existence & format.
        Has general format to allow for future file type implementation.       

    Arguments:
        name:   filename (with extension) to check existence from current dir
        marker: String value indicating whether file is of geno or pheno data
                (see Marker class)
    Returns:
        ext:  a String of the file extension if it passes name and ext check
        None: if the file name or extension are not valid
    """

    # Checking for case of no provided input (test dataset input)
    if (name == None):
        return None

    # Intializing placeholder for proper file extension list
    curr_list = None
    # Initalizing placeholder for file extension
    curr_ext = None

    # Validating existent file name
    valid_name = exists(name)

    # Intializing extension check as False by default
    valid_ext = False

    # Validating file extension through checking marker value (by Marker class)
    # Then assigning curr_list's value to one of the two Extension class lists
    # Only two cases possible, calls are within program
    
    # Feature data file case
    if (marker == Marker.GT):
        curr_list = Extension.GT_EXT
    # Target data file case
    elif (marker == Marker.PT):
        curr_list = Extension.PT_EXT

    # Carrying out the extension check logic using the now-updated curr_list
    for ext in curr_list:
        if name.endswith(ext):
            # If match is found, updating valid_ext to True
            valid_ext = True
            # Also assigning curr_ext's value
            curr_ext = ext
            # Breaking once match is found
            break

    # Accounting for both checks 
    if (valid_name and valid_ext):
        # if checks are valid, returning file extension
        return curr_ext
    else:
        # Invalid file name message
        if (not valid_name):
            print("File name '" + name + 
                "' is invalid (not found from current directory).")
        # Invalid file extension message
        if (not valid_ext):
            # Note: Future file type implementation must list all valid
            # extensions for classifier input or target, respectively 
            print("File extension of '" + name + "' is invalid. " + 
                "Expected extension(s): ", end="")
            print(curr_list)
        
        print()
        # Returning None in the case of invalid file name
        return None


def process(name, ext, num_samples):
    """ Function that processes valid file given its extension as an argument

    Arguments:
        name:   name of file to process
        ext:    extension of file to process (obtained from file_valid call)
        num_samples:  number of samples
    Returns:    returns a call to the appropriate helper function that contains
                actual logic for processing
    """

    if (ext == Extension.VCF_GZ_EXT):
        return gene_process(name, ext, num_samples)
    elif (ext == Extension.PHENE_EXT):
        return phene_process(name)


def filter_alts(ref, arr):
    """ Function that filters the list of all alleles seen among samples for an SNP 
        by outputting the most frequent non-reference alternative allele

    Arguments:
        ref:    reference allele as a string
        arr:    array of genotypes in the format ['G|C', 'C|G', 'G|G']
    Returns:    string of most frequent non-reference alternative allele, -1 if only ref was seen
    """

    # counts: {"G": 3, "A": 43, "TAC": 421}
    counts = {}

    # fill in counts to find most frequent alternate
    for alt in arr:
        first,second = alt.split("|", 1)
        # gt = set([first, second])
        
        if (first not in counts):
            counts[first] = 1
        else:
            counts[first] += 1

        if (second not in counts):
            counts[second] = 1
        else:
            counts[second] += 1

    curr_max = max(counts, key=counts.get)

    if (curr_max == ref and len(counts.items()) == 1):
        # we ignore this SNP as it only has one allele for every sample
        return -1
    elif (curr_max == ref):
        # return the most common allele that isn't the reference
        counts.pop(ref)
        return max(counts, key=counts.get)
    elif (len(counts.items()) > 2):
        # there are more than 2 alleles including ref and most common alt
        return -1
    else:
        # most common allele is an alternative
        return curr_max


def extract_sample_names(filename, num_samples):
    """ Function that takes all the relevant sample names in the genotype file

    Arguments:
        filename:     name of the genotype file to extract sample names from  
        num_samples:  number of samples for index range to extract from array
    Returns:    list of relevant sample names in alphanumeric order
    """
    sample_names = []

    # Reading until the last line of the header is found (is first line without ##)
    for line in gzip.open(filename, "rt"):

        if not line.startswith("##") and line.startswith("#"):
            line = line.replace("\n", "")
            arr = line.split("\t")
            break

    # This line has the sample names
    sample_names = arr[len(arr)-num_samples:]

    return sample_names


def gene_process(name, ext, num_samples):
    """ Function that processes valid file given its name and extension as an argument

    Arguments:
        name:         name of file to process
        ext:          extension of file to process (obtained from file_valid call)
        num_samples:  number of samples
    Returns:    dataframe of rows=SNPs and cols=samples with value of genotype (0=homo ref, 1=hetero, 2=homo alt) 
    """
    
    # Output df looks as follows:
    #              SAMPLE1 SAMPLE2 SAMPLE3 . . .
    #        SNP1    0        0       1
    #        SNP2    2        2       1
    #        SNP3    1        1       0
    #          .
    #          .
    #          .                  

    # Genotype file processing for vcf.gz file format
    if (ext == Extension.VCF_GZ_EXT):

        sample_names = extract_sample_names(name, num_samples)

        arr = []
        snp_ids = []
        
        # Iterating through each SNP
        for variant in VCF(name): # or VCF('some.bcf')
            # Reference allele
            ref = variant.REF

            # Most common alternative allele
            alt = filter_alts(ref, variant.gt_bases)

            # Checking for case of no alternative allele present, only ref, so skip SNP
            if (alt == -1):
                continue
            
            # Add SNP ID consisting of (CHR#,POSITION ON CHR) because this SNP will be considered
            snp_ids.append((int(variant.CHROM), variant.POS, variant.ID, alt))
            
            # Possible genotypes for this valid SNP
            homo_ref = set([ref,ref])
            hetero = set([ref,alt])
            homo_alt = set([alt,alt])

            # List of genotypes for this specific SNP
            snp_gts = []

            # For each sample's genotype of the SNP
            for gt in variant.gt_bases:
                first,second = gt.split("|", 1)

                curr_gt = set([first,second])

                if (curr_gt == homo_ref):
                    snp_gts.append(0)
                if (curr_gt == hetero):
                    snp_gts.append(1)
                if (curr_gt == homo_alt):
                    snp_gts.append(2)

            arr.append(snp_gts)
            
        df = pd.DataFrame(arr, index=snp_ids)

        df.columns = sample_names

        # Sorting columns to match .phen file alphabetical and numeric sorting
        df = df.reindex(sorted(df.columns), axis=1)

    return df


def phene_process(name):
    """ Function that processes a phenotype file into a dataframe

    Arguments:
        name:   name of file to process
    Returns:    dataframe of sample vs corresponding phenotype value
    """
    # Assumes that the phenotype file has two columns of sample names to skip
    df = pd.read_csv(name, usecols=[1,2], sep = '\s+', header=None)

    df.columns = ["Sample","Phenotype"]

    # Sorting by sample name
    df = df.reindex(sorted(df.columns), axis=1)

    return df


def assoc_result(chrs, snp_ids, bps, alts, test_labels, nmiss, beta, t_stat, pvals, filename):
    """ Function that creates and outputs the result assoc.linear file to the specified output directory 
    Arguments:
        chrs:        array of chromosome numbers
        snp_ids:     array of SNP ids
        bps:         array of SNP positions on chromosomes by base pairs
        alts:        array of each SNP's alternative allele
        test_labels: array of correlated found in data, is additive by default
        nmiss:       array of number of samples used for each SNP's linear regression model
        beta:        array of slopes for each SNP's linear regression model
        t_stat:      array of t-statistic for each SNP's linear regression model (is slope/slope_std_err)
        pvals:       array of p-values for each SNP's linear regression model, p-value of line's slope being 0 
        filename:    name of file with the output directory included, user can only adjust the output directory, is . by default
    Returns:    Does not return anything
    """

    f = open(filename, "w")

    f.write("CHR\tSNP\tBP\tA1\tTEST\tNMISS\tBETA\tSTAT\tP\n")

    for i in range(0,len(snp_ids)):
        f.write(chrs[i] + "\t" + snp_ids[i]+ "\t" + bps[i]+ "\t" + alts[i]+ "\t" + test_labels[i]
                + "\t" + nmiss[i]+ "\t" + beta[i]+ "\t" + t_stat[i]+ "\t" + pvals[i]+ "\n")


def plot(in_name, out_name):
    """ Function that plots a Manhattan plot and a qqplot plot of the data 

    Arguments:
        in_name:    Name of assoc.linear result file to obtain data from
        out_name:   Name to use as prefix for manhattan and qqplot png files
    Returns:        Does not return anything
    """
    data = pd.read_csv(in_name, delim_whitespace=True)
    fig, (ax0, ax1) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [2, 1]})
    fig.set_size_inches((15, 5))
    qqman.manhattan(data, ax=ax0)
    qqman.qqplot(data, ax=ax1, out=out_name+"_qqplot.png")
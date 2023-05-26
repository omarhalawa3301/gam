# Importing modules
from os.path import exists
import pandas as pd
from cyvcf2 import VCF

# Importing constant marker and extension strings
from Marker import *
from Extension import *

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


def process(name, ext):
    """ Function that processes valid file given its extension as an argument

    Arguments:
        name:   name of file to process
        ext:    extension of file to process (obtained from file_valid call)
    Returns:    returns a call to the appropriate helper function that contains
                actual logic for processing
    """

    if (ext == Extension.VCF_GZ_EXT):
        return gene_process(name, ext)
    elif (ext == Extension.PHENE_EXT):
        return phene_process(name)


def gene_process(name, ext):
    """ Function that processes valid file given its name and extension as an argument

    Arguments:
        name:   name of file to process
        ext:    extension of file to process (obtained from file_valid call)
    Returns:    dataframe of rows=sample and cols=snps with value of genotype (0=homo ref, 1=hetero, 2=homo alt) 
    """
    
    # Genotype file processing for vcf.gz file format
    if (ext == Extension.VCF_GZ_EXT):
        for variant in VCF(name): # or VCF('some.bcf')

            print(variant.ALT) # worst case scenario, process if 0,1,or 2 via variant.gt_bases
            print("HELLOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO")

            # TODO: Figure out what the fuck we want to do when there are multiple alternate alleles
            # TODO: When this ^ is done, transform data into pandas df for further processing in pipeline

            # variant.REF, variant.ALT # e.g. REF='A', ALT=['C', 'T']

            # variant.CHROM, variant.start, variant.end, variant.ID, \
            #             variant.FILTER, variant.QUAL

            # # numpy arrays of specific things we pull from the sample fields.
            # # gt_types is array of 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
            # variant.gt_types, variant.gt_ref_depths, variant.gt_alt_depths # numpy arrays
            # variant.gt_phases, variant.gt_quals, variant.gt_bases # numpy array

    # return df


def phene_process(name):
    """ Function that processes a phenotype file into a dataframe

    Arguments:
        name:   name of file to process
    Returns:    dataframe of sample vs corresponding phenotype value
    """
    # Assumes that the phenotype file has two columns of sample names to skip
    df = pd.read_csv(name, usecols=[1,2], sep = '\s+', header=None)

    df.columns = ["Sample","Phenotype"]

    return df


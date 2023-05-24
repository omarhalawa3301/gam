# Importing modules
from os.path import exists

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
    if (marker == Marker.FEAT):
        curr_list = Extension.FEAT_EXT
    # Target data file case
    elif (marker == Marker.TAR):
        curr_list = Extension.TAR_EXT

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
                "' is invalid (not found in current directory).")
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


#!/usr/bin/env python3

"""
    Name:          Pei Ting Chua Chai, Omar Halawa
    Email:         pchai@ucsd.edu, ohalawa@ucsd.edu
    File name:     Extension.py
    Project:       GAM (GWAS with Advanced Machine-learning) - CSE 185 Project
    Repository:    https://github.com/omarhalawa3301/gam
    Description:   GAM constant extension string python-containing script                   
    References:    https://www.cog-genomics.org/plink/1.9/assoc#linear
"""


# Initializing Extension class and its variables for file_valid logic
# Designed to allow for smooth implementation of other file types
class Extension:
    ASSOC_LINEAR_EXT = ".assoc.linear"
    VCF_GZ_EXT = ".vcf.gz"
    PHENE_EXT = ".phen"
    GT_EXT =  [VCF_GZ_EXT]
    PT_EXT = [PHENE_EXT]
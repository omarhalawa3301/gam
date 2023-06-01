#!/usr/bin/env python3

"""
    Name:          Pei Ting Chua Chai, Omar Halawa
    Email:         pchai@ucsd.edu, ohalawa@ucsd.edu
    File name:     Marker.py
    Project:       GAM (GWAS with Advanced Machine-learning) - CSE 185 Project
    Repository:    https://github.com/omarhalawa3301/gam
    Description:   GAM marker string python-containing script.
    References:    https://www.cog-genomics.org/plink/1.9/assoc#linear
"""


# Initializing Marker class and its variables for file_valid argument markers and constants
class Marker:
    GT = "GENOTYPE"
    PT = "PHENOTYPE"
    NA = 0
    P_STD = 0.05
    ADDITIVE = "ADD"
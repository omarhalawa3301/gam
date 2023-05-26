[Plink repo](https://github.com/chrchang/plink-ng/tree/master)


vcf.gz -> vcf -> process into pandas df

take in phen  -> process into pandas df 

May also implement .bed input later?? 

And then, perform linear regression on each SNP via a set of for-loops, etc..

Manhattan plot at end (get from lab qqman)

observed vs expected -log10(p) plot

When parse header in:

#filter in header, if below a certain threshold, do filtering pipeline

Verify that header information matches pandas df information (# of samples, # of snps, etc..)


Include output directory

Allow no sex (some samples dont have sex in dataset) Does this matter (if there was important snp in y chromosome?) 

if sex is allowed then qqman plot would have chr 23 too?

format after the header:
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE1 SAMPLE2 SAMPLE3...
where CHROM and POS tell the position of the snp (whose id is ID), 
REF and ALT are reference and alternative alleles,
FILTER will say PASS if this snp meets all of the requirements and will have an id which tells us how it is disqualified from the analysis
FORMAT is the format of each sample for that snp (ex. allele|allele:genotype quality)

From the format of the vcf, we can take the genotypes of all the samples for every snp and put it into a dataframe/numpy array:  
sample1, 2, 3
0|0 0|1 0|1 snp  
0|0 0|1 0|1 snp  
0|0 0|1 0|1 snp  

And then from the phenotype file, we can add the value of the phenotype at the end of the columns? 

MAF??

Take in vcf file via:

	import vcf
	vcf_reader = vcf.Reader(filename="in.vcf.gz")
	https://pyvcf.readthedocs.io/en/latest/

OR

	import pandas as pd
	import gzip

	def get_vcf_names(vcf_path):
	    with gzip.open(vcf_path, "rt") as ifile:
	          for line in ifile:
	            if line.startswith("#CHROM"):
	                  vcf_names = [x for x in line.split('\t')]
	                  break
	    ifile.close()
	    return vcf_names


	names = get_vcf_names('file.vcf.gz')
	vcf = pd.read_csv('file.vcf.gz', compression='gzip', comment='#', chunksize=10000, delim_whitespace=True, header=None, names=names)


OR

	#!/usr/bin/env python

	import io
	import os
	import pandas as pd


	def read_vcf(path):
	    with open(path, 'r') as f:
	        lines = [l for l in f if not l.startswith('##')]
	    return pd.read_csv(
	        io.StringIO(''.join(lines)),
	        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
	               'QUAL': str, 'FILTER': str, 'INFO': str},
	        sep='\t'
	    ).rename(columns={'#CHROM': 'CHROM'})

**NOTE:** 'QUAL': str ? Should be int?

Also take in a phen file
	python open file to read in, very simple format

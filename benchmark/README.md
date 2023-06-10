# Benchmarking GAM -l against plink --linear (uses small test)

## Timing

```
time gam -g ./example-files/test.vcf.gz -p ./example-files/test.phen -l --dir test_run_results --out run1 --covar ./example-files/lab3_gwas.eigenvec

time plink --linear hide-covar --covar example-files/lab3_gwas.eigenvec --allow-no-sex -g example-files/test.vcf.gz -p example-files/test.phen --maf 0.05 --out plink_run
```

## Memory usage (uses small test)

See: https://github.com/jhclark/memusg
Note, this command was failing on OSX but worked on Ubuntu

```
memusg gam -g ./example-files/test.vcf.gz -p ./example-files/test.phen -l --dir test_run_results --out run1 --covar ./example-files/lab3_gwas.eigenvec

memusg plink --linear hide-covar --covar example-files/lab3_gwas.eigenvec --allow-no-sex -g example-files/test.vcf.gz -p example-files/test.phen --maf 0.05 --out plink_run
```

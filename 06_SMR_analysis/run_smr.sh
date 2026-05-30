#!/bin/bash
# Run SMR analysis for all 6 colocalised loci
# Author: Ruthushree P Basappa
# Date: May 2026
# Reference: Zhu et al. 2016, Nature Genetics

# First run prepare_gwas_smr.R to generate gwas_smr_chrbp.txt
# and g1000_eur_chrbp/ directory

./smr_x86 \
  --bfile g1000_eur_chrbp/g1000_eur \
  --gwas-summary gwas_smr_chrbp.txt \
  --beqtl-summary FB_Brain_2 \
  --out smr_all_genes \
  --extract-probe all_probes.txt \
  --peqtl-smr 5e-8 \
  --diff-freq 0.4 \
  --diff-freq-prop 0.6 \
  --thread-num 4

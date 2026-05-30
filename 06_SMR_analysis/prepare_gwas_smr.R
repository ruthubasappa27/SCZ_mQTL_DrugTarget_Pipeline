# Prepare GWAS summary statistics for SMR analysis
# Converts rsID format to chr:bp format to match FB_Brain_2 BESD file
# Author: Ruthushree P Basappa
# Date: May 2026

library(data.table)

# Load QC'd GWAS file
gwas_qc <- fread("gwas_qc.csv")

# Create chr:bp SNP ID to match FB_Brain_2 format
gwas_qc[, SNP_chrbp := paste0(CHR, ":", POS)]

# Create SMR format file
gwas_smr_chrbp <- gwas_qc[, .(
  SNP = SNP_chrbp,
  A1 = EA,
  A2 = OA,
  freq = MAF,
  b = BETA_gwas,
  se = SE_gwas,
  p = P_gwas,
  N = NEFF
)]

fwrite(gwas_smr_chrbp, "gwas_smr_chrbp.txt", sep="\t")
cat("GWAS SMR file created:", nrow(gwas_smr_chrbp), "SNPs\n")

# Also convert LD reference BIM file to chr:bp format
system("mkdir -p g1000_eur_chrbp")
system("cp 'g1000_eur 2/g1000_eur.fam' g1000_eur_chrbp/g1000_eur.fam")
system("cp 'g1000_eur 2/g1000_eur.bed' g1000_eur_chrbp/g1000_eur.bed")
system("awk '{print $1\"\\t\"$1\":\"$4\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6}' 'g1000_eur 2/g1000_eur.bim' > g1000_eur_chrbp/g1000_eur.bim")
cat("LD reference converted to chr:bp format\n")

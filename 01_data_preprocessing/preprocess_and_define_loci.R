# SCZ mQTL Drug Target Prioritisation Pipeline
# Script 1: Data Preprocessing and Locus Definition
# Author: Ruthushree P Basappa
# Institution: University of Galway
# Supervisor: Prof. Derek Morris
# Date: May 2026

library(data.table)
library(dplyr)
library(here)

here::i_am("01_data_preprocessing/preprocess_and_define_loci.R")

dir.create(here("data"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("results"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("03_colocalisation"), showWarnings = FALSE, recursive = TRUE)

#Step 1: mQTL Preprocessing 

bim_all <- rbindlist(lapply(1:22, function(chr) {
  fread(
    here("references", "ld_reference", paste0("EUR_phase3_chr", chr, ".bim")),
    col.names = c("CHR", "SNP", "CM", "POS", "A1_ref", "A2_ref")
  )
}))
cat("Total SNPs in reference panel:", nrow(bim_all), "\n")

mqtl_full <- fread(here("data", "All_Imputed_BonfSignificant_mQTLs.csv"))

mqtl_full <- mqtl_full %>% rename(
  CHR = SNP_Chr,
  POS = SNP_BP,
  EA_mqtl = SNP_Allele,
  BETA_mqtl = beta,
  P_mqtl = p.value
)

mqtl_full <- inner_join(
  mqtl_full,
  bim_all %>% dplyr::select(CHR, POS, SNP, A1_ref, A2_ref),
  by = c("CHR", "POS")
)
cat("mQTL SNPs after adding rsIDs:", nrow(mqtl_full), "\n")

mqtl_full <- mqtl_full %>%
  mutate(
    Z_mqtl = qnorm(P_mqtl / 2, lower.tail = FALSE),
    SE_mqtl = abs(BETA_mqtl) / Z_mqtl
  ) %>%
  filter(is.finite(SE_mqtl), SE_mqtl > 0)

mqtl_full <- mqtl_full %>%
  group_by(SNP) %>%
  slice_min(P_mqtl, n = 1) %>%
  ungroup()
cat("Final mQTL SNPs after QC:", nrow(mqtl_full), "\n")

mqtl_full <- mqtl_full %>% mutate(N_mqtl = 166)

#Step 2: GWAS QC

gwas <- fread(
  here("data", "PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv"),
  skip = "CHROM"
)
cat("Raw GWAS SNPs loaded:", nrow(gwas), "\n")

gwas <- gwas %>% rename(
  CHR = CHROM,
  SNP = ID,
  EA = A1,
  OA = A2,
  BETA_gwas = BETA,
  SE_gwas = SE,
  P_gwas = PVAL
)

gwas <- gwas %>% filter(IMPINFO >= 0.8)
cat("SNPs after IMPINFO filter:", nrow(gwas), "\n")

gwas <- gwas %>%
  mutate(MAF = pmin(FCAS, 1 - FCAS)) %>%
  filter(MAF >= 0.01)
cat("SNPs after MAF filter:", nrow(gwas), "\n")

gwas <- gwas %>%
  filter(!(EA == "A" & OA == "T"),
         !(EA == "T" & OA == "A"),
         !(EA == "C" & OA == "G"),
         !(EA == "G" & OA == "C"))
cat("SNPs after removing palindromic:", nrow(gwas), "\n")

gwas <- gwas %>%
  filter(!duplicated(SNP)) %>%
  filter(SE_gwas > 0, abs(BETA_gwas) < 10)
cat("Final GWAS SNPs after QC:", nrow(gwas), "\n")

fwrite(gwas, here("data", "gwas_qc.csv"))
cat("Clean GWAS saved\n")

gwas_maf <- gwas %>% dplyr::select(SNP, MAF)
mqtl_full <- mqtl_full %>% left_join(gwas_maf, by = "SNP")
cat("mQTL SNPs with MAF:", sum(!is.na(mqtl_full$MAF)), "\n")
cat("mQTL SNPs without MAF:", sum(is.na(mqtl_full$MAF)), "\n")

fwrite(mqtl_full, here("data", "mqtl_preprocessed.csv"))
cat("Preprocessed mQTL saved\n")

#Step 3: Locus Definition

gwas_sig <- gwas %>% filter(P_gwas < 5e-8)
cat("Genome-wide significant SNPs:", nrow(gwas_sig), "\n")

define_loci_fast <- function(gwas_sig, window = 1000000) {
  df <- gwas_sig %>% arrange(P_gwas)
  lead_snps <- data.frame()
  while (nrow(df) > 0) {
    top <- df[1, ]
    lead_snps <- rbind(lead_snps, top)
    df <- df %>% filter(!(CHR == top$CHR & abs(POS - top$POS) <= window / 2))
  }
  lead_snps
}

lead_snps <- define_loci_fast(gwas_sig)
cat("Independent loci defined:", nrow(lead_snps), "\n")

fwrite(lead_snps, here("results", "lead_snps.csv"))

#Step 4: Create COLOC-reporter Input Files

regions_spargo <- lead_snps %>%
  mutate(
    start = ifelse(POS - 500000 < 0, 0, POS - 500000),
    end = POS + 500000,
    traits = "SCZ,mQTL",
    region = paste(CHR, start, end, sep = ",")
  ) %>%
  dplyr::select(traits, region)

write.table(
  regions_spargo,
  here("03_colocalisation", "set.regions.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
cat("set.regions.txt created with", nrow(regions_spargo), "regions\n")

gwas_samples <- data.frame(
  ID = c("SCZ", "mQTL"),
  type = c("cc", "quant"),
  prop = c(0.5, "NA"),
  traitSD = c("NA", "NA"),
  p_col = c("P_gwas", "P_mqtl"),
  stat_col = c("BETA_gwas", "BETA_mqtl"),
  N_col = c("NEFF", "N_mqtl"),
  chr_col = c("CHR", "CHR"),
  pos_col = c("POS", "POS"),
  se_col = c("SE_gwas", "SE_mqtl"),
  snp_col = c("SNP", "SNP"),
  A1_col = c("EA", "EA_mqtl"),
  A2_col = c("OA", "A2_ref"),
  freq_col = c("MAF", "MAF"),
  traitLabel = c("Schizophrenia", "Fetal_Brain_mQTL"),
  FILEPATH = c(
    here("data", "gwas_qc.csv"),
    here("data", "mqtl_preprocessed.csv")
  )
)

write.table(
  gwas_samples,
  here("03_colocalisation", "GWAS_samples.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
cat("GWAS_samples.txt created\n")



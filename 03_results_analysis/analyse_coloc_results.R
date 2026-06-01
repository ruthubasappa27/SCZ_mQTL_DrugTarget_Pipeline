
# SCZ mQTL Drug Target Prioritisation Pipeline
# Script 2: Results Analysis
# Author: Ruthushree P Basappa
# Institution: University of Galway
# Supervisor: Prof. Derek Morris
# Date: May 2026


library(data.table)
library(dplyr)
library(httr)
library(jsonlite)
library(here)

here::i_am("04_results_analysis/analyse_coloc_results.R")

dir.create(here("results"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("results", "credible_sets"), showWarnings = FALSE, recursive = TRUE)

#Step 1: Load Colocalisation Results 

all_results <- rbindlist(
  lapply(
    list.files(
      here("results", "coloc"),
      pattern = "results_summary_coloc_abf.csv",
      recursive = TRUE,
      full.names = TRUE
    ),
    fread
  ),
  fill = TRUE
)

all_results_unique <- all_results[!duplicated(region)]
cat("Total unique loci tested:", nrow(all_results_unique), "\n")
cat("Colocalised hits PP4 > 0.8:", sum(all_results_unique$PP.H4.abf > 0.8, na.rm = TRUE), "\n")

hits <- all_results_unique[PP.H4.abf > 0.8][order(-PP.H4.abf)]
print(hits[, .(region, nsnps, PP.H4.abf)])

fwrite(hits, here("results", "coloc_hits_pp4_gt_0.8.csv"))

#Step 2: Credible Set Construction

locus_folders <- c(
  "1_28532580_29532580",
  "2_212685645_213685645",
  "6_27859632_28859632",
  "8_9532894_10532894",
  "5_139833952_140833952",
  "2_172456449_173456449"
)

gene_names <- c("OPRD1", "ERBB4", "ZNF389", "PIK3C2B", "WDR55", "DYNC1I2")

top_snps <- c("rs61783570", "rs7607363", "rs3118359",
              "rs4633059", "rs2530240", "rs4667693")

for (i in seq_along(locus_folders)) {
  snpwise_file <- here(
    "results", "coloc",
    paste0(locus_folders[i], "_coloc"),
    "tables",
    "coloc.abf_snpwise_PP_H4_abf.csv"
  )

  if (!file.exists(snpwise_file)) next

  snpwise <- fread(snpwise_file)
  snpwise <- snpwise[order(-SNP.PP.H4)]
  snpwise[, cumPP := cumsum(SNP.PP.H4) / sum(SNP.PP.H4)]
  credible_set <- snpwise[cumPP <= 0.95]

  fwrite(
    credible_set,
    here("results", "credible_sets", paste0(gene_names[i], "_credible_set.csv"))
  )

  cat("=== Gene:", gene_names[i], "| Locus:", locus_folders[i], "===\n")
  cat("95% credible set size:", nrow(credible_set), "SNPs\n")
  cat("Top SNP:", snpwise$snp[1], "| SNP.PP.H4:", round(snpwise$SNP.PP.H4[1], 4), "\n\n")
}

#Step 3: Gene Annotation

mqtl <- fread(here("data", "mqtl_preprocessed.csv"))

for (i in seq_along(top_snps)) {
  snp_info <- mqtl %>%
    dplyr::filter(SNP == top_snps[i]) %>%
    dplyr::select(SNP, CHR, POS, ProbeID, DNAm_CHR, DNAm_BP, GeneAnnotation)

  cat("=== Gene:", gene_names[i], "| SNP:", top_snps[i], "===\n")
  print(snp_info)
  cat("\n")
}

#Step 4: Variant Location Annotation

get_snp_annotation <- function(rsid) {
  url <- paste0(
    "https://myvariant.info/v1/variant/", rsid,
    "?assembly=hg19&fields=cadd.gene,cadd.consequence,dbsnp.gene"
  )
  response <- GET(url)
  if (status_code(response) == 200) {
    return(fromJSON(content(response, "text", encoding = "UTF-8")))
  }
  NULL
}

for (i in seq_along(top_snps)) {
  cat("=== Gene:", gene_names[i], "| SNP:", top_snps[i], "===\n")
  result <- get_snp_annotation(top_snps[i])
  if (!is.null(result)) {
    tryCatch({
      cat("Gene:", result$cadd$gene$genename, "\n")
      cat("Consequence:", result$cadd$consequence, "\n")
    }, error = function(e) cat("Could not parse annotation\n"))
  }
  cat("\n")
}

#Step 5: Directionality Analysis

for (i in seq_along(top_snps)) {
  harm_file <- here(
    "results", "coloc",
    paste0(locus_folders[i], "_coloc"),
    "data", "datasets", "harmonised_sumstats.csv"
  )

  if (!file.exists(harm_file)) next

  harm <- fread(harm_file)
  dir_result <- harm %>%
    dplyr::filter(snp == top_snps[i]) %>%
    dplyr::select(snp, trait, beta, pvalues, cs)

  scz_beta <- dir_result$beta[dir_result$trait == "SCZ"]
  mqtl_beta <- dir_result$beta[dir_result$trait == "mQTL"]

  cat("=== Gene:", gene_names[i], "| SNP:", top_snps[i], "===\n")
  print(dir_result)

  if (length(scz_beta) > 0 && length(mqtl_beta) > 0) {
    if (sign(scz_beta) == sign(mqtl_beta)) {
      cat("Direction: CONCORDANT\n")
    } else {
      cat("Direction: DISCORDANT\n")
    }
  }
  cat("\n")
}


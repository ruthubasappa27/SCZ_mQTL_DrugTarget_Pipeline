# SCZ mQTL Drug Target Prioritisation Pipeline

## Overview

This repository contains the computational pipeline used to prioritise candidate therapeutic targets for schizophrenia (SCZ) by integrating genome-wide association study (GWAS) summary statistics with fetal brain methylation quantitative trait loci (mQTL) data.

The project uses SCZ-fetal brain mQTL colocalisation as the primary statistical analysis, followed by functional annotation, causal-support analysis, and drug-target prioritisation. The workflow is designed to produce a conservative, evidence-based shortlist of candidate genes for downstream translational follow-up.

**Author:** Ruthushree P Basappa  
**Institution:** University of Galway  
**Supervisor:** Prof. Derek Morris  
**Degree:** MSc Clinical Neuroscience  

## Study design

- **GWAS:** PGC3 schizophrenia GWAS (Trubetskoy et al. 2022), n=127,906
- **mQTL:** Hannon et al. fetal brain mQTL (2016), n=166
- **Genome build:** GRCh37 / hg19
- **Primary colocalisation framework:** COLOC-reporter (Spargo et al. 2024)

## Key results

Six loci showed strong evidence of colocalisation between SCZ GWAS and fetal brain mQTL signals using a PP4 threshold of 0.8.

| Gene | PP4 | p_SMR | p_HEIDI |
|------|-----:|------:|--------:|
| OPRD1 | 0.993 | 1.19e-05 | NA |
| ERBB4 | 0.990 | 5.61e-07 | 0.108 |
| ZNF389 | 0.945 | NA | NA |
| PIK3C2B | 0.913 | 1.13e-04 | 0.285 |
| WDR55 | 0.881 | 3.51e-05 | 0.604 |
| DYNC1I2 | 0.816 | 1.42e-04 | 0.825 |

HEIDI was used as a heterogeneity check in the SMR framework, and non-significant HEIDI values were treated as supportive of a shared underlying signal rather than definitive proof of causality.


## Repository structure

```text
SCZ_mQTL_DrugTarget_Pipeline/
├── data/
│   └── README.md
├── references/
│   └── README.md
├── results/
│   └── README.md
├── 01_data_preprocessing/
│   └── preprocess_and_define_loci.R
├── 03_colocalisation/
│   ├── GWAS_samples.txt
│   └── run_coloc_all.sh
├── 04_results_analysis/
│   └── analyse_coloc_results.R
├── 05_functional_validation/
│   └── cadd_annotation.R
├── 06_SMR_analysis/
│   ├── all_probes.txt
│   ├── prepare_gwas_smr.R
│   └── run_smr.sh
└── README.md
```

## Requirements

### Software

- R 4.5.1
- PLINK v1.90
- SMR v1.03

### R packages

- data.table
- dplyr
- coloc
- susieR
- httr
- jsonlite

### Operating system

- Tested on macOS Apple Silicon
- Linux-compatible with the appropriate SMR binary

## Data sources

All datasets used in this pipeline are publicly available.

1. **PGC3 SCZ GWAS:** [PGC download page](https://pgc.unc.edu/for-researchers/download-results/)
2. **Hannon fetal brain mQTL:** [Essex mQTL resource](https://epigenetics.essex.ac.uk/mQTL/) and [PubMed article](https://pubmed.ncbi.nlm.nih.gov/26619357/)
3. **Hannon BESD-format mQTL:** [SMR data resources](https://yanglab.westlake.edu.cn/software/smr/#DataResource)
4. **1000 Genomes Phase 3 EUR LD reference panel:** [PLINK format](https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe), [SMR format](https://yanglab.westlake.edu.cn/software/smr/)
5. **COLOC-reporter pipeline:** [GitHub repository](https://github.com/ThomasPSpargo/COLOC-reporter)

## Workflow

### 1. Data preprocessing

Run `01_data_preprocessing/preprocess_and_define_loci.R`.

This step performs:
- GWAS QC, including filtering by imputation quality, MAF, and palindromic SNP removal
- mQTL preprocessing, including column standardisation, standard error calculation, and rsID assignment
- Export of cleaned GWAS and mQTL summary-statistics files

**Outputs:**
- `gwas_qc.csv`
- `mqtl_preprocessed.csv`

### 2. Locus definition

This step is performed within `01_data_preprocessing/preprocess_and_define_loci.R`.

A distance-based greedy approach was used to define independent GWAS loci using 1 Mb windows centered on lead SNPs.

**Output:**
- `201` independent loci

### 3. Colocalisation

Run `03_colocalisation/run_coloc_all.sh`.

This step uses COLOC-reporter to test for shared causal variants at each locus.

- `coloc.susie` is used when SuSiE fine-mapping converges
- `coloc.abf` is used as a fallback when SuSiE does not converge
- Loci with PP4 >= 0.8 are considered colocalised

### 4. Results analysis

Run `04_results_analysis/analyse_coloc_results.R`.

This step includes:
- candidate credible set construction
- directionality analysis
- variant annotation
- locus-level prioritisation

### 5. Functional validation

Run `05_functional_validation/cadd_annotation.R`.

This step performs:
- CADD score annotation
- variant consequence annotation
- retrieval of functional information using the myvariant.info API

**Output:**
- `credible_set_CADD_scores_complete.csv`

### 6. SMR analysis

SMR was used as an orthogonal analysis to evaluate whether the SCZ GWAS signal and fetal brain mQTL signal at each colocalised locus are consistent with a shared underlying variant.

Run:

```bash
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
```

The Hannon fetal brain mQTL BESD file uses chr:bp SNP identifiers rather than rsIDs. GWAS summary statistics and the LD reference panel were converted from rsID to chr:bp format before SMR analysis to ensure consistent SNP matching.

## Notes

- This repository is under active development
- The current workflow prioritises specificity and reproducibility over maximal sensitivity
- Conservative harmonisation choices may exclude some borderline variants, but retained loci represent high-confidence candidates
- This repository is intended for computational prioritisation and does not claim experimental validation of drug targets

## References

1. Trubetskoy V, et al. Nature. 2022;604:502-508. [doi](https://doi.org/10.1038/s41586-022-04434-5)
2. Hannon E, et al. Nature Neuroscience. 2016;19:48-54. [doi](https://doi.org/10.1038/nn.4182)
3. Spargo TP, et al. eLife. 2024;12:RP88768. [doi](https://doi.org/10.7554/eLife.88768)
4. Giambartolomei C, et al. PLoS Genetics. 2014;10:e1004383. [doi](https://doi.org/10.1371/journal.pgen.1004383)
5. Wallace C. PLoS Genetics. 2021;17:e1009440. [doi](https://doi.org/10.1371/journal.pgen.1009440)
6. Zhu Z, et al. Nature Genetics. 2016;48:481-487. [doi](https://doi.org/10.1038/ng.3538)
7. Wang G, et al. J R Stat Soc Series B. 2020;82:1273-1300. [doi](https://doi.org/10.1111/rssb.12388)
8. 1000 Genomes Project Consortium. Nature. 2015;526:68-74. [doi](https://doi.org/10.1038/nature15393)
9. Rentzsch P, et al. Nucleic Acids Research. 2019;47:D886-D894. [doi](https://doi.org/10.1093/nar/gky1016)

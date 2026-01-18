# protein-biomarkers-for-bile-duct-cancer

# Data Folder

## 1. Raw input data

(1) Instrumental variables of plasma proteins used in MR analysis

The reference for the instrumental variables of plasma proteins are (DOI: 10.1186/s13073-023-01229-9)
The raw instrumental variables of plasma proteins is "Raw_input_data/Plasma.Protein_Instruments.Raw.xls"

The cis-pQTLs were selected based on the following criteria:
(a). SNPs linked to any protein with a significance level of P < 5×10-8 were included;
(b). SNPs and proteins located within the Major Histocompatibility Complex (MHC) region on chromosome 6 (25.5–34.0 Mb) were excluded because of the intricate linkage disequilibrium (LD) patterns in this area;
(c). subsequently, LD clumping was performed to identify independent protein quantitative trait loci (pQTLs) for each protein, ensuring that variants had an r² value less than 0.001;
(d). F-statistics genetic instruments were higher than 10

Ultimately, the instrumental variables file used for analysis is "Raw_input_data/Instrumental variables of plasma proteins used in MR analysis.xls"

(2) Three GWAS summary data of Bile tract cancer (BTC) patients

Three GWAS summary data of BTC patients and controls of European ancestry from three public cohorts.

<img width="2208" height="295" alt="image" src="https://github.com/user-attachments/assets/a80b892b-00ba-44d7-b142-86b8c1a9234e" />


Notes:

(a). The raw C3_BILIARY_GALLBLADDER_EXALLC data link: https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_C3_BILIARY_GALLBLADDER_EXALLC.gz

(b). The raw ieu-b-4915 VCF data link: https://opengwas.io/datasets/ieu-b-4915#

(c). The raw GCST90043859 data link: https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90043001-GCST90044000/GCST90043859/GCST90043859_buildGRCh37.tsv.gz

## 2. intermediate files

After obtaining the GWAS data and plasma protein instrumental variables, we extracted the SNP information of the exposure and outcome, harmonized the data, and ultimately ensured that the effect alleles of SNPs in the exposure and outcome data were consistent. 

The three obtained datasets are as follows：

(a). "intermediate files/Harmonising_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer.xls": C3_BILIARY_GALLBLADDER_EXALLC harmonized with  "Raw_input_data/Instrumental variables of plasma proteins used in MR analysis.xls"

(c). "intermediate files/Harmonising_exposure.cis-pQTLs_protein_outcome.ieu4915.bile.duct.cancer.xls": ieu-b-4915 harmonized with  "Raw_input_data/Instrumental variables of plasma proteins used in MR analysis.xls"

(c). "intermediate files/Harmonising_exposure.cis-pQTLs_protein_outcome.GCST90043859.bile.duct.cancer.xls": GCST90043859 harmonized with  "Raw_input_data/Instrumental variables of plasma proteins used in MR analysis.xls"

## 3.final results

### final results/Discovery cohort

"mr_result_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer_addOR.xls" : The results of the MR analysis for the Discovery cohort(C3_BILIARY_GALLBLADDER_EXALLC)

"mr_result_exposure.cis-PlasmaProtein_outcome.finn12_steiger_direction.xls" : The results of the Steiger test for the Discovery cohort(C3_BILIARY_GALLBLADDER_EXALLC)

### final results/Validation cohort 1

"mr_result_exposure.cis-pQTLs_protein_outcome.ieu4915.bile.duct.cancer_addOR.xls" : The results of the MR analysis for the Validation cohort 1(ieu-b-4915)

"mr_result_exposure.cis-PlasmaProtein_ieu4915_steiger_direction.xls" : The results of the Steiger test for the Validation cohort 1(ieu-b-4915)

### final results/Validation cohort 2

"mr_result_exposure.cis-pQTLs_protein_outcome.g3859.bile.duct.cancer_addOR.xls" : The results of the MR analysis for the Validation cohort 1(GCST90043859)

"mr_result_exposure.cis-PlasmaProtein_g3859_steiger_direction.xls" : The results of the Steiger test for the Validation cohort 1(GCST90043859)

# Prerequisites

## R dependencies

### R version 4.4.3 

Users running Platform: x86_64-pc-linux-gnu and Running under: Ubuntu 20.04.6 LTS, to install the latest version of TwoSampleMR_0.5.8

install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")

other attached packages:
 [1] rio_1.2.4                     reticulate_1.44.1             BiocParallel_1.40.2           patchwork_1.3.2                     
 [5] GenomicRanges_1.58.0          GenomeInfoDb_1.42.3           IRanges_2.40.1                MatrixGenerics_1.18.1                 
 [9] matrixStats_1.5.0             Rcpp_1.1.0                    labeling_0.4.3                lubridate_1.9.4     
[13] forcats_1.0.1                 stringr_1.5.2                 purrr_1.1.0                   tidyr_1.3.1                 
[17] tibble_3.3.0                  tidyverse_2.0.0               cowplot_1.2.0                 ggraph_2.2.2   
[21] irlba_2.3.5.1                  VGAM_1.1-14                  Biobase_2.66.0                Matrix_1.6-5  
[25] sp_2.2-0                      forestplot_3.1.7              abind_1.4-8                   ggpubr_0.6.2  
[29] ggsci_4.0.0                   ggplot2_4.0.0                 TwoSampleMR_0.5.8             readr_2.1.6    
[33] S4Vectors_0.44.0              BiocGenerics_0.52.0           plinkbinr_0.0.0.9000          MendelianRandomization_0.10.0                  
[37] gassocplot_0.0.2              gwasglue_0.0.0.9000           ieugwasr_1.0.1                gwasvcf_0.1.2                                
[41] coloc_5.2.3                   locuscomparer_1.0.0           MRPRESSO_1.0                  plyr_1.8.9                      
[45] readxl_1.4.5                  data.table_1.17.8             dplyr_1.1.4




















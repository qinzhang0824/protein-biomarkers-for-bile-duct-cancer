library(TwoSampleMR)
library(readr)
library(dplyr)
library(data.table)
library(readxl)
library(plyr)
library(MRPRESSO)
library("locuscomparer")
library(coloc)
library(gwasglue)
library(gassocplot)
#library(export)
library(TwoSampleMR)
library(data.table)
library(MendelianRandomization)
library(ieugwasr)
library(plinkbinr)
library(S4Vectors)
library(rio)
library(ggplot2)
library(ggsci)
library("ggpubr")
library(forestplot)

########################################################################################### MR analysis finngen.R12
mydata <- fread('intermediate_files/Harmonising_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer.txt',header = T)

result <- mr(mydata, method_list=c("mr_egger_regression", "mr_ivw","mr_weighted_median","mr_wald_ratio"))
MR.OR<-generate_odds_ratios(result)

export(MR.OR,"Final_results/Discovery_cohort/mr_result_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer_addOR.txt",format = "\t") 

mr_steiger_direction <-directionality_test(mydata)
export(mr_steiger_direction,"Final_results/Discovery_cohort/mr_result_exposure.cis-PlasmaProtein_outcome.finn12_steiger_direction.txt",format = "\t")

###########################################################################################  IEU-b-4915
mydata <- fread('intermediate_files/Harmonising_exposure.cis-pQTLs_protein_outcome.ieu4915.bile.duct.cancer.xls',header = T)

result <- mr(mydata, method_list=c("mr_egger_regression", "mr_ivw","mr_weighted_median","mr_wald_ratio"))
MR.OR<-generate_odds_ratios(result)

export(MR.OR,"Final_results/Validation_cohort_1/mr_result_exposure.cis-pQTLs_protein_outcome.ieu4915.bile.duct.cancer_addOR.xls",format = "\t") 

mr_steiger_direction <-directionality_test(mydata)
export(mr_steiger_direction,"Final_results/Validation_cohort_1/mr_result_exposure.cis-PlasmaProtein_ieu4915_steiger_direction.xls",format = "\t")

############################################################################################## GCST90043859
mydata <- fread('intermediate_files/Harmonising_exposure.cis-pQTLs_protein_outcome.GCST90043859.bile.duct.cancer.xls',header = T)

result <- mr(mydata, method_list=c("mr_egger_regression", "mr_ivw","mr_weighted_median","mr_wald_ratio"))
MR.OR<-generate_odds_ratios(result)

export(MR.OR,"Final_results/Validation_cohort_2/mr_result_exposure.cis-pQTLs_protein_outcome.g3859.bile.duct.cancer_addOR.xls",format = "\t") 

mr_steiger_direction <-directionality_test(mydata)
export(mr_steiger_direction,"Final_results/Validation_cohort_2/mr_result_exposure.cis-PlasmaProtein_g3859_steiger_direction.xls",format = "\t")


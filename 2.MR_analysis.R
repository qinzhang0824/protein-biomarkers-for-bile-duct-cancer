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

.libPaths(c('~/R/x86_64-pc-linux-gnu-library/4.3',
            '/refdir/Rlib_4.3',
            '/usr/local/lib/R/library'))


mydata <- fread('Harmonising_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer.xls',header = T)

result <- mr(mydata, method_list=c("mr_egger_regression", "mr_ivw","mr_weighted_median","mr_wald_ratio"))
export(result,"mr_result_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer.xls",format = "\t") 

MR.OR<-generate_odds_ratios(result)
export(MR.OR,"mr_result_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer_addOR.xls",format = "\t") 

MR.OR<-fread('mr_result_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer_addOR.xls',header = T)

het <- mr_heterogeneity(mydata)
export(het,"mr_result_exposure.cis-PlasmaProtein_outcome.finn12_het.xls",format = "\t")

res_single <- mr_singlesnp(mydata,parameters = default_parameters(),
                           single_method = "mr_wald_ratio",
                           all_method = c("mr_ivw", "mr_egger_regression"))
export(res_single,"mr_result_exposure.cis-PlasmaProtein_outcome.finn12_singleSNP.xls",format = "\t")

pleio <- mr_pleiotropy_test(mydata)
export(pleio,"mr_result_exposure.cis-PlasmaProtein_outcome.finn12_pleio.xls",format = "\t")

mr_steiger_direction <-directionality_test(mydata)
export(mr_steiger_direction,"mr_result_exposure.cis-PlasmaProtein_outcome.finn12_steiger_direction.xls",format = "\t")

single <- mr_leaveoneout(mydata)
export(single,"mr_result_exposure.cis-PlasmaProtein_outcome.finn12_Leave-one-out.xls",format = "\t")

###########################################################################################  IEU-b-4915
mydata <- fread('Harmonising_exposure.cis-pQTLs_protein_outcome.ieu4915.bile.duct.cancer.xls',header = T)

result <- mr(mydata, method_list=c("mr_egger_regression", "mr_ivw","mr_weighted_median","mr_wald_ratio"))
export(result,"mr_result_exposure.cis-pQTLs_protein_outcome.ieu4915.bile.duct.cancer.xls",format = "\t") 

MR.OR<-generate_odds_ratios(result)
export(MR.OR,"mr_result_exposure.cis-pQTLs_protein_outcome.ieu4915.bile.duct.cancer_addOR.xls",format = "\t") 
MR.OR<-fread('mr_result_exposure.cis-pQTLs_protein_outcome.ieu4915.bile.duct.cancer_addOR.xls',header = T)

het <- mr_heterogeneity(mydata)
export(het,"mr_result_exposure.cis-PlasmaProtein_ieu4915_het.xls",format = "\t")

res_single <- mr_singlesnp(mydata,parameters = default_parameters(),
                           single_method = "mr_wald_ratio",
                           all_method = c("mr_ivw", "mr_egger_regression"))
export(res_single,"mr_result_exposure.cis-PlasmaProtein_ieu4915_singleSNP.xls",format = "\t")

pleio <- mr_pleiotropy_test(mydata)
export(pleio,"mr_result_exposure.cis-PlasmaProtein_ieu4915_pleio.xls",format = "\t")

mr_steiger_direction <-directionality_test(mydata)
export(mr_steiger_direction,"mr_result_exposure.cis-PlasmaProtein_ieu4915_steiger_direction.xls",format = "\t")

single <- mr_leaveoneout(mydata)
export(single,"mr_result_exposure.cis-PlasmaProtein_ieu4915_Leave-one-out.xls",format = "\t")

############################################################################################################### GCST90043859
mydata <- fread('Harmonising_exposure.cis-pQTLs_protein_outcome.GCST90043859.bile.duct.cancer.xls',header = T)

result <- mr(mydata, method_list=c("mr_egger_regression", "mr_ivw","mr_weighted_median","mr_wald_ratio"))

result.t <- mr(mydata.T, method_list=c("mr_egger_regression", "mr_ivw","mr_weighted_median","mr_wald_ratio"))
export(result,"mr_result_exposure.cis-pQTLs_protein_outcome.g3859.bile.duct.cancer.xls",format = "\t")

MR.OR<-generate_odds_ratios(result)
MR.OR.t<-generate_odds_ratios(result.t)

export(MR.OR,"mr_result_exposure.cis-pQTLs_protein_outcome.g3859.bile.duct.cancer_addOR.xls",format = "\t") 
MR.OR<-fread('mr_result_exposure.cis-pQTLs_protein_outcome.g3859.bile.duct.cancer_addOR.xls',header = T)

het <- mr_heterogeneity(mydata)
export(het,"mr_result_exposure.cis-PlasmaProtein_g3859_het.xls",format = "\t")

res_single <- mr_singlesnp(mydata,parameters = default_parameters(),
                           single_method = "mr_wald_ratio",
                           all_method = c("mr_ivw", "mr_egger_regression"))
export(res_single,"mr_result_exposure.cis-PlasmaProtein_g3859_singleSNP.xls",format = "\t")

pleio <- mr_pleiotropy_test(mydata)
export(pleio,"mr_result_exposure.cis-PlasmaProtein_g3859_pleio.xls",format = "\t")

mr_steiger_direction <-directionality_test(mydata)
export(mr_steiger_direction,"mr_result_exposure.cis-PlasmaProtein_g3859_steiger_direction.xls",format = "\t")

single <- mr_leaveoneout(mydata)
export(single,"mr_result_exposure.cis-PlasmaProtein_g3859_Leave-one-out.xls",format = "\t")

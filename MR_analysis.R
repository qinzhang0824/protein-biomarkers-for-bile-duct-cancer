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

setwd('/home/data/t020412/Qinzhang_data/MR_drug.target_BileDuctCancer')

#########Liver & bile duct cancer
vcf <- VariantAnnotation::readVcf("ieu-b-4915.vcf", "hg19")
Bile_outcome <- gwasvcf_to_TwoSampleMR(vcf,type = "outcome")
export(Bile_outcome,"ieu-b-4915_outcome.MR.format.xls",format = "\t")
outcome_dat <- fread ("ieu-b-4915_outcome.MR.format.xls",sep='\t',header = T)

#########Biliary tract cancer East Asian
vcf <- VariantAnnotation::readVcf("bbj-a-92.vcf", "hg19")
Bile_outcome <- gwasvcf_to_TwoSampleMR(vcf,type = "outcome")
export(Bile_outcome,"bbj-a-92_outcome.MR.format.xls",format = "\t")
outcome_dat <- fread ("bbj-a-92_outcome.MR.format.xls",sep='\t',header = T)



#############血浆蛋白数据
#############血浆蛋白数据,排除掉sun_2 和Yao
protein <- fread("Plasma.Protein_Instruments.Ref_no_Sun.2_Yao.xls",sep='\t')
protein$phenotype <- "Proteins"
protein <- as.data.frame(protein)
protein.mr <- format_data(protein,type = "exposure",
                          header = TRUE,
                          phenotype_col = "phenotype",
                          snp_col = "SNP",
                          beta_col = "Beta",
                          se_col = "Se",
                          eaf_col = "Eaf",
                          effect_allele_col = "Effect allele",
                          other_allele_col = "Other allele",
                          pval_col = "P",
                          chr_col = "Chr",
                          samplesize_col = "N",
                          id_col = "Protein",
                          pos_col = "Pos")

#############血浆蛋白数据,排除掉sun_2 和Yao 挑出cis-pQTLs 总计1178个pro,1193个IVs

p.sig <- 0.05/1178

protein <- fread("Plasma.Protein_Instruments.Ref_no_Sun.2_Yao_cis.xls",sep='\t')
protein$phenotype <- "Proteins"
protein <- as.data.frame(protein)
protein.mr <- format_data(protein,type = "exposure",
                          header = TRUE,
                          phenotype_col = "phenotype",
                          snp_col = "SNP",
                          beta_col = "Beta",
                          se_col = "Se",
                          eaf_col = "Eaf",
                          effect_allele_col = "Effect allele",
                          other_allele_col = "Other allele",
                          pval_col = "P",
                          chr_col = "Chr",
                          samplesize_col = "N",
                          id_col = "Protein",
                          pos_col = "Pos")

################################################################################################################### 胆管癌FinnGen
finngen <-fread('finngen_R12_C3_BILIARY_GALLBLADDER_EXALLC',header = T,sep='\t')
finngen <-as.data.frame(finngen)
finngen$phenotype <- "Bile.duct.cancer"
finngen$samplesize <- as.numeric(381047)
bile.finn <- format_data(finngen,type = "outcome",
                         header = TRUE,
                         phenotype_col = "phenotype",
                         snp_col = "rsids",
                         beta_col = "beta",
                         se_col = "sebeta",
                         eaf_col = "af_alt_cases",
                         effect_allele_col = "alt",
                         other_allele_col = "ref",
                         pval_col = "pval",
                         chr_col = "#chrom",
                         samplesize_col = "samplesize",
                         id_col = "phenotype",
                         pos_col = "pos")

export(bile.finn,"finngen_R12_C3_Bile_outcome.MR.format.xls",format = "\t")
bile.finn <- fread ("finngen_R12_C3_Bile_outcome.MR.format.xls",sep='\t',header = T)

##################################################################################################### MR 分析cis-pQTL VS finngen R12

mydata <- harmonise_data(exposure_dat=protein.mr,outcome_dat=bile.finn,action= 2)

export(mydata,"Harmonising_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer.xls",format = "\t")

mydata <- fread('Harmonising_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer.xls',header = T)

result <- mr(mydata, method_list=c("mr_egger_regression", "mr_ivw","mr_weighted_median","mr_wald_ratio"))
export(result,"mr_result_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer.xls",format = "\t") 

MR.OR<-generate_odds_ratios(result)
export(MR.OR,"mr_result_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer_addOR.xls",format = "\t") 

MR.OR<-fread('mr_result_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer_addOR.xls',header = T)
##################################################################################################2025.03.03
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
ieu4915<-outcome_dat
mydata <- harmonise_data(exposure_dat=protein.mr,outcome_dat=ieu4915,action= 2)

export(mydata,"Harmonising_exposure.cis-pQTLs_protein_outcome.ieu4915.bile.duct.cancer.xls",format = "\t")
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


####################################################################################################  bbj-a-92
bbj92<-outcome_dat


mydata <- harmonise_data(exposure_dat=protein.mr,outcome_dat=bbj92,action= 2)

export(mydata,"Harmonising_exposure.cis-pQTLs_protein_outcome.bbj92.bile.duct.cancer.xls",format = "\t")
mydata <- fread('Harmonising_exposure.cis-pQTLs_protein_outcome.bbj92.bile.duct.cancer.xls',header = T)

result <- mr(mydata, method_list=c("mr_egger_regression", "mr_ivw","mr_weighted_median","mr_wald_ratio"))
export(result,"mr_result_exposure.cis-pQTLs_protein_outcome.bbj92.bile.duct.cancer.xls",format = "\t")

MR.OR<-generate_odds_ratios(result)
export(MR.OR,"mr_result_exposure.cis-pQTLs_protein_outcome.bbj92.bile.duct.cancer_addOR.xls",format = "\t") 
MR.OR<-fread('mr_result_exposure.cis-pQTLs_protein_outcome.bbj92.bile.duct.cancer_addOR.xls',header = T)

het <- mr_heterogeneity(mydata)
export(het,"mr_result_exposure.cis-PlasmaProtein_bbj92_het.xls",format = "\t")

res_single <- mr_singlesnp(mydata,parameters = default_parameters(),
                           single_method = "mr_wald_ratio",
                           all_method = c("mr_ivw", "mr_egger_regression"))
export(res_single,"mr_result_exposure.cis-PlasmaProtein_bbj92_singleSNP.xls",format = "\t")

pleio <- mr_pleiotropy_test(mydata)
export(pleio,"mr_result_exposure.cis-PlasmaProtein_bbj92_pleio.xls",format = "\t")

mr_steiger_direction <-directionality_test(mydata)
export(mr_steiger_direction,"mr_result_exposure.cis-PlasmaProtein_bbj92_steiger_direction.xls",format = "\t")

single <- mr_leaveoneout(mydata)
export(single,"mr_result_exposure.cis-PlasmaProtein_bbj92_Leave-one-out.xls",format = "\t")

##############################################################################################################GCST90013662
g3662 <-fread('GCST90013662_buildGRCh37.tsv',header = T,sep='\t')
g3662 <-as.data.frame(g3662)
g3662$phenotype <- "BileDuct.cancer"
g3662$samplesize <- as.numeric(408183)
g3662MR <- format_data(g3662,type = "outcome",
                       header = TRUE,
                       phenotype_col = "phenotype",
                       snp_col = "rsid",
                       beta_col = "beta_SAIGE",
                       se_col = "se_SAIGE",
                       #eaf_col = "af_alt",
                       effect_allele_col = "effect_allele_SAIGE",
                       other_allele_col = "other_allele_SAIGE",
                       pval_col = "pvalue_SAIGE",
                       chr_col = "chromosome",
                       samplesize_col = "samplesize",
                       pos_col = "base_pair_location")

export(g3662MR,"BileDuct.cancer_g3662MR.MRinput.txt",format = "\t") 

mydata <- harmonise_data(exposure_dat=protein.mr,outcome_dat=g3662MR,action= 2)

result <- mr(mydata, method_list=c("mr_egger_regression", "mr_ivw","mr_weighted_median","mr_wald_ratio"))
export(result,"mr_result_exposure.cis-pQTLs_protein_outcome.g3662.bile.duct.cancer.xls",format = "\t")

MR.OR<-generate_odds_ratios(result)
export(MR.OR,"mr_result_exposure.cis-pQTLs_protein_outcome.g3662.bile.duct.cancer_addOR.xls",format = "\t") 
MR.OR<-fread('mr_result_exposure.cis-pQTLs_protein_outcome.g3662.bile.duct.cancer_addOR.xls',header = T)

het <- mr_heterogeneity(mydata)
export(het,"mr_result_exposure.cis-PlasmaProtein_g3662_het.xls",format = "\t")

res_single <- mr_singlesnp(mydata,parameters = default_parameters(),
                           single_method = "mr_wald_ratio",
                           all_method = c("mr_ivw", "mr_egger_regression"))
export(res_single,"mr_result_exposure.cis-PlasmaProtein_g3662_singleSNP.xls",format = "\t")

pleio <- mr_pleiotropy_test(mydata)
export(pleio,"mr_result_exposure.cis-PlasmaProtein_g3662_pleio.xls",format = "\t")

mr_steiger_direction <-directionality_test(mydata)
export(mr_steiger_direction,"mr_result_exposure.cis-PlasmaProtein_g3662_steiger_direction.xls",format = "\t")

single <- mr_leaveoneout(mydata)
export(single,"mr_result_exposure.cis-PlasmaProtein_g3662_Leave-one-out.xls",format = "\t")

############################################################################################################### GCST90043859
g3859 <-fread('GCST90043859_buildGRCh37.tsv',header = T,sep='\t')
g3859 <-as.data.frame(g3859)
g3859$phenotype <- "BileDuct.cancer"
#g3859$samplesize <- as.numeric(456285)
g3859MR <- format_data(g3859,type = "outcome",
                       header = TRUE,
                       phenotype_col = "phenotype",
                       snp_col = "variant_id",
                       beta_col = "beta",
                       se_col = "standard_error",
                       eaf_col = "effect_allele_frequency",
                       effect_allele_col = "effect_allele",
                       other_allele_col = "other_allele",
                       pval_col = "p_value",
                       chr_col = "chromosome",
                       samplesize_col = "N",
                       pos_col = "base_pair_location")


g3859MR.T <- format_data(g3859,type = "outcome",
                         header = TRUE,
                         phenotype_col = "phenotype",
                         snp_col = "variant_id",
                         beta_col = "T",
                         se_col = "SE_T",
                         eaf_col = "effect_allele_frequency",
                         effect_allele_col = "effect_allele",
                         other_allele_col = "other_allele",
                         pval_col = "P_noSPA",
                         chr_col = "chromosome",
                         samplesize_col = "N",
                         pos_col = "base_pair_location")


export(g3859MR,"BileDuct.cancer_g3859MR.MRinput.txt",format = "\t") 

mydata <- harmonise_data(exposure_dat=protein.mr,outcome_dat=g3859MR,action= 2)

mydata.T <- harmonise_data(exposure_dat=protein.mr,outcome_dat=g3859MR.T,action= 2)

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


##############################################################################OOOOOOOOOOOOOOOLLLLLLLLLLLLLLLLLLLLLLDDDDDDDDDDDDD
####################################################################################################################结局黑色素瘤 ieu-b-4959
################################################################################################################################结局黑色素瘤 UKB-b-12915
################################################################################################################################  UKB-d-C43
################################################################################################################################ukb-c-c3-melanoma.ski
############## 	CHMP1A  Chromosome 16: 89710839-89724253
t.chr <-16
t.start <-89710839
t.end <- 89724253


##################################################################################################################### 火山图
library(limma)
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(ggthemes)
library(TwoSampleMR)
library(data.table)

mr.or <-fread("mr_result_exposure.PlasmaProtein_outcome.FinnGen.R11.melanoma_addOR.xls",sep='\t',header=T)

mr.or <-fread("mr_result_exposure.PlasmaProtein_outcome.ieu-b-4959.melanoma_addOR_cisPQTLs.xls",sep='\t',header=T)

mr.or <-fread('mr_result_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer_addOR.xls',header = T)
mr.or <-fread('mr_result_exposure.cis-pQTLs_protein_outcome.ieu4915.bile.duct.cancer_addOR.xls',header = T)
mr.or <-fread('mr_result_exposure.cis-pQTLs_protein_outcome.g3859.bile.duct.cancer_addOR.xls',header = T)

mr.or <- mr.or[mr.or$method!="Weighted median",]
mr.or <- mr.or[mr.or$method!="MR Egger",]

0.05/1753
0.05/1178
-log10(0.05/1753)
-log10(p.sig)
-log10(0.05)
-log10(0.7372437)


res <- mr.or %>% filter(!is.na(pval)) %>% 
  mutate( logP = -log10(pval) ) %>%
  mutate( OR = or ) %>%
  mutate( tag = "Validation Bile Duct Cancer")%>%
  mutate( Gene = rownames(id.exposure))

res <- mr.or %>% filter(!is.na(pval)) %>% 
  mutate( logP = -log10(pval) ) %>%
  mutate( OR = or ) %>%
  mutate( tag = "Discovery Bile Duct Cancer")%>%
  mutate( Gene = rownames(id.exposure))


sig.res <- res %>% mutate(group=case_when(
  (pval < 0.05/1178 & or > 1) ~ "Positively associated",
  (pval < 0.05/1178 & or < 1) ~ "Negatively associated",
  .default = "Not significant"))

sig.res <- res %>% mutate(group=case_when(
  (pval < 0.05 & or > 1) ~ "Positively associated",
  (pval < 0.05 & or < 1) ~ "Negatively associated",
  .default = "Not significant"))

sig.res <- sig.res %>% mutate(group=factor(group,levels = c("Positively associated","Negatively associated","Not significant")))


my_label <- paste0( "P<4.25 x 10-5 ; ",
                    "Positively associated:",table(sig.res$group)[1]," ; ",
                    "Negatively associated:",table(sig.res$group)[2])

my_label <- paste0( "P<0.05 ; ",
                    "Positively associated:",table(sig.res$group)[1]," ; ",
                    "Negatively associated:",table(sig.res$group)[2])


label <-subset(sig.res$id.exposure,sig.res$pval < (0.05/1178))

label <-subset(sig.res$id.exposure,sig.res$id.exposure %in% c("NCAN"))
label <-subset(sig.res$id.exposure,sig.res$id.exposure %in% c("SERPINA1"))
label <-c("NCAN","SERPINA1")


p <- ggscatter(sig.res,
               x = "OR", y = "logP",
               label = "id.exposure",
               label.select = label,
               color = "group", size = 2,
               main = paste0("Replication cohort2(GCST90043859)"), # ***
               xlab = "OR", ylab = "-log10(P.Value)",
               palette = c("#D01910","#00599F","#CCCCCC"),
               ylim = c(-0.1,5),xlim=c(0,3.0))+
  theme_base()+
  geom_hline(yintercept = -log10(0.05/1178), linetype="dashed", color = "#222222") +
  geom_vline(xintercept = 1 , linetype="dashed", color = "#222222") +
  #geom_vline(xintercept = -1, linetype="dashed", color = "#222222") +
  labs(subtitle = my_label)+
  theme(plot.background = element_blank())

p

ggsave("./MR_Protein_vs_GCST90043859_volcano.pdf", p, width = 12, height = 10)
ggsave("./MR_Protein_vs_FinngenR12_volcano.pdf", p, width = 12, height = 10)


p <- ggscatter(sig.res,
               x = "OR", y = "logP",
               label = "id.exposure",
               label.select = label,
               color = "group", size = 2,
               main = paste0("Validation Bile Duct Cancer"), # ***
               xlab = "OR", ylab = "-log10(P.Value)",
               palette = c("#D01910","#00599F","#CCCCCC"),
               ylim = c(-0.1, 7),xlim=c(0.995,1.005))+
  theme_base()+
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "#222222") +
  geom_vline(xintercept = 1 , linetype="dashed", color = "#222222") +
  #geom_vline(xintercept = -1, linetype="dashed", color = "#222222") +
  labs(subtitle = my_label)+
  theme(plot.background = element_blank())

ggsave("./MR_Protein_vs_IEU-b-4915_volcano.pdf", p, width = 12, height = 10)


#########################################################################################################################
mr.or <-fread("mr_result_exposure.PlasmaProtein_outcome.ieu-b-4969.melanoma_cisPQTLs.xls",sep='\t',header=T)
mr.or <- mr.or[mr.or$method!="Weighted median",]
mr.or <- mr.or[mr.or$method!="MR Egger",]

0.05/1753
-log10(0.05/1753)

res <- mr.or %>% filter(!is.na(pval)) %>% 
  mutate( logP = -log10(pval) ) %>%
  mutate( OR = or ) %>%
  mutate( tag = "Melanoma.TMR")%>%
  mutate( Gene = rownames(id.exposure))

sig.res <- res %>% mutate(group=case_when(
  (pval < 0.05/1753 & or > 1) ~ "Positively associated",
  (pval < 0.05/1753 & or < 1) ~ "Negatively associated",
  .default = "Not significant"))

sig.res <- sig.res %>% mutate(group=factor(group,levels = c("Positively associated","Negatively associated","Not significant")))


my_label <- paste0( "P<2.85x10-5 ; ",
                    "Positively associated:",table(sig.res$group)[1]," ; ",
                    "Negatively associated:",table(sig.res$group)[2])

label <-subset(sig.res$id.exposure,sig.res$pval < (0.05/1753))

p <- ggscatter(sig.res,
               x = "OR", y = "logP",
               label = "id.exposure",
               label.select = label,
               color = "group", size = 2,
               main = paste0("Melanoma"), # ***
               xlab = "OR", ylab = "-log10(P.Value)",
               palette = c("#D01910","#CCCCCC"),
               ylim = c(-0.1, 8),xlim=c(0.95,1.05))+
  theme_base()+
  geom_hline(yintercept = -log10(0.05/1753), linetype="dashed", color = "#222222") +
  geom_vline(xintercept = 1 , linetype="dashed", color = "#222222") +
  #geom_vline(xintercept = -1, linetype="dashed", color = "#222222") +
  labs(subtitle = my_label)+
  theme(plot.background = element_blank())

ggsave("./MR_Protein_vs_Skin.melanoma_4969.pdf", p, width = 12, height = 10)

###############################################################################################################CXCL8共定位
######################################################################################################  ASB9共定位 Chromosome X: 15253410-15288589
######################################################################################################  PITPNA共定位 Chromosome 17: 1421012-1466110
######################################################################################################  PHGDH共定位Chromosome 1: 120202421-120286838
######################################################################################################  NEU1共定位 Chromosome 6: 31825436-31830683

###################ASB9
t.chr <-x
t.start <-15253410
t.end <- 15288589

##############################################################################################################################PITPNA
#protein <- fread("/home/data/t020412/Qinzhang_data/MR_analysis/9934_29_PITPNA_PIPNA.txt.gz",sep='\t')

#protein <- fread("/home/data/t020412/Qinzhang_data/MR_analysis/15548_35_PHGDH_SERA.txt.gz",sep='\t')

#protein <- fread("/home/data/t020412/Qinzhang_data/MR_analysis/15426_5_NEU1_NEUR1.txt.gz",sep='\t')

protein <- fread("/home/data/t020412/Qinzhang_data/MR_drug.target_BileDuctCancer/Ferkingstad.15573_110_NCAN_CSPG3.txt",sep='\t')

#########   NCAN 共定位分析
out <- fread ("/home/data/t020412/Qinzhang_data/MR_drug.target_BileDuctCancer/finngen_R12_C3_Bile_outcome.MR.format.xls",sep='\t',header = T)
############### 2025.08 补
out <- fread ("/home/data/t020412/Qinzhang_data/MR_drug.target_BileDuctCancer/BileDuct.cancer_g3859MR.MRinput.txt",sep='\t',header = T)
out <- fread ("/home/data/t020412/Qinzhang_data/MR_drug.target_BileDuctCancer/ieu-b-4915_outcome.MR.format.xls",sep='\t',header = T)


t.chr <-'19' ###finngen
t.chr <-'chr19' #####protein
t.pos <-19329924

protein <- fread("/home/data/t020412/Qinzhang_data/MR_drug.target_BileDuctCancer/Ferkingstad.15573_110_NCAN_CSPG3.txt",sep='\t')
protein$phenotype <- "Proteins"
protein <- as.data.frame(protein)
data <- format_data(protein,type = "exposure",
                    header = TRUE,
                    phenotype_col = "phenotype",
                    snp_col = "rsids",
                    beta_col = "Beta",
                    se_col = "SE",
                    eaf_col = "ImpMAF",
                    effect_allele_col = "effectAllele",
                    other_allele_col = "otherAllele",
                    pval_col = "Pval",
                    chr_col = "Chrom",
                    samplesize_col = "N",
                    id_col = "Proteins",
                    pos_col = "Pos")
pQTL <-data[data$chr.exposure==t.chr,]
pQTL <- pQTL[pQTL$pos.exposure >t.pos-1000000 & pQTL$pos.exposure < t.pos+1000000,]
gwas <- out[out$chr.outcome==t.chr,]
gwas <- gwas[gwas$pos.outcome >t.pos-1000000 & gwas$pos.outcome < t.pos+1000000,]

########## finngen 队列case 2298例，总共381047例
gwas$s <-as.numeric(2298/gwas$samplesize.outcome)

############g3859 队列case 104例
gwas$s <-as.numeric(104/gwas$samplesize.outcome)

############ieu-b-4915 队列case 350例
gwas$s <-as.numeric(350/gwas$samplesize.outcome)

###################ASB9
t.chr <-'23'
t.chr <-'X'
t.start <-15253410
t.end <- 15288589
protein <- fread("/home/data/t020412/Qinzhang_data/MR_analysis/ASB9_19601_15.txt",sep='\t')
protein$phenotype <- "Proteins"
protein <- as.data.frame(protein)
data <- format_data(protein,type = "exposure",
                    header = TRUE,
                    phenotype_col = "phenotype",
                    snp_col = "rsid",
                    beta_col = "Effect",
                    se_col = "StdErr",
                    eaf_col = "MinFreq",
                    effect_allele_col = "Allele1",
                    other_allele_col = "Allele2",
                    pval_col = "Pvalue",
                    chr_col = "chr",
                    samplesize_col = "TotalSampleSize",
                    id_col = "Proteins",
                    pos_col = "pos")


out <- fread ("/home/data/t020412/Qinzhang_data/MR_analysis/MR_Plasma.Protein_vs_Cancer/ieu-b-4959_outcome.MR.format.xls",sep='\t',header = T)

out <- fread ("/home/data/t020412/Qinzhang_data/MR_analysis/MR_Plasma.Protein_vs_Cancer/ieu-b-4969_outcome.MR.format.xls",sep='\t',header = T)

out <- fread ("/home/data/t020412/Qinzhang_data/MR_analysis/MR_Plasma.Protein_vs_Cancer/ukb-b-12915_outcome.MR.format.xls",sep='\t',header = T)



#########SERPINA1

t.chr <-'14' ###finngen protein
t.pos <-94844947

protein <- fread("/home/data/t020412/Qinzhang_data/MR_drug.target_BileDuctCancer/Pietzner.a1_Antitrypsin_3580_25.txt.gz",sep='\t')
protein$phenotype <- "Proteins"
protein <- as.data.frame(protein)
data <- format_data(protein,type = "exposure",
                    header = TRUE,
                    phenotype_col = "phenotype",
                    snp_col = "rsid",
                    beta_col = "Effect",
                    se_col = "StdErr",
                    eaf_col = "MinFreq",
                    effect_allele_col = "Allele1",
                    other_allele_col = "Allele2",
                    pval_col = "Pvalue",
                    chr_col = "chr",
                    samplesize_col = "TotalSampleSize",
                    id_col = "Proteins",
                    pos_col = "pos")
pQTL <-data[data$chr.exposure==t.chr,]
pQTL <- pQTL[pQTL$pos.exposure >t.pos-1000000 & pQTL$pos.exposure < t.pos+1000000,]

gwas <- out[out$chr.outcome==t.chr,]
gwas <- gwas[gwas$pos.outcome >t.pos-1000000 & gwas$pos.outcome < t.pos+1000000,]



####################################################################################共定位

gwas$MAF <- ifelse(gwas$eaf.outcome<0.5,gwas$eaf.outcome,1-gwas$eaf.outcome)
pQTL$MAF <- ifelse(pQTL$eaf.exposure<0.5,pQTL$eaf.exposure,1-pQTL$eaf.exposure)

gwas <- gwas[gwas$MAF != 'NA',]
pQTL <- pQTL[pQTL$MAF !='NA',]

sameSNP <- intersect(pQTL$SNP,gwas$SNP)

pQTL <- pQTL[pQTL$SNP %in% sameSNP,]
gwas <- gwas[gwas$SNP %in% sameSNP,]

##### finngen
gwas$s <-as.numeric(2298/gwas$samplesize.outcome)
####
#result <- coloc.abf(dataset1=list(pvalues=gwas$pval.outcome, snp=gwas$SNP,MAF=gwas$MAF,beta=gwas$beta.outcome, varbeta=gwas$se.outcome^2,type="cc", s=gwas$s[1], N=gwas$samplesize.outcome),dataset2=list(pvalues=pQTL$pval.exposure, snp=pQTL$SNP,MAF=pQTL$MAF,beta=pQTL$beta.exposure, varbeta=pQTL$se.exposure^2, type="quant", N=pQTL$samplesize.exposure), MAF=pQTL$MAF)


result <- coloc.abf(dataset1=list(pvalues=gwas$pval.outcome, snp=gwas$SNP,MAF=gwas$MAF,beta=gwas$beta.outcome, varbeta=gwas$se.outcome^2,type="cc", s=gwas$s[1], N=gwas$samplesize.outcome),dataset2=list(pvalues=pQTL$pval.exposure, snp=pQTL$SNP,MAF=pQTL$MAF,beta=pQTL$beta.exposure, varbeta=pQTL$se.exposure, type="quant", N=pQTL$samplesize.exposure), MAF=pQTL$MAF)

#result <- coloc.abf(dataset1=list(pvalues=gwas$pval.outcome, snp=gwas$SNP,MAF=gwas$MAF,beta=gwas$beta.outcome, varbeta=gwas$se.outcome,type="cc", s=gwas$s[1], N=gwas$samplesize.outcome),dataset2=list(pvalues=pQTL$pval.exposure, snp=pQTL$SNP,MAF=pQTL$MAF,beta=pQTL$beta.exposure, varbeta=pQTL$se.exposure, type="quant", N=pQTL$samplesize.exposure), MAF=pQTL$MAF)


need_result <- result$results %>% dplyr::arrange(desc(SNP.PP.H4))

export(need_result,"NCAN_finngen_coloc_result.xls",format = "\t")
export(need_result,"NCAN_g3859_coloc_result.xls",format = "\t")
export(need_result,"NCAN_ieu-b-4915_coloc_result.xls",format = "\t")

export(need_result,"SERPINA1_finngen_coloc_result.xls",format = "\t")
export(need_result,"SERPINA1_g3859_coloc_result.xls",format = "\t")

gwas_fn <- gwas[,c('SNP','pval.outcome')] %>% dplyr::rename(rsid = SNP, pval = pval.outcome)
pQTL_fn <- pQTL[,c('SNP','pval.exposure')] %>% dplyr::rename(rsid = SNP, pval = pval.exposure)

locuscompare(in_fn1 =pQTL_fn , in_fn2 = gwas_fn, title1 = 'pQTL',title2 = 'GWAS',snp="rs2228603")
locuscompare(in_fn1 =pQTL_fn , in_fn2 = gwas_fn, title1 = 'pQTL',title2 = 'GWAS',snp="rs28929474")

ggsave('NCAN_finngen_Coloc.pdf',width=9,height = 4.5)
ggsave('NCAN_g3859_Coloc.pdf',width=9,height = 4.5)
ggsave('NCAN_ieu-b-4915_Coloc.pdf',width=9,height = 4.5)

ggsave('SERPINA1_finngen_Coloc.pdf',width=9,height = 4.5)
ggsave('SERPINA1_g3859_Coloc.pdf',width=9,height = 4.5)



#######1000G

bfile="/home/data/t020412/R/x86_64-pc-linux-gnu-library/4.3/plinkbinr/EUR"

####################################################################################################################Discovery森林图
mr.or <-fread("mr_result_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer_addOR.xls",sep='\t',header=T)
mr.or <-fread('mr_result_exposure.cis-pQTLs_protein_outcome.ieu4915.bile.duct.cancer_addOR.xls',header = T)
mr.or <-fread('mr_result_exposure.cis-pQTLs_protein_outcome.g3859.bile.duct.cancer_addOR.xls',header = T)

out_multi.finn <-mr.or[mr.or$pval<(0.05/1178),]
out_multi.4915 <-mr.or[mr.or$id.exposure %in% c("NCAN","SERPINA1"),]
out_multi.g3859 <-mr.or[mr.or$id.exposure %in% c("NCAN","SERPINA1"),]

out_multi <-rbind(out_multi.finn,out_multi.4915,out_multi.g3859)
out_multi$cancer.type <-ifelse(out_multi$outcome=="ieu-b-4959","Melanoma","Skin.melanoma")

hz <- paste(round(out_multi$or,3),
            "(",round(out_multi$or_lci95,3),
            "-",round(out_multi$or_uci95,3),")",sep = "")


tabletext <- cbind(c(NA,"Exposure",out_multi$id.exposure),
                   c(NA,"Outcome",out_multi$outcome),
                   c(NA,"Methods",out_multi$method),
                   c(NA,"P value",ifelse(out_multi$pval<0.001,"P < 0.001",round(out_multi$pval,3))),
                   #c(NA,"P value",out_multi$pval),
                   c(NA,"OR(95% CI)",hz))


pforest <- forestplot(labeltext=tabletext, 
                      graph.pos=4,  #为Pvalue箱线图所在的位置
                      col=fpColors(box="#00468B", lines="black", zero = "red"),
                      mean=c(NA,NA,out_multi$or),
                      lower=c(NA,NA,out_multi$or_lci95), #95%置信区间下限
                      upper=c(NA,NA,out_multi$or_uci95), #95%置信区间上限
                      xlab="Odds Ratio(95%CI)",
                      boxsize=0.2,lwd.ci=0.8,   #箱子大小，线的宽度
                      ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
                      zero=1,lwd.zero=1,      #zero线宽 基准线的位置
                      colgap=unit(14,"mm"),    #列间隙
                      xticks = c(0.3,0.5,0.6,0.7,0.8,0.9,1,1.1,1.3,1.5,5,7), #横坐标刻度
                      lwd.xaxis=2,            #X轴线宽
                      lineheight = unit(1,"cm"), #固定行高
                      graphwidth = unit(0.2,"npc"), #图在表中的宽度比例
                      #cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
                      hrzl_lines=list("2" = gpar(lwd=2, col="black"),
                                      "3" = gpar(lwd=2, col="black")),           
                      #第三行顶部加黑线，引号内数字标记行位置
                      #mar=unit(rep(0.1, times = 4), "mm"),#图形页边距
                      clip=c(0,3),
                      txt_gp=fpTxtGp(label=gpar(cex=1.25),
                                     ticks=gpar(cex=1.1),
                                     xlab=gpar(cex = 1.4),
                                     title=gpar(cex = 1.2))
)
pdf("./Plasma_Protein_vs_Discovery_and_validation_Forest.pdf", width=12, height=6)
pforest
dev.off()

####################################################################################################################### Validiation

#mr.or <-fread("mr_result_exposure.PlasmaProtein_outcome.ieu-b-4959.melanoma_addOR_cisPQTLs.xls",sep='\t',header=T)
#mr.or.skin <-fread("mr_result_exposure.PlasmaProtein_outcome.ieu-b-4969.melanoma_cisPQTLs.xls",sep='\t',header=T)
#out_multi.4959 <-mr.or[mr.or$pval<(0.05/1753),]
#out_multi.4969 <-mr.or.skin[mr.or.skin$pval<(0.05/1753),]

out_multi <-fread("Validation_forest_input.xls",sep='\t',header=T)
#out_multi <-rbind(out_multi.4959,out_multi.4969)
#out_multi$cancer.type <-ifelse(out_multi$outcome=="ieu-b-4959","Melanoma","Skin.melanoma")

hz <- paste(round(out_multi$or,3),
            "(",round(out_multi$or_lci95,3),
            "-",round(out_multi$or_uci95,3),")",sep = "")


tabletext <- cbind(c(NA,"Exposure",out_multi$id.exposure),
                   c(NA,"Outcome",out_multi$outcome),
                   c(NA,"Methods",out_multi$method),
                   c(NA,"P value",ifelse(out_multi$pval<0.001,"P < 0.001",round(out_multi$pval,3))),
                   #c(NA,"P value",out_multi$pval),
                   c(NA,"OR(95% CI)",hz))


pforest <- forestplot(labeltext=tabletext, 
                      graph.pos=4,  #为Pvalue箱线图所在的位置
                      col=fpColors(box="#00468B", lines="black", zero = "red"),
                      mean=c(NA,NA,out_multi$or),
                      lower=c(NA,NA,out_multi$or_lci95), #95%置信区间下限
                      upper=c(NA,NA,out_multi$or_uci95), #95%置信区间上限
                      xlab="Odds Ratio(95%CI)",
                      boxsize=0.2,lwd.ci=0.8,   #箱子大小，线的宽度
                      ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
                      zero=1,lwd.zero=1,      #zero线宽 基准线的位置
                      colgap=unit(14,"mm"),    #列间隙
                      xticks = c(0.6,0.7,0.8,0.9,1,1.1,1.3,1.5), #横坐标刻度
                      lwd.xaxis=2,            #X轴线宽
                      lineheight = unit(1,"cm"), #固定行高
                      graphwidth = unit(0.2,"npc"), #图在表中的宽度比例
                      #cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
                      hrzl_lines=list("2" = gpar(lwd=2, col="black"),
                                      "3" = gpar(lwd=2, col="black")),           
                      #第三行顶部加黑线，引号内数字标记行位置
                      #mar=unit(rep(0.1, times = 4), "mm"),#图形页边距
                      clip=c(0,3),
                      txt_gp=fpTxtGp(label=gpar(cex=1.25),
                                     ticks=gpar(cex=1.1),
                                     xlab=gpar(cex = 1.4),
                                     title=gpar(cex = 1.2))
)
pdf("./Plasma_Protein_vs_Validaiton_FinnGen_Forest.pdf", width=12, height=6)
pforest
dev.off()


#################  SMR数据分析
#1.step  gwas数据准备
g <- fread ("/home/data/t020412/Qinzhang_data/MR_analysis/MR_Plasma.Protein_vs_Cancer/ieu-b-4959_outcome.MR.format.xls",sep='\t',header = T)
g <- fread ("/home/data/t020412/Qinzhang_data/MR_analysis/MR_Plasma.Protein_vs_Cancer/ieu-b-4969_outcome.MR.format.xls",sep='\t',header = T)
g <- fread ("/home/data/t020412/Qinzhang_data/MR_drug.target_BileDuctCancer/finngen_R12_C3_Bile_outcome.MR.format.xls",sep='\t',header = T)

f<-select(g,SNP,effect_allele.outcome,other_allele.outcome,eaf.outcome,beta.outcome,se.outcome,pval.outcome,samplesize.outcome)#之前的列名
colnames(f)<-c("SNP","A1","A2","freq","b","se","p","n")

f_clean <- f[!apply(f, 1, function(row) {any(is.na(row) | row == "")}), , drop = FALSE]

f_clean[, "SNP"] <- sub(",.*", "", f[, "SNP"])

export(f_clean,"SMR_GWAS_FinngenR12_input.txt",format = "\t")

################### smr软件位置
#/home/data/t020412/Qinzhang_data/MR_analysis/MR_Plasma.Protein_vs_Cancer/SMR_analysis/smr-1.3.1-linux-x86_64/smr


##############################meta 分析
library(meta)      # 用于meta分析
library(metafor)   # 用于meta回归
data<-read.table("meta_analysis_NCAN.xls",header=T,sep="\t")
data <-subset(data,data$ncase!="104")

beta <- log(data[,"or"])
loguci <- log(data[,"or_uci95"])
loglci <- log(data[,"or_lci95"])

metadata=metagen(TE=data$b,seTE=se,data=data,sm="OR",
                 n.e = ncase,n.c = ncontrol,pval = data$pval,
                 random = TRUE,common = FALSE,
                 studlab=data$Dataset)

pdf('NCAN_meta_analysis.pdf',width=9,height =2.5)
forest(metadata)
dev.off()


data<-read.table("meta_analysis_SERPINA1.xls",header=T,sep="\t")

beta <- log(data[,"or"])
loguci <- log(data[,"or_uci95"])
loglci <- log(data[,"or_lci95"])

metadata=metagen(TE=data$b,seTE=se,data=data,sm="OR",
                 n.e = ncase,n.c = ncontrol,pval = data$pval,
                 random = TRUE,common = FALSE,
                 studlab=data$Dataset)

pdf('SERPINA1_meta_analysis.pdf',width=9,height =2.5)
forest(metadata)
dev.off()

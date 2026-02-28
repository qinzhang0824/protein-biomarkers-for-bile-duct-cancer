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
library(TwoSampleMR)
library(data.table)
library(MendelianRandomization)
library(ieugwasr)
library(plinkbinr)
library(S4Vectors)
library(rio)
library(ggplot2)
library(ggsci)
library(forestplot)
library(VariantAnnotation)
library(limma)
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(ggthemes)
library(grid)

############################################################################################################ read Instruments variables of plasma proteins
protein <- fread("Raw_input_data/Instrumental_variables_of_plasma_proteins_used_in_MR_analysis.txt",sep='\t')
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

##########################################################################Read FinnGen R12 download data and harmonise with Instruments variables of plasma proteins
finngen <-fread('Raw_input_data/finngen_R12_C3_BILIARY_GALLBLADDER_EXALLC',header = T,sep='\t')
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
export(bile.finn,"intermediate_files/finngen_R12_C3_Bile_outcome.MR.format.xls",format = "\t")
mydata <- harmonise_data(exposure_dat=protein.mr,outcome_dat=bile.finn,action= 2)
export(mydata,"intermediate_files/Harmonising_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer.xls",format = "\t")

###################################################################### Read IEU-b-4915 download data (vcf) and harmonise with Instruments variables of plasma proteins
vcf <- VariantAnnotation::readVcf("Raw_input_data/ieu-b-4915.vcf", "hg19")
outcome_dat <- gwasvcf_to_TwoSampleMR(vcf,type = "outcome")

ieu4915<-outcome_dat
mydata <- harmonise_data(exposure_dat=protein.mr,outcome_dat=ieu4915,action= 2)

export(mydata,"intermediate_files/Harmonising_exposure.cis-pQTLs_protein_outcome.ieu4915.bile.duct.cancer.xls",format = "\t")
########################################################################################### Read GCST90043859 download data and harmonise with Instruments variables of plasma proteins
g3859 <-fread('Raw_input_data/GCST90043859_buildGRCh37.tsv',header = T,sep='\t')
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

mydata <- harmonise_data(exposure_dat=protein.mr,outcome_dat=g3859MR,action= 2)
export(mydata,"intermediate_files/Harmonising_exposure.cis-pQTLs_protein_outcome.GCST90043859.bile.duct.cancer.xls",format = "\t")

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

################################################################################################ Figure2 plot

mr.or <-fread('final results/Discovery cohort/mr_result_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer_addOR.xls',header = T)
mr.or <- mr.or[mr.or$method!="Weighted median",]
mr.or <- mr.or[mr.or$method!="MR Egger",]

res <- mr.or %>% filter(!is.na(pval)) %>% 
  mutate( logP = -log10(pval) ) %>%
  mutate( OR = or ) %>%
  mutate( tag = "Validation Bile Duct Cancer")%>%
  mutate( Gene = rownames(id.exposure))

sig.res <- res %>% mutate(group=case_when(
  (pval < 0.05/1193 & or > 1) ~ "Positively associated",
  (pval < 0.05/1193 & or < 1) ~ "Negatively associated",
  .default = "Not significant"))


sig.res <- sig.res %>% mutate(group=factor(group,levels = c("Positively associated","Negatively associated","Not significant")))

label <-subset(sig.res$id.exposure,sig.res$pval < (0.05/1193))
label <-subset(sig.res$id.exposure,sig.res$id.exposure %in% c("NCAN"))
label <-subset(sig.res$id.exposure,sig.res$id.exposure %in% c("SERPINA1"))
label <-c("NCAN","SERPINA1")

my_label <- paste0( "P<4.19 x 10-5 ; ",
                    "Positively associated:",table(sig.res$group)[1]," ; ",
                    "Negatively associated:",table(sig.res$group)[2])
p <- ggscatter(sig.res,
               x = "OR", y = "logP",
               label = "id.exposure",
               label.select = label,
               color = "group", size = 2,
               main = paste0("Discovery cohort(FinngenR12)"), # ***
               xlab = "OR", ylab = "-log10(P.Value)",
               palette = c("#D01910","#00599F","#CCCCCC"),
               ylim = c(-0.1,5),xlim=c(0,3.0))+
  theme_base()+
  geom_hline(yintercept = -log10(0.05/1193), linetype="dashed", color = "#222222") +
  geom_vline(xintercept = 1 , linetype="dashed", color = "#222222") +
  #geom_vline(xintercept = -1, linetype="dashed", color = "#222222") +
  labs(subtitle = my_label)+
  theme(plot.background = element_blank())
p
ggsave("final results/Figure2A_MR_Protein_vs_FinngenR12_volcano.pdf", p, width = 12, height = 10)
##########################################################################
mr.or <-fread('final results/Validation cohort 1/mr_result_exposure.cis-pQTLs_protein_outcome.ieu4915.bile.duct.cancer_addOR.xls',header = T)
mr.or <- mr.or[mr.or$method!="Weighted median",]
mr.or <- mr.or[mr.or$method!="MR Egger",]

res <- mr.or %>% filter(!is.na(pval)) %>% 
  mutate( logP = -log10(pval) ) %>%
  mutate( OR = or ) %>%
  mutate( tag = "Discovery Bile Duct Cancer")%>%
  mutate( Gene = rownames(id.exposure))
sig.res <- res %>% mutate(group=case_when(
  (pval < 0.05 & or > 1) ~ "Positively associated",
  (pval < 0.05 & or < 1) ~ "Negatively associated",
  .default = "Not significant"))

label <-c("NCAN")
my_label <- paste0( "P<0.05 ; ",
                    "Positively associated:",table(sig.res$group)[1]," ; ",
                    "Negatively associated:",table(sig.res$group)[2])
p <- ggscatter(sig.res,
               x = "OR", y = "logP",
               label = "id.exposure",
               label.select = label,
               color = "group", size = 2,
               main = paste0("Replication cohort 1(IEU--b-4915)"), # ***
               xlab = "OR", ylab = "-log10(P.Value)",
               palette = c("#D01910","#00599F","#CCCCCC"),
               ylim = c(-0.1, 7),xlim=c(0.995,1.005))+
  theme_base()+
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "#222222") +
  geom_vline(xintercept = 1 , linetype="dashed", color = "#222222") +
  #geom_vline(xintercept = -1, linetype="dashed", color = "#222222") +
  labs(subtitle = my_label)+
  theme(plot.background = element_blank())

ggsave("final results/Figure2B_MR_Protein_vs_IEU-b-4915_volcano.pdf", p, width = 12, height = 10)
###########################################################
mr.or <-fread('final results/Validation cohort 2/mr_result_exposure.cis-pQTLs_protein_outcome.g3859.bile.duct.cancer_addOR.xls',header = T)

mr.or <- mr.or[mr.or$method!="Weighted median",]
mr.or <- mr.or[mr.or$method!="MR Egger",]

res <- mr.or %>% filter(!is.na(pval)) %>% 
  mutate( logP = -log10(pval) ) %>%
  mutate( OR = or ) %>%
  mutate( tag = "Discovery Bile Duct Cancer")%>%
  mutate( Gene = rownames(id.exposure))
sig.res <- res %>% mutate(group=case_when(
  (pval < 0.05 & or > 1) ~ "Positively associated",
  (pval < 0.05 & or < 1) ~ "Negatively associated",
  .default = "Not significant"))

label <-c("SERPINA1")
my_label <- paste0( "P<0.05 ; ",
                    "Positively associated:",table(sig.res$group)[1]," ; ",
                    "Negatively associated:",table(sig.res$group)[2])
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
  geom_hline(yintercept = -log10(0.05/1193), linetype="dashed", color = "#222222") +
  geom_vline(xintercept = 1 , linetype="dashed", color = "#222222") +
  #geom_vline(xintercept = -1, linetype="dashed", color = "#222222") +
  labs(subtitle = my_label)+
  theme(plot.background = element_blank())
p

ggsave("final results/Figure2C_MR_Protein_vs_GCST90043859_volcano.pdf", p, width = 12, height = 10)

####################################################################################################### Figure3 plot
mr.or <-fread("final results/Discovery cohort/mr_result_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer_addOR.xls",sep='\t',header=T)
mr.or1 <-fread('final results/Validation cohort 1/mr_result_exposure.cis-pQTLs_protein_outcome.ieu4915.bile.duct.cancer_addOR.xls',header = T)
mr.or2 <-fread('final results/Validation cohort 2/mr_result_exposure.cis-pQTLs_protein_outcome.g3859.bile.duct.cancer_addOR.xls',header = T)

out_multi.finn <-mr.or[mr.or$pval<(0.05/1193),]
out_multi.4915 <-mr.or1[mr.or1$id.exposure %in% c("NCAN","SERPINA1"),]
out_multi.g3859 <-mr.or2[mr.or2$id.exposure %in% c("NCAN","SERPINA1"),]

out_multi <-rbind(out_multi.finn,out_multi.4915,out_multi.g3859)

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
                      graph.pos=4,
                      col=fpColors(box="#00468B", lines="black", zero = "red"),
                      mean=c(NA,NA,out_multi$or),
                      lower=c(NA,NA,out_multi$or_lci95), 
                      upper=c(NA,NA,out_multi$or_uci95), 
                      xlab="Odds Ratio(95%CI)",
                      boxsize=0.2,lwd.ci=0.8,  
                      ci.vertices.height = 0.08,ci.vertices=TRUE,
                      zero=1,lwd.zero=1,     
                      colgap=unit(14,"mm"),   
                      xticks = c(0.3,0.5,0.6,0.7,0.8,0.9,1,1.1,1.3,1.5,5,7), 
                      lwd.xaxis=2,         
                      lineheight = unit(1,"cm"), 
                      graphwidth = unit(0.2,"npc"), 
                      #cex=0.9, fn.ci_norm = fpDrawCircleCI, 
                      hrzl_lines=list("2" = gpar(lwd=2, col="black"),
                                      "3" = gpar(lwd=2, col="black")),           
                      clip=c(0,3),
                      txt_gp=fpTxtGp(label=gpar(cex=1.25),
                                     ticks=gpar(cex=1.1),
                                     xlab=gpar(cex = 1.4),
                                     title=gpar(cex = 1.2))
)
pdf("final results/Figure3_Plasma_Protein_vs_Discovery_and_validation_Forest.pdf", width=12, height=6)
pforest
dev.off()
######################################################################################################### Figure4 plot 
#####################################################################  NCAN
t.chr <-'19'
t.pos <-19329924

out <- fread ("/intermediate_files/finngen_R12_C3_Bile_outcome.MR.format_NCAN.xls",sep='\t',header = T)
gwas <- out[out$chr.outcome=="19",]
gwas <- gwas[gwas$pos.outcome >t.pos-1000000 & gwas$pos.outcome < t.pos+1000000,]

protein <- fread("/Raw_input_data/Ferkingstad.15573_110_NCAN_CSPG3.txt",sep='\t')
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
pQTL <-data[data$chr.exposure=="chr19",]
pQTL <- pQTL[pQTL$pos.exposure >t.pos-1000000 & pQTL$pos.exposure < t.pos+1000000,]

########## finngen case 2298
gwas$s <-as.numeric(2298/gwas$samplesize.outcome)
gwas$MAF <- ifelse(gwas$eaf.outcome<0.5,gwas$eaf.outcome,1-gwas$eaf.outcome)
pQTL$MAF <- ifelse(pQTL$eaf.exposure<0.5,pQTL$eaf.exposure,1-pQTL$eaf.exposure)

gwas <- gwas[gwas$MAF != 'NA',]
pQTL <- pQTL[pQTL$MAF !='NA',]

sameSNP <- intersect(pQTL$SNP,gwas$SNP)

pQTL <- pQTL[pQTL$SNP %in% sameSNP,]
gwas <- gwas[gwas$SNP %in% sameSNP,]

result <- coloc.abf(dataset1=list(pvalues=gwas$pval.outcome, snp=gwas$SNP,MAF=gwas$MAF,beta=gwas$beta.outcome, varbeta=gwas$se.outcome^2,type="cc", s=gwas$s[1], N=gwas$samplesize.outcome),dataset2=list(pvalues=pQTL$pval.exposure, snp=pQTL$SNP,MAF=pQTL$MAF,beta=pQTL$beta.exposure, varbeta=pQTL$se.exposure^2, type="quant", N=pQTL$samplesize.exposure), MAF=pQTL$MAF)

gwas_fn <- gwas[,c('SNP','pval.outcome')] %>% dplyr::rename(rsid = SNP, pval = pval.outcome)
pQTL_fn <- pQTL[,c('SNP','pval.exposure')] %>% dplyr::rename(rsid = SNP, pval = pval.exposure)

locuscompare(in_fn1 =pQTL_fn , in_fn2 = gwas_fn, title1 = 'pQTL',title2 = 'GWAS',snp="rs2228603")
ggsave('/final_results/Figure4_NCAN_finngen_Coloc.png',width=9,height = 4.5)

####################################################################### SERPINA1

t.chr <-'14' ###finngen protein
t.pos <-94844947

out <- fread ("/intermediate_files/finngen_R12_C3_Bile_outcome.MR.format_SERPINA1.xls",sep='\t',header = T)
gwas <- out[out$chr.outcome=="14",]
gwas <- gwas[gwas$pos.outcome >t.pos-1000000 & gwas$pos.outcome < t.pos+1000000,]

protein <- fread("/Raw_input_data/Pietzner.a1_Antitrypsin_3580_25.txt",sep='\t')
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
pQTL <-data[data$chr.exposure=="chr14",]
pQTL <- pQTL[pQTL$pos.exposure >t.pos-1000000 & pQTL$pos.exposure < t.pos+1000000,]

gwas$MAF <- ifelse(gwas$eaf.outcome<0.5,gwas$eaf.outcome,1-gwas$eaf.outcome)
pQTL$MAF <- ifelse(pQTL$eaf.exposure<0.5,pQTL$eaf.exposure,1-pQTL$eaf.exposure)

gwas <- gwas[gwas$MAF != 'NA',]
pQTL <- pQTL[pQTL$MAF !='NA',]

sameSNP <- intersect(pQTL$SNP,gwas$SNP)

pQTL <- pQTL[pQTL$SNP %in% sameSNP,]
gwas <- gwas[gwas$SNP %in% sameSNP,]

##### finngen
gwas$s <-as.numeric(2298/gwas$samplesize.outcome)

result <- coloc.abf(dataset1=list(pvalues=gwas$pval.outcome, snp=gwas$SNP,MAF=gwas$MAF,beta=gwas$beta.outcome, varbeta=gwas$se.outcome^2,type="cc", s=gwas$s[1], N=gwas$samplesize.outcome),dataset2=list(pvalues=pQTL$pval.exposure, snp=pQTL$SNP,MAF=pQTL$MAF,beta=pQTL$beta.exposure, varbeta=pQTL$se.exposure^2, type="quant", N=pQTL$samplesize.exposure), MAF=pQTL$MAF)

gwas_fn <- gwas[,c('SNP','pval.outcome')] %>% dplyr::rename(rsid = SNP, pval = pval.outcome)
pQTL_fn <- pQTL[,c('SNP','pval.exposure')] %>% dplyr::rename(rsid = SNP, pval = pval.exposure)

locuscompare(in_fn1 =pQTL_fn , in_fn2 = gwas_fn, title1 = 'pQTL',title2 = 'GWAS',snp="rs28929474")

ggsave('/final_results/Figure4_SERPINA1_finngen_Coloc.png',width=9,height = 4.5)
###############################################################################################################################  SMR input

g <- fread ("intermediate_files/finngen_R12_C3_Bile_outcome.MR.format.xls",sep='\t',header = T)

f<-select(g,SNP,effect_allele.outcome,other_allele.outcome,eaf.outcome,beta.outcome,se.outcome,pval.outcome,samplesize.outcome)
colnames(f)<-c("SNP","A1","A2","freq","b","se","p","n")

f_clean <- f[!apply(f, 1, function(row) {any(is.na(row) | row == "")}), , drop = FALSE]

f_clean[, "SNP"] <- sub(",.*", "", f[, "SNP"])

export(f_clean,"SMR_analysis/SMR_GWAS_FinngenR12_input.txt",format = "\t")



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
mydata.finn <- harmonise_data(exposure_dat=protein.mr,outcome_dat=bile.finn,action= 2)
export(mydata.finn,"intermediate_files/Harmonising_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer.xls",format = "\t")

###################################################################### Read IEU-b-4915 download data (vcf) and harmonise with Instruments variables of plasma proteins
vcf <- VariantAnnotation::readVcf("Raw_input_data/ieu-b-4915.vcf", "hg19")
ieu4915 <- gwasvcf_to_TwoSampleMR(vcf,type = "outcome")
mydata.ieu4915 <- harmonise_data(exposure_dat=protein.mr,outcome_dat=ieu4915,action= 2)

export(mydata.ieu4915,"intermediate_files/Harmonising_exposure.cis-pQTLs_protein_outcome.ieu4915.bile.duct.cancer.xls",format = "\t")
########################################################################################### Read GCST90043859 download data and harmonise with Instruments variables of plasma proteins
g3859 <-fread('Raw_input_data/GCST90043859_buildGRCh37.tsv',header = T,sep='\t')
g3859 <-as.data.frame(g3859)
g3859$phenotype <- "BileDuct.cancer"
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

mydata.g3859 <- harmonise_data(exposure_dat=protein.mr,outcome_dat=g3859MR,action= 2)
export(mydata.g3859,"intermediate_files/Harmonising_exposure.cis-pQTLs_protein_outcome.GCST90043859.bile.duct.cancer.xls",format = "\t")

########################################################################################### MR analysis finngen.R12
#### mydata <- fread('intermediate_files/Harmonising_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer.txt',header = T)

result.finn <- mr(mydata.finn, method_list=c("mr_egger_regression", "mr_ivw","mr_weighted_median","mr_wald_ratio"))
MR.OR.finn<-generate_odds_ratios(result.finn)

export(MR.OR.finn,"Final_results/Discovery_cohort/mr_result_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer_addOR.txt",format = "\t") 

mr_steiger_direction.finn <-directionality_test(mydata.finn)
export(mr_steiger_direction.finn,"Final_results/Discovery_cohort/mr_result_exposure.cis-PlasmaProtein_outcome.finn12_steiger_direction.txt",format = "\t")

###########################################################################################  IEU-b-4915
#### mydata <- fread('intermediate_files/Harmonising_exposure.cis-pQTLs_protein_outcome.ieu4915.bile.duct.cancer.xls',header = T)

result.ieu4915 <- mr(mydata.ieu4915, method_list=c("mr_egger_regression", "mr_ivw","mr_weighted_median","mr_wald_ratio"))
MR.OR.ieu4915<-generate_odds_ratios(result.ieu4915)

export(MR.OR.ieu4915,"Final_results/Validation_cohort_1/mr_result_exposure.cis-pQTLs_protein_outcome.ieu4915.bile.duct.cancer_addOR.xls",format = "\t") 

mr_steiger_direction.ieu4915 <-directionality_test(mydata.ieu4915)
export(mr_steiger_direction.ieu4915,"Final_results/Validation_cohort_1/mr_result_exposure.cis-PlasmaProtein_ieu4915_steiger_direction.xls",format = "\t")

############################################################################################## GCST90043859
### mydata <- fread('intermediate_files/Harmonising_exposure.cis-pQTLs_protein_outcome.GCST90043859.bile.duct.cancer.xls',header = T)

result.g3859 <- mr(mydata.g3859, method_list=c("mr_egger_regression", "mr_ivw","mr_weighted_median","mr_wald_ratio"))
MR.OR.g3859<-generate_odds_ratios(result.g3859)

export(MR.OR.g3859,"Final_results/Validation_cohort_2/mr_result_exposure.cis-pQTLs_protein_outcome.g3859.bile.duct.cancer_addOR.xls",format = "\t") 

mr_steiger_direction.g3859 <-directionality_test(mydata.g3859)
export(mr_steiger_direction.g3859,"Final_results/Validation_cohort_2/mr_result_exposure.cis-PlasmaProtein_g3859_steiger_direction.xls",format = "\t")

################################################################################################ Figure2 plot

mr.or.finn <- MR.OR.finn[MR.OR.finn$method!="Weighted median",]
mr.or.finn <- mr.or.finn[mr.or.finn$method!="MR Egger",]

res.finn <- mr.or.finn %>% filter(!is.na(pval)) %>% 
  mutate( logP = -log10(pval) ) %>%
  mutate( OR = or ) %>%
  mutate( tag = "Discovery Bile Duct Cancer")%>%
  mutate( Gene = rownames(id.exposure))

sig.res.finn <- res.finn %>% mutate(group=case_when(
  (pval < 0.05/1193 & or > 1) ~ "Positively associated",
  (pval < 0.05/1193 & or < 1) ~ "Negatively associated",
  .default = "Not significant"))


sig.res.finn <- sig.res.finn %>% mutate(group=factor(group,levels = c("Positively associated","Negatively associated","Not significant")))

label <-subset(sig.res.finn$id.exposure,sig.res.finn$pval < (0.05/1193))
label <-subset(sig.res.finn$id.exposure,sig.res.finn$id.exposure %in% c("NCAN"))
label <-subset(sig.res.finn$id.exposure,sig.res.finn$id.exposure %in% c("SERPINA1"))
label <-c("NCAN","SERPINA1")

my_label <- paste0( "P<4.19 x 10-5 ; ",
                    "Positively associated:",table(sig.res$group)[1]," ; ",
                    "Negatively associated:",table(sig.res$group)[2])

p <- ggscatter(sig.res.finn,
               x = "OR", y = "logP",
               label = "id.exposure",
               label.select = label,
               color = "group", size = 2,
               main = paste0("Discovery Bile Duct Cancer"), # ***
               xlab = "OR", ylab = "-log10(P.Value)",
               palette = c("#00599F","#CCCCCC"),
               ylim = c(-0.1,15),xlim=c(0,3.0))+
  theme_base()+
  geom_hline(yintercept = -log10(0.05/1193), linetype="dashed", color = "#222222") +
  geom_vline(xintercept = 1 , linetype="dashed", color = "#222222") +
  #geom_vline(xintercept = -1, linetype="dashed", color = "#222222") +
  labs(subtitle = my_label)+
  theme(plot.background = element_blank())
p
ggsave("Figure2A_MR_Protein_vs_FinngenR12_volcano.pdf", p, width = 12, height = 10)
##########################################################################
mr.or.ieu4915 <- MR.OR.ieu4915[MR.OR.ieu4915$method!="Weighted median",]
mr.or.ieu4915 <- mr.or.ieu4915[mr.or.ieu4915$method!="MR Egger",]

res.ieu4915 <- mr.or.ieu4915 %>% filter(!is.na(pval)) %>% 
  mutate( logP = -log10(pval) ) %>%
  mutate( OR = or ) %>%
  mutate( tag = "Replication cohort 1(IEU-b-4915)")%>%
  mutate( Gene = rownames(id.exposure))
sig.res.ieu4915 <- res.ieu4915 %>% mutate(group=case_when(
  (pval < 0.05 & or > 1) ~ "Positively associated",
  (pval < 0.05 & or < 1) ~ "Negatively associated",
  .default = "Not significant"))

sig.res.ieu4915 <- sig.res.ieu4915 %>% mutate(group=factor(group,levels = c("Positively associated","Negatively associated","Not significant")))

label <-c("NCAN")
my_label <- paste0( "P<0.05 ; ",
                    "Positively associated:",table(sig.res.ieu4915$group)[1]," ; ",
                    "Negatively associated:",table(sig.res.ieu4915$group)[2])
p <- ggscatter(sig.res.ieu4915,
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

ggsave("Figure2B_MR_Protein_vs_IEU-b-4915_volcano.pdf", p, width = 12, height = 10)

###########################################################

mr.or.g3859 <- MR.OR.g3859[MR.OR.g3859$method!="Weighted median",]
mr.or.g3859 <- mr.or.g3859[mr.or.g3859$method!="MR Egger",]

res.g3859 <- mr.or.g3859 %>% filter(!is.na(pval)) %>% 
  mutate( logP = -log10(pval) ) %>%
  mutate( OR = or ) %>%
  mutate( tag = "Discovery Bile Duct Cancer")%>%
  mutate( Gene = rownames(id.exposure))
sig.res.g3859 <- res.g3859 %>% mutate(group=case_when(
  (pval < 0.05 & or > 1) ~ "Positively associated",
  (pval < 0.05 & or < 1) ~ "Negatively associated",
  .default = "Not significant"))

sig.res.g3859 <- sig.res.g3859 %>% mutate(group=factor(group,levels = c("Positively associated","Negatively associated","Not significant")))

label <-c("SERPINA1")
my_label <- paste0( "P<0.05 ; ",
                    "Positively associated:",table(sig.res.g3859$group)[1]," ; ",
                    "Negatively associated:",table(sig.res.g3859$group)[2])
p <- ggscatter(sig.res.g3859,
               x = "OR", y = "logP",
               label = "id.exposure",
               label.select = label,
               color = "group", size = 2,
               main = paste0("Replication cohort2(GCST90043859)"), # ***
               xlab = "OR", ylab = "-log10(P.Value)",
               palette = c("#D01910","#00599F","#CCCCCC"),
               ylim = c(-0.1,5),xlim=c(0,3.0))+
  theme_base()+
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color = "#222222") +
  geom_vline(xintercept = 1 , linetype="dashed", color = "#222222") +
  #geom_vline(xintercept = -1, linetype="dashed", color = "#222222") +
  labs(subtitle = my_label)+
  theme(plot.background = element_blank())
p

ggsave("Figure2C_MR_Protein_vs_GCST90043859_volcano.pdf", p, width = 12, height = 10)

####################################################################################################### Figure3 plot
##mr.or <-fread("Final_results/Discovery_cohort/mr_result_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer_addOR.xls",sep='\t',header=T)
##mr.or1 <-fread('Final_results/Validation_cohort_1/mr_result_exposure.cis-pQTLs_protein_outcome.ieu4915.bile.duct.cancer_addOR.xls',header = T)
##mr.or2 <-fread('Final_results/Validation_cohort_2/mr_result_exposure.cis-pQTLs_protein_outcome.g3859.bile.duct.cancer_addOR.xls',header = T)

out_multi.finn <-MR.OR.finn[MR.OR.finn$pval<(0.05/1193),]
out_multi.4915 <-MR.OR.ieu4915[MR.OR.ieu4915$id.exposure %in% c("NCAN","SERPINA1"),]
out_multi.g3859 <-MR.OR.g3859[MR.OR.g3859$id.exposure %in% c("NCAN","SERPINA1"),]

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
pdf("Final_results/Figures/Figure3_Plasma_Protein_vs_Discovery_and_validation_Forest.pdf", width=12, height=6)
pforest
dev.off()
######################################################################################################### Figure4 plot 
#####################################################################  NCAN
t.chr <-'19'
t.pos <-19329924

gwas.NCAN <- bile.finn[bile.finn$chr.outcome=="19",]
gwas.NCAN <- gwas.NCAN[gwas.NCAN$pos.outcome >t.pos-1000000 & gwas.NCAN$pos.outcome < t.pos+1000000,]

protein.NCAN <- fread("/Raw_input_data/Ferkingstad.15573_110_NCAN_CSPG3.txt",sep='\t')
protein.NCAN$phenotype <- "Proteins"
protein.NCAN <- as.data.frame(protein.NCAN)
data.NCAN <- format_data(protein.NCAN,type = "exposure",
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
pQTL.NCAN <-data.NCAN[data.NCAN$chr.exposure=="chr19",]
pQTL.NCAN <- pQTL.NCAN[pQTL.NCAN$pos.exposure >t.pos-1000000 & pQTL.NCAN$pos.exposure < t.pos+1000000,]

########## finngen case 2298
gwas.NCAN$s <-as.numeric(2298/gwas.NCAN$samplesize.outcome)
gwas.NCAN$MAF <- ifelse(gwas.NCAN$eaf.outcome<0.5,gwas.NCAN$eaf.outcome,1-gwas.NCAN$eaf.outcome)
pQTL.NCAN$MAF <- ifelse(pQTL.NCAN$eaf.exposure<0.5,pQTL.NCAN$eaf.exposure,1-pQTL.NCAN$eaf.exposure)

gwas.NCAN <- gwas.NCAN[gwas.NCAN$MAF != 'NA',]
pQTL.NCAN <- pQTL.NCAN[pQTL.NCAN$MAF !='NA',]

sameSNP.NCAN <- intersect(pQTL.NCAN$SNP,gwas.NCAN$SNP)

pQTL.NCAN <- pQTL.NCAN[pQTL.NCAN$SNP %in% sameSNP.NCAN,]
gwas.NCAN <- gwas.NCAN[gwas.NCAN$SNP %in% sameSNP.NCAN,]

result.NCAN <- coloc.abf(dataset1=list(pvalues=gwas.NCAN$pval.outcome, snp=gwas.NCAN$SNP,MAF=gwas.NCAN$MAF,beta=gwas.NCAN$beta.outcome, varbeta=gwas.NCAN$se.outcome^2,type="cc", s=gwas.NCAN$s[1], N=gwas.NCAN$samplesize.outcome),dataset2=list(pvalues=pQTL.NCAN$pval.exposure, snp=pQTL.NCAN$SNP,MAF=pQTL.NCAN$MAF,beta=pQTL.NCAN$beta.exposure, varbeta=pQTL.NCAN$se.exposure^2, type="quant", N=pQTL.NCAN$samplesize.exposure), MAF=pQTL.NCAN$MAF)

gwas_fn.NCAN <- gwas.NCAN[,c('SNP','pval.outcome')] %>% dplyr::rename(rsid = SNP, pval = pval.outcome)
pQTL_fn.NCAN <- pQTL.NCAN[,c('SNP','pval.exposure')] %>% dplyr::rename(rsid = SNP, pval = pval.exposure)

locuscompare(in_fn1 =pQTL_fn.NCAN , in_fn2 = gwas_fn.NCAN, title1 = 'pQTL',title2 = 'GWAS',snp="rs2228603")
ggsave('/Final_results/Figures/Figure4_NCAN_finngen_Coloc.png',width=9,height = 4.5)

####################################################################### SERPINA1
t.chr <-'14' ###finngen protein
t.pos <-94844947

gwas.SERPINA1 <- bile.finn[bile.finn$chr.outcome=="14",]
gwas.SERPINA1 <- gwas.SERPINA1[gwas.SERPINA1$pos.outcome >t.pos-1000000 & gwas.SERPINA1$pos.outcome < t.pos+1000000,]

protein.SERPINA1 <- fread("Pietzner.a1_Antitrypsin_3580_25.txt.gz",sep='\t')
protein.SERPINA1$phenotype <- "Proteins"
protein.SERPINA1 <- as.data.frame(protein.SERPINA1)
data.SERPINA1 <- format_data(protein.SERPINA1,type = "exposure",
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
pQTL.SERPINA1 <-data.SERPINA1[data.SERPINA1$chr.exposure=="14",]
pQTL.SERPINA1 <- pQTL.SERPINA1[pQTL.SERPINA1$pos.exposure >t.pos-1000000 & pQTL.SERPINA1$pos.exposure < t.pos+1000000,]

gwas.SERPINA1$MAF <- ifelse(gwas.SERPINA1$eaf.outcome<0.5,gwas.SERPINA1$eaf.outcome,1-gwas.SERPINA1$eaf.outcome)
pQTL.SERPINA1$MAF <- ifelse(pQTL.SERPINA1$eaf.exposure<0.5,pQTL.SERPINA1$eaf.exposure,1-pQTL.SERPINA1$eaf.exposure)

gwas.SERPINA1 <- gwas.SERPINA1[gwas.SERPINA1$MAF != 'NA',]
pQTL.SERPINA1 <- pQTL.SERPINA1[pQTL.SERPINA1$MAF !='NA',]

sameSNP.SERPINA1 <- intersect(pQTL.SERPINA1$SNP,gwas.SERPINA1$SNP)

pQTL.SERPINA1 <- pQTL.SERPINA1[pQTL.SERPINA1$SNP %in% sameSNP.SERPINA1,]
gwas.SERPINA1 <- gwas.SERPINA1[gwas.SERPINA1$SNP %in% sameSNP.SERPINA1,]

##### finngen
gwas.SERPINA1$s <-as.numeric(2298/gwas.SERPINA1$samplesize.outcome)

result.SERPINA1 <- coloc.abf(dataset1=list(pvalues=gwas.SERPINA1$pval.outcome, snp=gwas.SERPINA1$SNP,MAF=gwas.SERPINA1$MAF,beta=gwas.SERPINA1$beta.outcome, varbeta=gwas.SERPINA1$se.outcome^2,type="cc", s=gwas.SERPINA1$s[1], N=gwas.SERPINA1$samplesize.outcome),dataset2=list(pvalues=pQTL.SERPINA1$pval.exposure, snp=pQTL.SERPINA1$SNP,MAF=pQTL.SERPINA1$MAF,beta=pQTL.SERPINA1$beta.exposure, varbeta=pQTL.SERPINA1$se.exposure^2, type="quant", N=pQTL.SERPINA1$samplesize.exposure), MAF=pQTL.SERPINA1$MAF)

gwas_fn.SERPINA1 <- gwas.SERPINA1[,c('SNP','pval.outcome')] %>% dplyr::rename(rsid = SNP, pval = pval.outcome)
pQTL_fn.SERPINA1 <- pQTL.SERPINA1[,c('SNP','pval.exposure')] %>% dplyr::rename(rsid = SNP, pval = pval.exposure)

locuscompare(in_fn1 =pQTL_fn.SERPINA1 , in_fn2 = gwas_fn.SERPINA1, title1 = 'pQTL',title2 = 'GWAS',snp="rs28929474")

ggsave('Figure4_SERPINA1_finngen_Coloc.png',width=9,height = 4.5)

###############################################################################################################################  SMR input

f<-select(bile.finn,SNP,effect_allele.outcome,other_allele.outcome,eaf.outcome,beta.outcome,se.outcome,pval.outcome,samplesize.outcome)
colnames(f)<-c("SNP","A1","A2","freq","b","se","p","n")

f_clean <- f[!apply(f, 1, function(row) {any(is.na(row) | row == "")}), , drop = FALSE]

f_clean[, "SNP"] <- sub(",.*", "", f[, "SNP"])

export(f_clean,"SMR_analysis/SMR_GWAS_FinngenR12_input.txt",format = "\t")



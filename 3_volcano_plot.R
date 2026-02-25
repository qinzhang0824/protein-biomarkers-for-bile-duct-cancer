##################################################################################################################### 
library(limma)
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(ggthemes)
library(TwoSampleMR)
library(data.table)

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
######################################################################################################
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
##############################################################################################################
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






#1.step  gwas
library(dplyr)
library(data.table)
library(readxl)
library(plyr)


g <- fread ("/finngen_R12_C3_Bile_outcome.MR.format.xls",sep='\t',header = T)

f<-select(g,SNP,effect_allele.outcome,other_allele.outcome,eaf.outcome,beta.outcome,se.outcome,pval.outcome,samplesize.outcome)
colnames(f)<-c("SNP","A1","A2","freq","b","se","p","n")

f_clean <- f[!apply(f, 1, function(row) {any(is.na(row) | row == "")}), , drop = FALSE]

f_clean[, "SNP"] <- sub(",.*", "", f[, "SNP"])

export(f_clean,"SMR_GWAS_FinngenR12_input.txt",format = "\t")

#1.step  gwas
g <- fread ("/home/data/t020412/Qinzhang_data/MR_drug.target_BileDuctCancer/finngen_R12_C3_Bile_outcome.MR.format.xls",sep='\t',header = T)

f<-select(g,SNP,effect_allele.outcome,other_allele.outcome,eaf.outcome,beta.outcome,se.outcome,pval.outcome,samplesize.outcome)
colnames(f)<-c("SNP","A1","A2","freq","b","se","p","n")

f_clean <- f[!apply(f, 1, function(row) {any(is.na(row) | row == "")}), , drop = FALSE]

f_clean[, "SNP"] <- sub(",.*", "", f[, "SNP"])

export(f_clean,"SMR_GWAS_FinngenR12_input.txt",format = "\t")

########################################
SMR_analysis/smr-1.3.1-linux-x86_64/smr --bfile R/x86_64-pc-linux-gnu-library/4.3/plinkbinr/EUR --gwas-summary SMR_GWAS_FinngenR12_input_NoD.txt --beqtl-summary MR_analysis/GTEx_V8_cis_eqtl_summary_lite/eQTL_besd_lite/Whole_Blood.lite --out FinngenR12.cis.eQTL.GTEX.v8

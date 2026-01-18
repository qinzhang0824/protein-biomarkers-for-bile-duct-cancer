library(TwoSampleMR)
library(readr)
library(dplyr)
library(data.table)
library(readxl)
library(plyr)
library(MRPRESSO)
library("locuscomparer")
library(coloc)

out <- fread ("finngen_R12_C3_Bile_outcome.MR.format.xls",sep='\t',header = T)

| chr.outcome | pos.outcome | other_allele.outcome | effect_allele.outcome | SNP | pval.outcome | beta.outcome | se.outcome | eaf.outcome | samplesize.outcome | outcome |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 13668 | G | A | rs2691328 | 0.314029 | -0.365242 | 0.362774 | 0.00541903 | 381047 | Bile.duct.cancer |
| 1 | 14506 | G | A | rs1240557819 | 0.0440494 | 0.765322 | 0.380072 | 0.00589771 | 381047 | Bile.duct.cancer |
| 1 | 14521 | C | T | rs1378626194 | 0.159142 | 1.09296 | 0.77627 | 0.00143347 | 381047 | Bile.duct.cancer |


t.chr <-'chr19'
t.pos <-19329924
##############NCAN
protein <- fread("Ferkingstad.15573_110_NCAN_CSPG3.txt",sep='\t')
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

########## finngen case 2298
gwas$s <-as.numeric(2298/gwas$samplesize.outcome)


#########SERPINA1

t.chr <-'14' ###finngen protein
t.pos <-94844947

protein <- fread("Pietzner.a1_Antitrypsin_3580_25.txt.gz",sep='\t')
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


gwas$MAF <- ifelse(gwas$eaf.outcome<0.5,gwas$eaf.outcome,1-gwas$eaf.outcome)
pQTL$MAF <- ifelse(pQTL$eaf.exposure<0.5,pQTL$eaf.exposure,1-pQTL$eaf.exposure)

gwas <- gwas[gwas$MAF != 'NA',]
pQTL <- pQTL[pQTL$MAF !='NA',]

sameSNP <- intersect(pQTL$SNP,gwas$SNP)

pQTL <- pQTL[pQTL$SNP %in% sameSNP,]
gwas <- gwas[gwas$SNP %in% sameSNP,]

##### finngen
gwas$s <-as.numeric(2298/gwas$samplesize.outcome)

result <- coloc.abf(dataset1=list(pvalues=gwas$pval.outcome, snp=gwas$SNP,MAF=gwas$MAF,beta=gwas$beta.outcome, varbeta=gwas$se.outcome^2,type="cc", s=gwas$s[1], N=gwas$samplesize.outcome),dataset2=list(pvalues=pQTL$pval.exposure, snp=pQTL$SNP,MAF=pQTL$MAF,beta=pQTL$beta.exposure, varbeta=pQTL$se.exposure, type="quant", N=pQTL$samplesize.exposure), MAF=pQTL$MAF)

gwas_fn <- gwas[,c('SNP','pval.outcome')] %>% dplyr::rename(rsid = SNP, pval = pval.outcome)
pQTL_fn <- pQTL[,c('SNP','pval.exposure')] %>% dplyr::rename(rsid = SNP, pval = pval.exposure)

locuscompare(in_fn1 =pQTL_fn , in_fn2 = gwas_fn, title1 = 'pQTL',title2 = 'GWAS',snp="rs2228603")
locuscompare(in_fn1 =pQTL_fn , in_fn2 = gwas_fn, title1 = 'pQTL',title2 = 'GWAS',snp="rs28929474")

ggsave('NCAN_finngen_Coloc.pdf',width=9,height = 4.5)
ggsave('NCAN_g3859_Coloc.pdf',width=9,height = 4.5)
ggsave('NCAN_ieu-b-4915_Coloc.pdf',width=9,height = 4.5)

ggsave('SERPINA1_finngen_Coloc.pdf',width=9,height = 4.5)
ggsave('SERPINA1_g3859_Coloc.pdf',width=9,height = 4.5)

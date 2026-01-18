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
library("ggpubr")
library(forestplot)

############################################################################################################ read Instruments variables of plasma proteins
protein <- fread("Raw_input_data/Instrumental variables of plasma proteins used in MR analysis.xls",sep='\t')
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
p.sig <- 0.05/1178

##########################################################################Read FinnGen R12 download data and harmonise with Instruments variables of plasma proteins
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
mydata <- harmonise_data(exposure_dat=protein.mr,outcome_dat=bile.finn,action= 2)
export(mydata,"intermediate files/Harmonising_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer.xls",format = "\t")

###################################################################### Read IEU-b-4915 download data (vcf) and harmonise with Instruments variables of plasma proteins
vcf <- VariantAnnotation::readVcf("ieu-b-4915.vcf", "hg19")
outcome_dat <- gwasvcf_to_TwoSampleMR(vcf,type = "outcome")

ieu4915<-outcome_dat
mydata <- harmonise_data(exposure_dat=protein.mr,outcome_dat=ieu4915,action= 2)

export(mydata,"intermediate files/Harmonising_exposure.cis-pQTLs_protein_outcome.ieu4915.bile.duct.cancer.xls",format = "\t")
########################################################################################### Read GCST90043859 download data and harmonise with Instruments variables of plasma proteins
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

mydata <- harmonise_data(exposure_dat=protein.mr,outcome_dat=g3859MR,action= 2)
export(mydata,"Harmonising_exposure.cis-pQTLs_protein_outcome.GCST90043859.bile.duct.cancer.xls",format = "\t")

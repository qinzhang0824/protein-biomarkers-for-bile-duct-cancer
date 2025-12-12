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

.libPaths(c('~/R/x86_64-pc-linux-gnu-library/4.3',
            '/refdir/Rlib_4.3',
            '/usr/local/lib/R/library'))


#########Instruments
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
p.sig <- 0.05/1178

############################################################################################## FinnGen R12 
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

############################################################################################# MR 分析cis-pQTL VS finngen R12

mydata <- harmonise_data(exposure_dat=protein.mr,outcome_dat=bile.finn,action= 2)
export(mydata,"Harmonising_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer.xls",format = "\t")
mydata <- fread('Harmonising_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer.xls',header = T)

###########################################################################################  IEU-b-4915
vcf <- VariantAnnotation::readVcf("ieu-b-4915.vcf", "hg19")
outcome_dat <- gwasvcf_to_TwoSampleMR(vcf,type = "outcome")

ieu4915<-outcome_dat
mydata <- harmonise_data(exposure_dat=protein.mr,outcome_dat=ieu4915,action= 2)

export(mydata,"Harmonising_exposure.cis-pQTLs_protein_outcome.ieu4915.bile.duct.cancer.xls",format = "\t")
mydata <- fread('Harmonising_exposure.cis-pQTLs_protein_outcome.ieu4915.bile.duct.cancer.xls',header = T)

############################################################################################ GCST90043859
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


export(g3859MR,"BileDuct.cancer_g3859MR.MRinput.txt",format = "\t") 

mydata <- harmonise_data(exposure_dat=protein.mr,outcome_dat=g3859MR,action= 2)
export(mydata,"Harmonising_exposure.cis-pQTLs_protein_outcome.GCST90043859.bile.duct.cancer.xls",format = "\t")

sessionInfo()

R version 4.4.3 (2025-02-28)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 20.04.6 LTS

Matrix products: default
BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so;  LAPACK version 3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Asia/Shanghai
tzcode source: system (glibc)

attached base packages:
 [1] splines   grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] rio_1.2.4                     CellChat_1.6.1                bigmemory_4.6.4               igraph_2.2.1                 
 [5] scDblFinder_1.20.2            SingleCellExperiment_1.28.1   DoubletFinder_2.0.6           magrittr_2.0.4               
 [9] reticulate_1.44.1             BiocParallel_1.40.2           SCP_0.5.6                     scCustomize_3.0.1            
[13] patchwork_1.3.2               celldex_1.16.0                SingleR_2.8.0                 SummarizedExperiment_1.36.0  
[17] GenomicRanges_1.58.0          GenomeInfoDb_1.42.3           IRanges_2.40.1                MatrixGenerics_1.18.1        
[21] matrixStats_1.5.0             scRNAtoolVis_0.1.0            jjAnno_0.0.3                  ggunchull_1.0.1              
[25] harmony_1.2.4                 Rcpp_1.1.0                    labeling_0.4.3                lubridate_1.9.4              
[29] forcats_1.0.1                 stringr_1.5.2                 purrr_1.1.0                   tidyr_1.3.1                  
[33] tibble_3.3.0                  tidyverse_2.0.0               cowplot_1.2.0                 clustree_0.5.1               
[37] ggraph_2.2.2                  monocle_2.34.0                DDRTree_0.1.5                 irlba_2.3.5.1                
[41] VGAM_1.1-14                   Biobase_2.66.0                Matrix_1.6-5                  Seurat_5.3.1                 
[45] SeuratObject_5.2.0            sp_2.2-0                      forestplot_3.1.7              abind_1.4-8                  
[49] checkmate_2.3.3               ggpubr_0.6.2                  ggsci_4.0.0                   ggplot2_4.0.0                
[53] S4Vectors_0.44.0              BiocGenerics_0.52.0           plinkbinr_0.0.0.9000          MendelianRandomization_0.10.0
[57] gassocplot_0.0.2              gwasglue_0.0.0.9000           ieugwasr_1.0.1                gwasvcf_0.1.2                
[61] coloc_5.2.3                   locuscomparer_1.0.0           MRPRESSO_1.0                  plyr_1.8.9                   
[65] readxl_1.4.5                  data.table_1.17.8             dplyr_1.1.4                   readr_2.1.6                  
[69] TwoSampleMR_0.5.8 

# protein-biomarkers-for-bile-duct-cancer

# Data Folder

## 1. Raw input data

(1) Instrumental variables of plasma proteins used in MR analysis

The reference for the instrumental variables of plasma proteins are (DOI: 10.1186/s13073-023-01229-9)
The raw instrumental variables of plasma proteins is "Raw_input_data/Plasma.Protein_Instruments.Raw.xls"

The cis-pQTLs were selected based on the following criteria:
(a). SNPs linked to any protein with a significance level of P < 5×10-8 were included;
(b). SNPs and proteins located within the Major Histocompatibility Complex (MHC) region on chromosome 6 (25.5–34.0 Mb) were excluded because of the intricate linkage disequilibrium (LD) patterns in this area;
(c). subsequently, LD clumping was performed to identify independent protein quantitative trait loci (pQTLs) for each protein, ensuring that variants had an r² value less than 0.001;
(d). F-statistics genetic instruments were higher than 10

Ultimately, the instrumental variables file used for analysis is "Raw_input_data/Instrumental variables of plasma proteins used in MR analysis.xls"

(2) Three GWAS summary data of Bile tract cancer (BTC) patients

Three GWAS summary data of BTC patients and controls of European ancestry from three public cohorts.

| GWAS ID | Cohort | Trait | Consortium | Sample size | ncase | ncontrol | Population | Number of SNPs |
| :--- | :--- | :--- | :--- | ---: | ---: | ---: | :--- | ---: |
| C3_BILIARY_GALLBLADDER_EXALLC | Discovery cohort | Malignant neoplasm of intrahepatic ducts, biliary tract and gallbladder | FinnGen R12 | 381,047 | 2,298 | 378,749 | European | 21,304,529 |
| ieu-b-4915 | Validation cohort 1 | Liver & bile duct cancer | Burrows | 372,366 | 350 | 372,016 | European | 7,687,713 |
| GCST90043859 | Validation cohort 2 | Intrahepatic bile duct carcinoma | Jiang L | 456,348 | 104 | 456,244 | European | 11,842,647 |


Notes:

(a). The raw C3_BILIARY_GALLBLADDER_EXALLC data link: https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_C3_BILIARY_GALLBLADDER_EXALLC.gz

(b). The raw ieu-b-4915 VCF data link: https://opengwas.io/datasets/ieu-b-4915#

(c). The raw GCST90043859 data link: https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90043001-GCST90044000/GCST90043859/GCST90043859_buildGRCh37.tsv.gz

## 2. intermediate files

After obtaining the GWAS data and plasma protein instrumental variables, we extracted the SNP information of the exposure and outcome, harmonized the data, and ultimately ensured that the effect alleles of SNPs in the exposure and outcome data were consistent. 

The three obtained datasets are as follows：

(a). "intermediate files/Harmonising_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer.xls": C3_BILIARY_GALLBLADDER_EXALLC harmonized with  "Raw_input_data/Instrumental variables of plasma proteins used in MR analysis.xls"

(c). "intermediate files/Harmonising_exposure.cis-pQTLs_protein_outcome.ieu4915.bile.duct.cancer.xls": ieu-b-4915 harmonized with  "Raw_input_data/Instrumental variables of plasma proteins used in MR analysis.xls"

(c). "intermediate files/Harmonising_exposure.cis-pQTLs_protein_outcome.GCST90043859.bile.duct.cancer.xls": GCST90043859 harmonized with  "Raw_input_data/Instrumental variables of plasma proteins used in MR analysis.xls"

## 3.final results

### final results/Discovery cohort

"mr_result_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer_addOR.xls" : The results of the MR analysis for the Discovery cohort(C3_BILIARY_GALLBLADDER_EXALLC)

"mr_result_exposure.cis-PlasmaProtein_outcome.finn12_steiger_direction.xls" : The results of the Steiger test for the Discovery cohort(C3_BILIARY_GALLBLADDER_EXALLC)

### final results/Validation cohort 1

"mr_result_exposure.cis-pQTLs_protein_outcome.ieu4915.bile.duct.cancer_addOR.xls" : The results of the MR analysis for the Validation cohort 1(ieu-b-4915)

"mr_result_exposure.cis-PlasmaProtein_ieu4915_steiger_direction.xls" : The results of the Steiger test for the Validation cohort 1(ieu-b-4915)

### final results/Validation cohort 2

"mr_result_exposure.cis-pQTLs_protein_outcome.g3859.bile.duct.cancer_addOR.xls" : The results of the MR analysis for the Validation cohort 1(GCST90043859)

"mr_result_exposure.cis-PlasmaProtein_g3859_steiger_direction.xls" : The results of the Steiger test for the Validation cohort 1(GCST90043859)

# Prerequisites

## R dependencies

### R version 4.4.3 

Users running Platform: x86_64-pc-linux-gnu and Running under: Ubuntu 20.04.6 LTS, to install the latest version of TwoSampleMR_0.5.8

```r
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")
```

**other attached packages and versions**

| 1 | 2 | 3 | 4 |
|---|---|---|---|
| TwoSampleMR_0.5.8 | plyr_1.8.9 | gassocplot_0.0.2 | ggplot2_4.0.0 |
| readr_2.1.6 | MRPRESSO_1.0 | MendelianRandomization_0.10.0 | ggsci_4.0.0 |
| dplyr_1.1.4 | locuscomparer_1.0.0 | ieugwasr_1.0.1 | ggpubr_0.6.2 |
| data.table_1.17.8 | coloc_5.2.3 | S4Vectors_0.44.0 | forestplot_3.1.7 |
| readxl_1.4.5 | gwasglue_0.0.0.9000 | rio_1.2.4 | ggthemes_5.2.0 |

## Summary‑data‑based MR (SMR) and HEIDI test

The summary data-based Mendelian randomization (SMR) and HEIDI analyses were conducted using SMR software version 1.3.1. For the SMR test

**SMR software was upload to SMA_analysis folder**

The SMR software path: SMR_analysis/smr-1.3.1-linux-x86_64.zip

### 1000G data

[1000G_EUR.bed](https://zenodo.org/records/13738569/files/1kg_hg38_filtered.bed.gz?download=1)

[1000G_EUR.bim](https://zenodo.org/records/13738569/files/1kg_hg38_filtered.bim.gz?download=1)

[1000G_EUR.fam](https://zenodo.org/records/13738569/files/1kg_hg38_filtered.fam.gz?download=1)

**Download the above three files to the "SMR_analysis" folder to prepare for the subsequent analysis.**

### eQTL data link

GTEx_V8 data link: https://yanglab.westlake.edu.cn/data/SMR/GTEx_V8_cis_eqtl_summary_lite.tar

eQTL data link: https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/SMR_formatted/cis-eQTL-SMR_20191212.tar.gz

**Download the above three files to the "SMR_analysis" folder to prepare for the subsequent analysis.**


# Getting Started

## Step 1: Data preparation

summary-level GWAS and cis-pQTL data that consist of at least 11 required columns for GWAS and 11 required columns for cis-pQTL. Other columns will be ignored during processing.

Notice: the names of required columns MUST be consistent with the following headers, while the order can be inconsistent. The required data format is presented as follows:

GWAS summary stats:

| chr.outcome | pos.outcome | other_allele.outcome | effect_allele.outcome | SNP | pval.outcome | beta.outcome | se.outcome | eaf.outcome | samplesize.outcome | outcome |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 13668 | G | A | rs2691328 | 0.314029 | -0.365242 | 0.362774 | 0.00541903 | 381047 | Bile.duct.cancer |
| 1 | 14506 | G | A | rs1240557819 | 0.0440494 | 0.765322 | 0.380072 | 0.00589771 | 381047 | Bile.duct.cancer |
| 1 | 14521 | C | T | rs1378626194 | 0.159142 | 1.09296 | 0.77627 | 0.00143347 | 381047 | Bile.duct.cancer |

cis-pQTL summary stats:

| Protein | N | SNP | Chr | Pos | Effect allele | Other allele | Eaf | Beta | Se | P |
| :--- | ---: | :--- | :--- | :--- | :--- | :--- | ---: | ---: | ---: | ---: |
| A1BG | 10,708 | rs893184 | 19 | 58864479 | T | C | 0.05 | -1.07 | 0.03 | 1.00E-200 |
| A2M | 10,708 | rs226384 | 12 | 9263647 | T | C | 0.35 | -0.14 | 0.01 | 7.14E-27 |
| A2ML1 | 10,708 | rs1558526 | 12 | 9009820 | A | G | 0.26 | -0.3 | 0.02 | 9.50E-90 |

Data preparation as follows:

```r
Rscript 1.data_process.R
```

## Step 2: MR analysis

MR analysis as follows:

```r
Rscript 2.MR_analysis.R
```

## Step 3: Figure2 plot

MR analysis as follows:

```r
Rscript 3. volcano_plot.R
```

## Step 4: Figure3 plot

Figure3 forest plot as follows:

```r
Rscript 4.Forest_plot.R
```

## Step 5: Figure4 plot

SERPINA1 and NCAN plasma proteins data link:

[Pietzner.a1_Antitrypsin_3580_25.txt.gz](https://download.decode.is/form/folder/proteomics)

[Ferkingstad.15573_110_NCAN_CSPG3.txt](https://omicscience.org/data_usage_agreement.php)

Figure4 forest plot as follows:

```r
Rscript 5. colocalization_analysis.R
```

## Step 6: Figure5 plot

Figure 5 is derived from the online analysis results of STRING, and the website address is as follows: 

[STRING](https://cn.string-db.org/)

NCAN STRING analysis result: https://cn.string-db.org/cgi/network?taskId=bzERFKo1o0Px&sessionId=b2H7LbHajq0R                        

SERPINA1 STRING analysis result: https://cn.string-db.org/cgi/network?taskId=bBPX6TOgzExM&sessionId=bGBRvdvhhDib

Figure 5C NCAN enrichment analysis related parameters:

<img width="1027" height="727" alt="image" src="https://github.com/user-attachments/assets/5464f620-b996-4eb3-878d-c055b1c205a0" />

Figure 5D SERPINA1 enrichment analysis related parameters:

<img width="1026" height="727" alt="image" src="https://github.com/user-attachments/assets/94b4dbdc-ccdb-45db-8cbf-fda88974be59" />

## Step 7: SMR and HEIDI tests

eQTLGen consortium analysis:

```r
SMR_analysis/smr-1.3.1-linux-x86_64/smr --bfile SMR_analysis/1000G_EUR --gwas-summary SMR_analysis/SMR_GWAS_FinngenR12_input.txt --beqtl-summary SMR_analysis/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense --out ./final results/FinngenR12.cis.eQTL
```

eQTL.GTEX.v8 analysis:

```r
SMR_analysis/smr-1.3.1-linux-x86_64/smr --bfile SMR_analysis/1000G_EUR --gwas-summary SMR_analysis/SMR_GWAS_FinngenR12_input.txt --beqtl-summary SMR_analysis/GTEx_V8_cis_eqtl_summary_lite/eQTL_besd_lite/Whole_Blood.lite --out ./final results/FinngenR12.cis.eQTL.GTEX.v8
```












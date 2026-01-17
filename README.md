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

<img width="2208" height="295" alt="image" src="https://github.com/user-attachments/assets/a80b892b-00ba-44d7-b142-86b8c1a9234e" />


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

(a) final results/Discovery cohort

"mr_result_exposure.cis-pQTLs_protein_outcome.finngen.R12.bile.duct.cancer_addOR.xls" : The results of the MR analysis for the Discovery cohort(C3_BILIARY_GALLBLADDER_EXALLC)

"mr_result_exposure.cis-PlasmaProtein_outcome.finn12_steiger_direction.xls" : The results of the Steiger test for the Discovery cohort(C3_BILIARY_GALLBLADDER_EXALLC)

1. The "Data" folder contains the original or preliminary input data used in the MR analysis process.
2. The "Analysis_data" folder contains some intermediate files or result files generated during the MR analysis process.
3. For SMR and HEIDI analysis, refer to https://yanglab.westlake.edu.cn/software/smr/#SMR&HEIDIanalysis.

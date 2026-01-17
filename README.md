# protein-biomarkers-for-bile-duct-cancer

1. Raw input data
   
(1)Instrumental variables of plasma proteins used in MR analysis

The reference for the instrumental variables of plasma proteins are DOI: 10.1186/s13073-023-01229-9

The cis-pQTLs were selected based on the following criteria:
(a). SNPs linked to any protein with a significance level of P < 5×10-8 were included;

(b). SNPs and proteins located within the Major Histocompatibility Complex (MHC) region on chromosome 6 (25.5–34.0 Mb) were excluded because of the intricate linkage disequilibrium (LD) patterns in this area;

(c). subsequently, LD clumping was performed to identify independent protein quantitative trait loci (pQTLs) for each protein, ensuring that variants had an r² value less than 0.001;

(d). F-statistics genetic instruments were higher than 10

Ultimately, the instrumental variables file used for analysis is "Instrumental variables of plasma proteins used in MR analysis.xls"

(2)Three GWAS summary data of BTC patients

Three GWAS summary data of BTC patients and controls of European ancestry from three public cohorts.
<img width="1936" height="221" alt="image" src="https://github.com/user-attachments/assets/640a03c3-b481-44fe-be62-4791a0e022f3" />




1. The "Data" folder contains the original or preliminary input data used in the MR analysis process.
2. The "Analysis_data" folder contains some intermediate files or result files generated during the MR analysis process.
3. For SMR and HEIDI analysis, refer to https://yanglab.westlake.edu.cn/software/smr/#SMR&HEIDIanalysis.

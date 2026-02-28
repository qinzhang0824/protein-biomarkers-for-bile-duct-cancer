####### Note: This "download_inputs.sh" script performs retrieval, decompression, renaming, and placement where possible, and otherwise fails with a clear message stating exactly what to download, from where, and under what filename to place it

################################################### Download three GWAS summary data of Bile tract cancer (BTC) patients
mkdir /Raw_input_data

####### 1.Download Discovery cohort(C3_BILIARY_GALLBLADDER_EXALLC/FinnGen R12)
wget -O /Raw_input_data/finngen_R12_C3_BILIARY_GALLBLADDER_EXALLC.gz "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_C3_BILIARY_GALLBLADDER_EXALLC.gz"
gunzip /Raw_input_data/finngen_R12_C3_BILIARY_GALLBLADDER_EXALLC.gz

####### 2.Download Validation cohort 1(ieu-b-4915)
### The raw ieu-b-4915 VCF data link: https://opengwas.io/datasets/ieu-b-4915#
### Notes: The Download link(s) will be valid for 2 hours. Click on a link to copy to clipboard.Replace the connection with the following https link

wget -O /Raw_input_data/ieu-b-4915.vcf.gz "https://ieup4.objectstorage.uk-london-1.oci.customer-oci.com/p/Pa8OpFsBhWNAhx1BYH5Nmei_y4PeuO4H5pZ9imYOdcr5ynGtvE7JRActIghlac1a/n/ieup4/b/igd/o/ieu-b-4915/ieu-b-4915.vcf.gz"
gunzip /Raw_input_data/ieu-b-4915.vcf.gz

####### 3.Download Validation cohort 2(GCST90043859)
wget -O /Raw_input_data/GCST90043859_buildGRCh37.tsv.gz "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90043001-GCST90044000/GCST90043859/GCST90043859_buildGRCh37.tsv.gz" 
gunzip /Raw_input_data/GCST90043859_buildGRCh37.tsv.gz

##################################################### Download Summary‑data‑based MR (SMR) and HEIDI test analysis-related files
mkdir /SMR_analysis

####### 4.Download SMR software
wget -O /SMR_analysis/smr-1.4.0-linux-x86_64.zip "https://yanglab.westlake.edu.cn/software/smr/download/smr-1.4.0-linux-x86_64.zip"
cd /SMR_analysis
unzip smr-1.4.0-linux-x86_64.zip

####### 5.Download PLINK GRCh37 1000 Genomes config files
wget -O /SMR_analysis/1000G_EUR.bed "https://zenodo.org/records/6614170/files/1000G_EUR.bed?download=1"
wget -O /SMR_analysis/1000G_EUR.bim "https://zenodo.org/records/6614170/files/1000G_EUR.bim?download=1"
wget -O /SMR_analysis/1000G_EUR.fam "https://zenodo.org/records/6614170/files/1000G_EUR.fam?download=1"

####### Download eQTL data files
#### 6.Download eQTL GTEx_V8
wget -O /SMR_analysis/GTEx_V8_cis_eqtl_summary_lite.tar "https://yanglab.westlake.edu.cn/data/SMR/GTEx_V8_cis_eqtl_summary_lite.tar"
tar -xf /SMR_analysis/GTEx_V8_cis_eqtl_summary_lite.tar
cd /SMR_analysis/GTEx_V8_cis_eqtl_summary_lite
unzip eQTL_besd_lite.zip

#### 7.Download eQTL "cis-eQTL-SMR_20191212.tar.gz"
wget -O /SMR_analysis/cis-eQTL-SMR_20191212.tar.gz "https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/SMR_formatted/cis-eQTL-SMR_20191212.tar.gz"
tar -xzf /SMR_analysis/cis-eQTL-SMR_20191212.tar.gz
gunzip /SMR_analysis/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.epi.gz
gunzip /SMR_analysis/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.besd.gz
gunzip /SMR_analysis/cis-eQTLs-full_eQTLGen_AF_incl_nr_formatted_20191212.new.txt_besd-dense.esi.gz

##################################################### Download Bayesian colocalization analysis-related files.
#### 8.Download target protein GWAS raw data files
## 8.1 The method for downloading the original GWAS data of the protein NCAN is as follows:
# data link : https://www.synapse.org/Synapse:syn51824537
# step1: Create one SYNAPSE account
# step2: Find a file named "CSPG3_15573_110.txt.gz" (The name of the gene NCAN is also known as CSPG3)
# step3: Click “Download” and save the file to the folder named “Raw_input_data”
# step4: Change the file name
mv /Raw_input_data/CSPG3_15573_110.txt.gz /Raw_input_data/Ferkingstad.15573_110_NCAN_CSPG3.txt.gz
gunzip /Raw_input_data/Ferkingstad.15573_110_NCAN_CSPG3.txt.gz











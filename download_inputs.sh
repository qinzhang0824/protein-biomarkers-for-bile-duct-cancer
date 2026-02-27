####### Note: This "download_inputs.sh" script performs retrieval, decompression, renaming, and placement where possible, and otherwise fails with a clear message stating exactly what to download, from where, and under what filename to place it

## Download three GWAS summary data of Bile tract cancer (BTC) patients
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

## Download SMR relatated files
mkdir /SMR_analysis

####### Download SMR software

####### Note: This "download_inputs.sh" script performs retrieval, decompression, renaming, and placement where possible, and otherwise fails with a clear message stating exactly what to download, from where, and under what filename to place it

####### Download three GWAS summary data of Bile tract cancer (BTC) patients
### 1.Download Discovery cohort(C3_BILIARY_GALLBLADDER_EXALLC/FinnGen R12)
wget https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_C3_BILIARY_GALLBLADDER_EXALLC.gz /Raw_input_data/finngen_R12_C3_BILIARY_GALLBLADDER_EXALLC.gz
gunzip /Raw_input_data/finngen_R12_C3_BILIARY_GALLBLADDER_EXALLC.gz

### 2.Download Validation cohort 1(ieu-b-4915)
wget 

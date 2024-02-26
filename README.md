# MICB 475: Data Science Research in Microbiology - Team 2 - Lab Notebook
This document contains summaries of what was done throughout the project

## February 26, 2024
<img width="970" alt="image" src="https://github.com/Xpado-star/MICB-475-Team-2/assets/158794595/585a1732-3b90-4d0a-81b5-dffb9b7bf307">


## February 12, 2024
David - Edited the R file, moved into R folder

John - fixed R file path script

**Ali** - Changed 'treatment' column from numeric to factor; saved metadata as .txt file; created new RProject for phyloseq object; created script for phyloseq object 'pd_phyloseq' and saved as .RData file  

## February 11, 2024
John - Created README file. Imported and demultiplexed parkinsons sequences. Denoising, clusterng, taxonomic training, filtering, and export conducted.

Ali - Removed NAs from PD drug columns. Added "treatment" column numbered 1-7 for the following treatments:  
control = 1   
pd_untreated = 2  
pd_entac = 3  
pd_prami = 4  
pd_rasag = 5  
pd_amant = 6  
pd_combo = 7  

### Visualization of sequences on view.qiime.org
John - All high quality reads throughout. Setting a standard of 30 for median quality score, no trimming is required. 251 is the trimming parameter and the length of all of the sequences.

## February 12, 2024
John - Uploaded new pd_metadata_treatment.tsv to server. Compiled table, taxa bar plots, and no mitochondria no chloroplast table as viewable .qzv files with new pd_metadata_treatment.tsv file.

John - Using table-no-mitochondria-no-chloroplast_treatment.qzv, we found a maximum sampling depth of 5421 in order to retain at least 6 samples in all treatment categories. Alpha-rarefaction curve displayed that 5421 falls under the plateau.

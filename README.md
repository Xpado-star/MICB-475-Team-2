# MICB 475: Data Science Research in Microbiology - Team 2 - Lab Notebook  
Team Members - Ali, Alicia, Cayden, David, John

# Project Summary
Parkinson’s Disease (PD) is a neurodegenerative disorder associated with dopaminergic neuron loss, leading to dopamine dysregulation. Dopaminergic therapeutics are often administered to restore dopamine levels and have been associated with changes to the gut microbiota. Through a secondary data analysis of a cross-sectional cohort of PD patients, we aimed to explore the use of four dopaminergic drugs (entacapone, pramipexole, rasagiline, amantadine) and the associated changes in the gut microbiome. Although the use of dopaminergic therapeutics was not associated with compositional alterations to the microbial diversity of PD patients, we observed changes to specific taxa. Amantadine and pramipexole therapeutics were both associated with a core microbiome that contains *Faecalibacterium* – a genus contained in the core microbiome of healthy individuals but absent in untreated PD patients. Furthermore, entacapone and amantadine use was associated with taxa that are indicative of a healthy gut microbiome, including *Lachnospiraceae* and *Colidextribacter*. We also identified three genera that were differentially abundant with dopaminergic drug use. Dopaminergic therapeutic use was generally associated with increased *Bifidobacterium*, decreased *Prevotella*, and increased *Akkermansia*. While increased *Bifidobacterium* is associated with a healthier gut microbiome and *Akkermansia* is associated with gut dysbiosis, the effects of *Prevotella* remain unclear. Our findings suggest that dopaminergic therapeutics are associated with alterations in the gut microbiome of PD patients that provide an overall benefit to the host. Future studies could incorporate higher resolution analysis at the species level, and explore causational effects of dopaminergic drugs in a prospective study.

# Links to Project Aims
[**Aim 1A**](https://github.com/Xpado-star/MICB-475-Team-2/tree/main/qiime_files) - QIIME2 processing  
[**Aim 1B**](https://github.com/Xpado-star/MICB-475-Team-2/tree/main/Metadata) - Metadata wrangling  
[**Aim 2**](https://github.com/Xpado-star/MICB-475-Team-2/tree/main/alpha_beta_diversity) - Alpha and beta diversity analysis  
[**Aim 3**](https://github.com/Xpado-star/MICB-475-Team-2/tree/main/core_microbiome) - Core microbiome analysis  
[**Aim 4**](https://github.com/Xpado-star/MICB-475-Team-2/tree/main/Indicator%20Taxa) - Indicator species analysis  
[**Aim 5**](https://github.com/Xpado-star/MICB-475-Team-2/tree/main/DESEq2) - Differential abundance analysis

# Lab notebook
Includes brief summaries of what was completed throughout the project duration.
## April 15, 2024
Alicia - updated core microbiome

## April 11, 2024
Cayden - updated volcano plot titles

## April 6, 2024
Ali - updated the diversity plots

## April 3, 2024
Alicia - updated the color codes for core microbiome figures

## March 31, 2024
Cayden - updated DESeq plots

## March 30, 2024
Ali - Linear regression analysis of identified genera from DESeq

## March 25, 2024
David - Added meeting notes of all previous meetings up to this point
Ali - Added additional files
Alicia - Added abundances for core microbiome along with a 3 way venn diagram

## March 24, 2024
David - Minor changes to ISA

## March 23, 2024
Cayden - minor changes to plots

## March 22, 2024
Cayden - Updated the genus bar plots

## March 21, 2024
David - Changed ISA output, adjusted and fixed some code
Cayden - Added genus sums for DESeq analysis

## March 20, 2024
Cayden - updated DESeq barplot titles

## March 18, 2024
Cayden - Reformatted barplots
Alicia - reconducted core microbiome analysis at 1% detection

## March 17, 2024
Cayden: 
- uploaded optional code for volcano plots including labels
- further DESeq analysis results
- including excel sheets (csvs) for comparison across all the genuses

John - changed formatting

## March 16, 2024
Cayden - conducted DESeq analysis and uploaded code for DESeq, generated figures and plots
John - fixing project files and file paths for DESeq

## March 12, 2024
Ali - conducted pairwise beta diversity analysis and comparisons

## March 11, 2024
David - Changed code for ISA, added excel sheets to generate tables for manuscript

## March 10, 2024
Alicia - adding to core microbiome code and providing venn diagram images as potential figures

## March 9, 2024
Ali - addded the alpha and beta diversity metrics, along with code

## March 7, 2024
David - Updated script for indicator species analysis, added R script for ISA, edited some code

## March 4, 2024
Ali - Added Meeting agenda, core microbiome analysis

## March 3, 2024
Alicia - Added TSV and code for rough heatmap

## February 26, 2024
<img width="970" alt="image" src="https://github.com/Xpado-star/MICB-475-Team-2/assets/158794595/585a1732-3b90-4d0a-81b5-dffb9b7bf307">

## February 12, 2024
David - Edited the R file, moved into R folder

John - fixed R file path script

Ali - Changed 'treatment' column from numeric to factor; saved metadata as .txt file; created new RProject for phyloseq object; created script for phyloseq object 'pd_phyloseq' and saved as .RData file  

John - Uploaded new pd_metadata_treatment.tsv to server. Compiled table, taxa bar plots, and no mitochondria no chloroplast table as viewable .qzv files with new pd_metadata_treatment.tsv file.

John - Using table-no-mitochondria-no-chloroplast_treatment.qzv, we found a maximum sampling depth of 5421 in order to retain at least 6 samples in all treatment categories. Alpha-rarefaction curve displayed that 5421 falls under the plateau.

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



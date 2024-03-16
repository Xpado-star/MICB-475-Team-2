######### creating the phyloseq object ###########
setwd("~/Documents/GitHub/MICB-475")

library(tidyverse)
library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)

#### Load data ####

#load metadata RData from Git and rename to 'meta'
load("../R_files/Metadata/pd_metadata_treatment.Rdata")
assign("meta", pd_data_treatment)


otufp <- "qiime_files/parkinsons_export/table_export/feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "qiime_files/parkinsons_export/taxonomy_export/taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "qiime_files/parkinsons_export/rooted_tree_export/tree.nwk"
phylotree <- read.tree(phylotreefp)


#### Format OTU table ####
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)
# [1] "phyloseq"

#### Format metadata ####
samp_df <- as.data.frame(meta[,-1])
rownames(samp_df)<- meta$`#SampleID`
SAMP <- sample_data(samp_df)
class(SAMP)
# [1] "phyloseq"


#### format taxonomy table #####
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$`Feature ID`
TAX <- tax_table(tax_mat)
class(TAX)
# [1] "phyloseq"


#### create phyloseq object #####
pd_phyloseq <- phyloseq(OTU, SAMP, TAX, phylotree)


#### filter and process ####
pd_filter <- subset_taxa(pd_phyloseq,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
pd_filter <- filter_taxa(pd_filter, function(x) sum(x)>5, prune = TRUE)

pd_phyloseq <- pd_filter

#### save as RData file ####
save(pd_phyloseq, file = "phyloseq/pd_phyloseq.RData")


#### rarefy phyloseq object ####
pd_rare <- rarefy_even_depth(pd_phyloseq, rngseed = 1, sample.size = 5421)

#### save as RData file ####
save(pd_rare, file = "phyloseq/pd_rare.RData")

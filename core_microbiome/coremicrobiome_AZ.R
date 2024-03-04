library(tidyverse)
library(phyloseq)
library(microbiome)
library(RColorBrewer)
library(reshape)

setwd("~/Documents/GitHub/MICB-475")

##### load phyloseq object #####
load("phyloseq/pd_phyloseq.RData")

##### convert to relative abundance #####
pd_RA <- transform_sample_counts(pd_phyloseq, fun=function(x) x/sum(x))


#### filter by drug treatment group ####
pd_control <- subset_samples(pd_RA, treatment =="1")
pd_untreat <- subset_samples(pd_RA, treatment =="2")
pd_entac <- subset_samples(pd_RA, treatment =="3")
pd_prami <- subset_samples(pd_RA, treatment =="4")
pd_rasag <- subset_samples(pd_RA, treatment =="5")
pd_amant <- subset_samples(pd_RA, treatment =="6")
pd_combo <- subset_samples(pd_RA, treatment =="7")


#### make the heatmaps #####
prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(0.01), log10(.2), length = 9), 3)

treatment_group = pd_amant
genus_names <- as.character(tax_table(treatment_group)[, "Genus"])
species_names <- as.character(tax_table(treatment_group)[, "Species"])
genus_species <- paste(genus_names, species_names)

p <- plot_core(treatment_group,
               plot.type = "heatmap", 
               colours = rev(brewer.pal(5, "Spectral")),
               prevalences = prevalences, 
               detections = detections, 
               min.prevalence = prevalence(treatment_group, sort = TRUE)[100]) +
  labs(x = "Detection Threshold\n(Relative Abundance (%))") +
  scale_y_discrete(labels = genus_species) +
  
  #Adjusts axis text size and legend bar height
  theme(axis.text.y= element_text(size=8, face="italic"),
        axis.text.x.bottom=element_text(size=8),
        axis.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))

print(p)

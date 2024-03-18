library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library("sf")
install.packages("gplots")
library(gplots)
install.packages("VennDiagram")
library(VennDiagram)

setwd("~/Desktop/MICB-475-Team-2")

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

#### set detection threshold at 2% and prevalence at 50%
control_ASVs <- core_members(pd_control, detection=0.01, prevalence = 0.5)
untreat_ASVs <- core_members(pd_untreat, detection=0.01, prevalence = 0.5)
entac_ASVs <- core_members(pd_entac, detection=0.01, prevalence = 0.5)
prami_ASVs <- core_members(pd_prami, detection=0.01, prevalence = 0.5)
rasag_ASVs <- core_members(pd_rasag, detection=0.01, prevalence = 0.5)
amant_ASVs <- core_members(pd_amant, detection=0.01, prevalence = 0.5)
combo_ASVs <- core_members(pd_combo, detection=0.01, prevalence = 0.5)

#### compare 2 groups
compare_1_3 <- list(Control = control_ASVs, Entacapone = entac_ASVs)
compare_1_4 <- list(Control = control_ASVs, Pramipexole = prami_ASVs)
compare_1_5 <- list(Control = control_ASVs, Rasagiline = rasag_ASVs)
compare_1_6 <- list(Control = control_ASVs, Amantadine = amant_ASVs)
compare_1_7 <- list(Control = control_ASVs, Combo = combo_ASVs)
compare_2_3 <- list(Untreated = untreat_ASVs, Entacapone = entac_ASVs)
compare_2_4 <- list(Untreated = untreat_ASVs, Pramipexole = prami_ASVs)
compare_2_5 <- list(Untreated = untreat_ASVs, Rasagiline = rasag_ASVs)
compare_2_6 <- list(Untreated = untreat_ASVs, Amantadine = amant_ASVs)
compare_2_7 <- list(Untreated = untreat_ASVs, Combo = combo_ASVs)

#### generate venn diagrams
venn_1_3 <- ggVennDiagram(x = compare_1_3)
venn_1_4 <- ggVennDiagram(x = compare_1_4)
venn_1_5 <- ggVennDiagram(x = compare_1_5)
venn_1_6 <- ggVennDiagram(x = compare_1_6)
venn_1_7 <- ggVennDiagram(x = compare_1_7)
venn_2_3 <- ggVennDiagram(x = compare_2_3)
venn_2_4 <- ggVennDiagram(x = compare_2_4)
venn_2_5 <- ggVennDiagram(x = compare_2_5)
venn_2_6 <- ggVennDiagram(x = compare_2_6)
venn_2_7 <- ggVennDiagram(x = compare_2_7)

ggsave("untreated_vs_combo.png", venn_2_7)

#### back up code to pull species
treatment_group = untreat_ASVs
genus_names <- as.character(tax_table(treatment_group)[, "Genus"])
species_names <- as.character(tax_table(treatment_group)[, "Species"])
genus_species <- paste(genus_names, species_names)

#### code used for the species excel form
#need to load object tax_mat from phyloseq code first, location: phyloseq/phyloseq.R
species <- tax_mat[untreat_ASVs, ]

#### list ASV ids
ItemsList <- venn(compare_1_3, show.plot = FALSE)
attributes(ItemsList)$intersections # list common ASV id
venn_list <- VennDiagram::get.venn.partitions(compare_2_3)


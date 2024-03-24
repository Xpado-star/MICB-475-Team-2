library(tidyverse)
library(phyloseq)
library(ggplot2)


##### load phyloseq object #####
load("../phyloseq/pd_phyloseq.RData")

####  Bifidobacterium graph ####
taxa = c("g__Bifidobacterium")
taxa_level = "Genus"
#Subset the phylsoeq to the taxa of interest
sub_beta_RA <- subset_taxa(pd_phyloseq, Genus == taxa)
#Collapse the ASVs of similar taxa rank (EX. Make all species of the same genus into one row)
sub_beta_order_RA <-tax_glom(sub_beta_RA, taxrank = taxa_level, NArm = FALSE)


#Pull out the ASV matrix (relative abundance transformed) 
ASV_df  =  data.frame(t(otu_table(sub_beta_order_RA)))
colnames(ASV_df) = "abs"
#make a new column of the patient ID
ASV_df$ID = rownames(ASV_df)

#pull out the metadata 
meta = data.frame(sample_data(sub_beta_order_RA))
#make a new column with patient ID
meta$ID = rownames(meta)

#Merge the ASV_df and metadata by the patient ID
df_joined = inner_join(meta, ASV_df , by = "ID")

#Making a vector with how many individuals per group for normalization later
individuals_per_treatment = c(nrow(filter(df_joined,df_joined$treatment == 1)),
                              nrow(filter(df_joined,df_joined$treatment == 2)),
                              nrow(filter(df_joined,df_joined$treatment == 3)),
                              nrow(filter(df_joined,df_joined$treatment == 4)),
                              nrow(filter(df_joined,df_joined$treatment == 5)),
                              nrow(filter(df_joined,df_joined$treatment == 6)))
#change the treatments to meaningful information
df_joined$treatment = dplyr::recode(df_joined$treatment,
                                    "1" = "Non-PD",
                                    "2" = "PD-untreated",
                                    "3" = "Entacapone",
                                    "4" = "Pramipexole",
                                    "5" = "Rasagiline",
                                    "6" = "Amantadine")
df_joined <- df_joined %>%
  filter(treatment != 7) %>%
  rename(Treatment = treatment)


#Calculate the sum of relative abundance for each treatment
df_sumed = df_joined %>%
  group_by(Treatment) %>%
  summarise(sum_abs = sum(abs))

#calculating the relative abundance of the group 
df_sumed$rel_abs = df_sumed$sum_abs/individuals_per_treatment


#generating the plot
Bifidobacterium_bar <- ggplot(df_sumed, aes(Treatment, rel_abs, fill = Treatment)) +
  geom_col(color = "black") + # Plotting bars
  scale_fill_manual(values = c("Non-PD" = "sienna1", "PD-untreated" = "indianred2", 
                               "Entacapone" = "#9ACD32", "Pramipexole" = "#569866", 
                               "Rasagiline" = "#104E8B", "Amantadine" = "steelblue2")) + # Setting custom colors
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + # Adjusting expansion of y-axis
  theme_classic() + # Changing plot theme to a white background
  theme(
    axis.title = element_text(face = "bold", size = rel(2)), # Changing axis label font
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black", face = "bold", size = rel(2)), # Changing y axis text font
    legend.text = element_text(size = rel(1.2)), # Changing legend text font
    legend.title = element_text(face = "bold", size = rel(1.5)) # Changing legend title font
  ) +
  labs(x = NULL, y = "Abundance (normalized)", title = "") +
  ggtitle("Bifidobacterium\n Abundance") +
  theme(plot.title = element_text(face = "bold", size = rel(2.5), hjust = 0.5)) 



#### Prevotella graph ####


taxa = c("g__Prevotella")
taxa_level = "Genus"
#Subset the phylsoeq to the taxa of interest
sub_beta_RA <- subset_taxa(pd_phyloseq, Genus == taxa)
#Collapse the ASVs of similar taxa rank (EX. Make all species of the same genus into one row)
sub_beta_order_RA <-tax_glom(sub_beta_RA, taxrank = taxa_level, NArm = FALSE)


#Pull out the ASV matrix (relative abundance transformed) 
ASV_df  =  data.frame(t(otu_table(sub_beta_order_RA)))
colnames(ASV_df) = "abs"
#make a new column of the patient ID
ASV_df$ID = rownames(ASV_df)

#pull out the metadata 
meta = data.frame(sample_data(sub_beta_order_RA))
#make a new column with patient ID
meta$ID = rownames(meta)

#Merge the ASV_df and metadata by the patient ID
df_joined = inner_join(meta, ASV_df , by = "ID")

#Making a vector with how many individuals per group for normalization later
individuals_per_treatment = c(nrow(filter(df_joined,df_joined$treatment == 1)),
                              nrow(filter(df_joined,df_joined$treatment == 2)),
                              nrow(filter(df_joined,df_joined$treatment == 3)),
                              nrow(filter(df_joined,df_joined$treatment == 4)),
                              nrow(filter(df_joined,df_joined$treatment == 5)),
                              nrow(filter(df_joined,df_joined$treatment == 6)))
#change the treatments to meaningful information
df_joined$treatment = dplyr::recode(df_joined$treatment,
                                    "1" = "Non-PD",
                                    "2" = "PD-untreated",
                                    "3" = "Entacapone",
                                    "4" = "Pramipexole",
                                    "5" = "Rasagiline",
                                    "6" = "Amantadine")
df_joined <- df_joined <- df_joined %>%
  filter(treatment != 7) %>%
  rename(Treatment = treatment)


#Calculate the sum of relative abundance for each treatment
df_sumed = df_joined %>%
  group_by(Treatment) %>%
  summarise(sum_abs = sum(abs))

#calculating the relative abundance of the group 
df_sumed$rel_abs = df_sumed$sum_abs/individuals_per_treatment

#generating the plot
Prevotella_bar <- ggplot(df_sumed, aes(Treatment, rel_abs, fill = Treatment)) +
  geom_col(color = "black") + # Plotting bars
  scale_fill_manual(values = c("Non-PD" = "sienna1", "PD-untreated" = "indianred2", 
                               "Entacapone" = "#9ACD32", "Pramipexole" = "#569866", 
                               "Rasagiline" = "#104E8B", "Amantadine" = "steelblue2")) + # Setting custom colors
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + # Adjusting expansion of y-axis
  theme_classic() + # Changing plot theme to a white background
  theme(
    axis.title = element_text(face = "bold", size = rel(2)), # Changing axis label font
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black", face = "bold", size = rel(2)), # Changing y axis text font
    legend.text = element_text(size = rel(1.2)), # Changing legend text font
    legend.title = element_text(face = "bold", size = rel(1.5)) # Changing legend title font
  ) +
  labs(x = NULL, y = "Abundance (normalized)", title = "") +
  ggtitle("Prevotella\n Abundance") +
  theme(plot.title = element_text(face = "bold", size = rel(2.5), hjust = 0.5)) 


####  Lachnospiraceae graph ####
taxa = c("g__Lachnospiraceae_UCG-001")
taxa_level = "Genus"
sub_beta_RA <- subset_taxa(pd_phyloseq, Genus == taxa)
sub_beta_order_RA <-tax_glom(sub_beta_RA, taxrank = taxa_level, NArm = FALSE)
ASV_df  =  data.frame(t(otu_table(sub_beta_order_RA)))
colnames(ASV_df) = "abs"
ASV_df$ID = rownames(ASV_df)
meta = data.frame(sample_data(sub_beta_order_RA))
meta$ID = rownames(meta)
df_joined = inner_join(meta, ASV_df , by = "ID")
individuals_per_treatment = c(nrow(filter(df_joined,df_joined$treatment == 1)),
                              nrow(filter(df_joined,df_joined$treatment == 2)),
                              nrow(filter(df_joined,df_joined$treatment == 3)),
                              nrow(filter(df_joined,df_joined$treatment == 4)),
                              nrow(filter(df_joined,df_joined$treatment == 5)),
                              nrow(filter(df_joined,df_joined$treatment == 6)))
df_joined$treatment = dplyr::recode(df_joined$treatment,
                                    "1" = "Non-PD",
                                    "2" = "PD-untreated",
                                    "3" = "Entacapone",
                                    "4" = "Pramipexole",
                                    "5" = "Rasagiline",
                                    "6" = "Amantadine")
df_joined <- df_joined <- df_joined %>%
  filter(treatment != 7) %>%
  rename(Treatment = treatment)
df_sumed = df_joined %>%
  group_by(Treatment) %>%
  summarise(sum_abs = sum(abs))
df_sumed$rel_abs = df_sumed$sum_abs/individuals_per_treatment
Lachnospiraceae_bar <- ggplot(df_sumed, aes(Treatment, rel_abs, fill = Treatment)) +
  geom_col(color = "black") + # Plotting bars
  scale_fill_manual(values = c("Non-PD" = "sienna1", "PD-untreated" = "indianred2", 
                               "Entacapone" = "#9ACD32", "Pramipexole" = "#569866", 
                               "Rasagiline" = "#104E8B", "Amantadine" = "steelblue2")) + # Setting custom colors
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + # Adjusting expansion of y-axis
  theme_classic() + # Changing plot theme to a white background
  theme(
    axis.title = element_text(face = "bold", size = rel(2)), # Changing axis label font
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black", face = "bold", size = rel(2)), # Changing y axis text font
    legend.text = element_text(size = rel(1.2)), # Changing legend text font
    legend.title = element_text(face = "bold", size = rel(1.5)) # Changing legend title font
  ) +
  labs(x = NULL, y = "Abundance (normalized)", title = "") +
  ggtitle("Lachnospiraceae_UCG-001\n Abundance") +
  theme(plot.title = element_text(face = "bold", size = rel(2.5), hjust = 0.5)) 



####  Akkermansia graph ####
taxa = c("g__Akkermansia")
taxa_level = "Genus"
sub_beta_RA <- subset_taxa(pd_phyloseq, Genus == taxa)
sub_beta_order_RA <-tax_glom(sub_beta_RA, taxrank = taxa_level, NArm = FALSE)
ASV_df  =  data.frame(t(otu_table(sub_beta_order_RA)))
colnames(ASV_df) = "abs"
ASV_df$ID = rownames(ASV_df)
meta = data.frame(sample_data(sub_beta_order_RA))
meta$ID = rownames(meta)
df_joined = inner_join(meta, ASV_df , by = "ID")
individuals_per_treatment = c(nrow(filter(df_joined,df_joined$treatment == 1)),
                              nrow(filter(df_joined,df_joined$treatment == 2)),
                              nrow(filter(df_joined,df_joined$treatment == 3)),
                              nrow(filter(df_joined,df_joined$treatment == 4)),
                              nrow(filter(df_joined,df_joined$treatment == 5)),
                              nrow(filter(df_joined,df_joined$treatment == 6)))
df_joined$treatment = dplyr::recode(df_joined$treatment,
                                    "1" = "Non-PD",
                                    "2" = "PD-untreated",
                                    "3" = "Entacapone",
                                    "4" = "Pramipexole",
                                    "5" = "Rasagiline",
                                    "6" = "Amantadine")
df_joined <- df_joined <- df_joined %>%
  filter(treatment != 7) %>%
  rename(Treatment = treatment)
df_sumed = df_joined %>%
  group_by(Treatment) %>%
  summarise(sum_abs = sum(abs))
df_sumed$rel_abs = df_sumed$sum_abs/individuals_per_treatment
Akkermansia_bar <- ggplot(df_sumed, aes(Treatment, rel_abs, fill = Treatment)) +
  geom_col(color = "black") + # Plotting bars
  scale_fill_manual(values = c("Non-PD" = "sienna1", "PD-untreated" = "indianred2", 
                               "Entacapone" = "#9ACD32", "Pramipexole" = "#569866", 
                               "Rasagiline" = "#104E8B", "Amantadine" = "steelblue2")) + # Setting custom colors
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + # Adjusting expansion of y-axis
  theme_classic() + # Changing plot theme to a white background
  theme(
    axis.title = element_text(face = "bold", size = rel(2)), # Changing axis label font
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black", face = "bold", size = rel(2)), # Changing y axis text font
    legend.text = element_text(size = rel(1.2)), # Changing legend text font
    legend.title = element_text(face = "bold", size = rel(1.5)) # Changing legend title font
  ) +
  labs(x = NULL, y = "Abundance (normalized)", title = "") +
  ggtitle("Akkermansia\n Abundance") +
  theme(plot.title = element_text(face = "bold", size = rel(2.5), hjust = 0.5)) 

Akkermansia_bar








####  Saving  ####
Bifidobacterium_bar
Prevotella_bar
Lachnospiraceae_bar
Akkermansia_bar

ggsave(filename="Bifidobacterium_levels.png", plot = Bifidobacterium_bar, width = 7, height = 6)
ggsave(filename="Prevotella_levels.png", plot = Prevotella_bar, width = 7, height = 6)
ggsave(filename="Lachnospiraceae_levels????.png", plot = Lachnospiraceae_bar, width = 7, height = 6)
ggsave(filename="Akkermansia_levels.png", plot = Akkermansia_bar, width = 7, height = 6)

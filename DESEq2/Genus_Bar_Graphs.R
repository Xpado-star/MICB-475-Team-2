library(tidyverse)
library(phyloseq)
library(ggplot2)
library(ggsignif)

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
ggplot(df_sumed, aes(Treatment, rel_abs, fill = Treatment)) +
  geom_col() + # Plotting bars
  scale_fill_manual(values = c("Non-PD" = "#104E8B", "PD-untreated" = "#256681", 
                               "Entacapone" = "#3A7E78", "Pramipexole" = "#569866", 
                               "Rasagiline" = "#78B24C", "Amantadine" = "#9ACD32")) + # Setting custom colors
  theme_classic() + # Changing plot theme to a white background
  theme(
    axis.title = element_text(face = "bold", size = rel(1.5)), # Changing axis label font
    axis.text.x = element_text(size = rel(1.5), angle = 45, hjust = 1), # Changing x axis text font
    axis.text.y = element_text(size = rel(1.5)), # Changing y axis text font
    legend.text = element_text(size = rel(1.2)), # Changing legend text font
    legend.title = element_text(size = rel(1.5)) # Changing legend title font
  ) +
  labs(x = NULL, y = "Abundance (normalized)", title = "") +
  ggtitle("Bifidobacterium \n Abundance") +
  theme(plot.title = element_text(face = "bold", size = rel(1.8), hjust = 0.5)) # Adjusting title size and centering



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
ggplot(df_sumed, aes(Treatment, rel_abs, fill = Treatment)) +
  geom_col() + # Plotting bars
  scale_fill_manual(values = c("Non-PD" = "#104E8B", "PD-untreated" = "#256681", 
                               "Entacapone" = "#3A7E78", "Pramipexole" = "#569866", 
                               "Rasagiline" = "#78B24C", "Amantadine" = "#9ACD32")) + # Setting custom colors
  theme_classic() + # Changing plot theme to a white background
  theme(
    axis.title = element_text(face = "bold", size = rel(1.5)), # Changing axis label font
    axis.text.x = element_text(size = rel(1.5), angle = 45, hjust = 1), # Changing x axis text font
    axis.text.y = element_text(size = rel(1.5)), # Changing y axis text font
    legend.text = element_text(size = rel(1.2)), # Changing legend text font
    legend.title = element_text(size = rel(1.5)) # Changing legend title font
  ) +
  labs(x = NULL, y = "Abundance (normalized)", title = "") +
  ggtitle("Prevotella \n Abundance") +
  theme(plot.title = element_text(face = "bold", size = rel(1.8), hjust = 0.5)) # Adjusting title size and centering




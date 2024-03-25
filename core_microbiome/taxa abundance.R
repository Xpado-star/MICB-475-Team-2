library(tidyverse)
library(phyloseq)
library(ggplot2)


#set your working dir
setwd()

##### load phyloseq object #####
load("phyloseq/pd_phyloseq.RData")


taxa = c("g__Bacteroides")
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
                              nrow(filter(df_joined,df_joined$treatment == 6)),
                              nrow(filter(df_joined,df_joined$treatment == 7)))
#change the treatments to meaningful information
df_joined$treatment = dplyr::recode(df_joined$treatment,
                                    "1" = "non parkinsons",
                                    "2" = "parkinsons",
                                    "3" = "entac",
                                    "4" = "prami",
                                    "5" = "rasag",
                                    "6" = "amant",
                                    "7" = "combo")


#Calculate the sum of relative abundance for each treatment
df_sumed = df_joined %>%
  group_by(treatment) %>%
  summarise(sum_abs = sum(abs))

#calculating the relative abundance of the group 
df_sumed$rel_abs = df_sumed$sum_abs/individuals_per_treatment


#generating the plot
ggplot(df_sumed,aes(treatment,rel_abs, fill = treatment)) +
  geom_col() + #plotting bars
  theme_classic() + #changeing plot theme to a white background
  theme(axis.title = element_text(face  = "bold", size = rel(1.2)), #Changing axis label font
        axis.text = element_text(face  = "bold", size = rel(1.2), angle = 90), #changing axis text font
        legend.text = element_text(face  = "bold", size = rel(1.2)), #changing legend text font
        legend.title = element_text(face  = "bold", size = rel(1.5))) +# changing legend title font 
  labs(x = "Treatments", y= "Bacteroides Abundance (normalized)") #Adding axis labels




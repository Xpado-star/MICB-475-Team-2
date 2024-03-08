setwd("D:/MICB-475/phyloseq")

library(writexl)
library(phyloseq)
library(tidyverse)
library(indicspecies)

load("pd_phyloseq.Rdata")




#filtering process

#for only groups 1 and 2 first to compare disease status.

df_1 = subset_samples(pd_phyloseq, treatment ==1 | treatment == 2)

#Then consider all groups including drug treatments (unfiltered)

df_2 = subset_samples(pd_phyloseq, treatment !=7)




# Run indicator species analysis using this is only groups 1&2 excluding all other groups
isa_df_1 <- multipatt(t(otu_table(df_1)), cluster = sample_data(df_1)$'treatment')
# Look at results
summary(isa_df_1)

#indicator species analysis for all groups?
isa_df_2 <- multipatt(t(otu_table(df_2)), cluster = sample_data(df_2)$'treatment')
# Look at results
summary(isa_df_2)

# Extract taxonomy table 
taxtable <- tax_table(pd_phyloseq) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# Merge taxonomy table with phyloseq object and filter by significant p-value
res <- isa_df_1$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value < 0.05) 

# Extract taxonomy table 
taxtable <- tax_table(pd_phyloseq) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# Merge taxonomy table with phyloseq object and filter by significant p-value
res <- isa_df_1$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value < 0.05) 

# View results
View(res)


# Merge taxonomy table with phyloseq object and filter by significant p-value
res2 <- isa_df_2$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value < 0.05) 

View(res2)

#making the table
group2_indi_table = res %>%
  filter(s.2 == 1)

group2_indi_table = group2_indi_table[,c(5,6,11,12)]
group2_indi_table$Genus = gsub("g__","",group2_indi_table$Genus)

write.csv(group2_indi_table, file = "Group1 indicator list.csv")





# create a list of dataframes to be written to separate sheets
df_list <- list(group1_indi_table, group2_indi_table)
names(df_list) <- c("Sheet1", "Sheet2")

# write the list to an Excel file
write_xlsx(df_list, "output.xlsx")

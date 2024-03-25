# Load necessary libraries
library(phyloseq)
library(ggplot2)
library(dplyr)

# Load the phyloseq object
load("../phyloseq/pd_phyloseq.RData")

# Filter the phyloseq object to include only Bifidobacterium genus
filtered_phyloseq <- subset_taxa(pd_phyloseq, Genus == "g__Bifidobacterium")

# Get the abundance matrix
genus_abundance <- as.data.frame(otu_table(filtered_phyloseq))
genus_abundance <- t(genus_abundance)
genus_abundance <- data.frame(genus_abundance, treatment = pd_phyloseq@sam_data[["treatment"]])
genus_abundance <- genus_abundance[genus_abundance$treatment != 7, ]
genus_abundance$treatment <- factor(genus_abundance$treatment, levels = c(1, 2, 3, 4, 5, 6), labels = c("Healthy", "PD Untreated", "Entacapone", "Pramipexole", "Rasagiline", "Amantadine"))
rownames(genus_abundance) <- NULL
genus_abundance <- genus_abundance %>%
  mutate(sum = rowSums(select(., -treatment), na.rm = TRUE))

Bifidobacterium_sums <- genus_abundance %>%
  group_by(treatment) %>%
  summarise(sum = sum(sum))

ggplot(Bifidobacterium_sums, aes(x = treatment, y = sum)) +
  geom_bar(stat = "identity", fill = "steelblue4") +
  labs(x = "Treatment", y = "Sum", title = "Bifidobacterium Abundances") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

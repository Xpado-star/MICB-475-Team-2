library(tidyverse)
library(phyloseq)
library(DESeq2)
library(patchwork)
library(ggrepel)
library(dplyr)
#### DESeq ####
# Convert phyloseq object to DESeq object
load("../phyloseq/pd_phyloseq.RData")

pd_plus1 <- transform_sample_counts(pd_phyloseq, function(x) x+1)
pd_deseq <- phyloseq_to_deseq2(pd_plus1, ~`treatment`)
DESEQ_pd <- DESeq(pd_deseq)


table(pd_phyloseq@sam_data[["treatment"]])

res_3_1 <- results(DESEQ_pd, tidy=TRUE, contrast = c("treatment","3","1"))
res_4_1 <- results(DESEQ_pd, tidy=TRUE, contrast = c("treatment","4","1"))
res_5_1 <- results(DESEQ_pd, tidy=TRUE, contrast = c("treatment","5","1"))
res_6_1 <- results(DESEQ_pd, tidy=TRUE, contrast = c("treatment","6","1"))

res_3_2 <- results(DESEQ_pd, tidy=TRUE, contrast = c("treatment","3","2"))
res_4_2 <- results(DESEQ_pd, tidy=TRUE, contrast = c("treatment","4","2"))
res_5_2 <- results(DESEQ_pd, tidy=TRUE, contrast = c("treatment","5","2"))
res_6_2 <- results(DESEQ_pd, tidy=TRUE, contrast = c("treatment","6","2"))

####  Volcano plots ####
vol_3_1 <- res_3_1 %>% 
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() + 
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) + 
  ggtitle("Group 3 vs 1") +
  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) + 
  theme(legend.position = "none")

vol_4_1 <- res_4_1 %>% 
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() + 
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) + 
  ggtitle("Group 4 vs 1") +
  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) + 
  theme(legend.position = "none")

vol_5_1 <- res_5_1 %>% 
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() + 
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) + 
  ggtitle("Group 5 vs 1") +
  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) + 
  theme(legend.position = "none")

vol_6_1 <- res_6_1 %>% 
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() + 
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) + 
  ggtitle("Group 6 vs 1") +
  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) + 
  theme(legend.position = "none")


vol_3_2 <- res_3_2 %>% 
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() + 
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) + 
  ggtitle("Group 3 vs 2") +
  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) + 
  theme(legend.position = "none")

vol_4_2 <- res_4_2 %>% 
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() + 
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) + 
  ggtitle("Group 4 vs 2") +
  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) + 
  theme(legend.position = "none")

vol_5_2 <- res_5_2 %>% 
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() + 
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) + 
  ggtitle("Group 5 vs 2") +
  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) + 
  theme(legend.position = "none")

vol_6_2 <- res_6_2 %>% 
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() + 
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) + 
  ggtitle("Group 6 vs 2") +
  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) + 
  theme(legend.position = "none")

wrap_plots(vol_3_1, vol_4_1, vol_5_1, vol_6_1, vol_3_2, vol_4_2, vol_5_2, vol_6_2, ncol = 4)

####  Bar Plots ####
#3_1
# To get table of results
sigASVs_3_1 <- res_3_1 %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row) %>%
  arrange(log2FoldChange) %>%
  group_by(sign(log2FoldChange)) %>%
  slice_head(n = 11)

# Get only asv names
sigASVs_vec_3_1 <- sigASVs_3_1 %>% pull(ASV)

# Prune phyloseq file
pd_DESeq_3_1 <- prune_taxa(sigASVs_vec_3_1,pd_phyloseq)
sigASVs_3_1 <- tax_table(pd_DESeq_3_1) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_3_1) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
sigASVs_3_1 <- sigASVs_3_1 %>%
  filter(!(is.na(Genus) | grepl("^NA\\.", Genus)))
view(sigASVs_3_1)

#3_1
#bar_3_1 <- ggplot(sigASVs_3_1) +
#  geom_bar(aes(x=Genus, y=log2FoldChange, fill = ifelse(log2FoldChange < 0, "blue", "red")), stat = "identity")+
#  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
#  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) + ggtitle("Group 3 vs 1") +
#  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) + theme(legend.position = "none")

bar_3_1 <- ggplot(sigASVs_3_1) +
  geom_bar(aes(x=log2FoldChange, y=Genus, fill = factor(log2FoldChange < 0)), stat = "identity") +
  geom_errorbar(aes(xmin=log2FoldChange-lfcSE, xmax=log2FoldChange+lfcSE, y=Genus)) +
  scale_fill_manual(values = c("TRUE" = "sienna1", "FALSE" = "#9ACD32")) +
  theme(axis.text.x = element_text(angle=0, hjust=1, vjust=0.5)) +
  ggtitle("Entacapone vs Non-PD") +
  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) +
  theme(legend.position = "none")

#4_1
sigASVs_4_1 <- res_4_1 %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row) %>%
  arrange(log2FoldChange) %>%
  group_by(sign(log2FoldChange)) %>%
  slice_head(n = 11)
sigASVs_vec_4_1 <- sigASVs_4_1 %>% pull(ASV)
pd_DESeq_4_1 <- prune_taxa(sigASVs_vec_4_1,pd_phyloseq)
sigASVs_4_1 <- tax_table(pd_DESeq_4_1) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_4_1) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
sigASVs_4_1 <- sigASVs_4_1 %>%
  filter(!(is.na(Genus) | grepl("^NA\\.", Genus)))
bar_4_1 <- ggplot(sigASVs_4_1) +
  geom_bar(aes(x=log2FoldChange, y=Genus, fill = factor(log2FoldChange < 0)), stat = "identity") +
  geom_errorbar(aes(xmin=log2FoldChange-lfcSE, xmax=log2FoldChange+lfcSE, y=Genus)) +
  scale_fill_manual(values = c("TRUE" = "sienna1", "FALSE" = "#569866")) +
  theme(axis.text.x = element_text(angle=0, hjust=1, vjust=0.5)) +
  ggtitle("Pramipexole vs Non-PD") +
  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) +
  theme(legend.position = "none")

#5_1
sigASVs_5_1 <- res_5_1 %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row) %>%
  arrange(log2FoldChange) %>%
  group_by(sign(log2FoldChange)) %>%
  slice_head(n = 11)
sigASVs_vec_5_1 <- sigASVs_5_1 %>% pull(ASV)
pd_DESeq_5_1 <- prune_taxa(sigASVs_vec_5_1,pd_phyloseq)
sigASVs_5_1 <- tax_table(pd_DESeq_5_1) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_5_1) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
sigASVs_5_1 <- sigASVs_5_1 %>%
  filter(!(is.na(Genus) | grepl("^NA\\.", Genus)))
bar_5_1 <- ggplot(sigASVs_5_1) +
  geom_bar(aes(x=log2FoldChange, y=Genus, fill = factor(log2FoldChange < 0)), stat = "identity") +
  geom_errorbar(aes(xmin=log2FoldChange-lfcSE, xmax=log2FoldChange+lfcSE, y=Genus)) +
  scale_fill_manual(values = c("TRUE" = "sienna1", "FALSE" = "#104E8B")) +
  theme(axis.text.x = element_text(angle=0, hjust=1, vjust=0.5)) +
  ggtitle("Rasagiline vs Non-PD") +
  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) +
  theme(legend.position = "none")

#6_1
sigASVs_6_1 <- res_6_1 %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row) %>%
  arrange(log2FoldChange) %>%
  group_by(sign(log2FoldChange)) %>%
  slice_head(n = 11)
sigASVs_vec_6_1 <- sigASVs_6_1 %>% pull(ASV)
pd_DESeq_6_1 <- prune_taxa(sigASVs_vec_6_1,pd_phyloseq)
sigASVs_6_1 <- tax_table(pd_DESeq_6_1) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_6_1) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
sigASVs_6_1 <- sigASVs_6_1 %>%
  filter(!(is.na(Genus) | grepl("^NA\\.", Genus)))
bar_6_1 <- ggplot(sigASVs_6_1) +
  geom_bar(aes(x=log2FoldChange, y=Genus, fill = factor(log2FoldChange < 0)), stat = "identity") +
  geom_errorbar(aes(xmin=log2FoldChange-lfcSE, xmax=log2FoldChange+lfcSE, y=Genus)) +
  scale_fill_manual(values = c("TRUE" = "sienna1", "FALSE" = "steelblue2")) +
  theme(axis.text.x = element_text(angle=0, hjust=1, vjust=0.5)) +
  ggtitle("Amantadine vs Non-PD") +
  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) +
  theme(legend.position = "none")

#3_2
sigASVs_3_2 <- res_3_2 %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row) %>%
  arrange(log2FoldChange) %>%
  group_by(sign(log2FoldChange)) %>%
  slice_head(n = 11)
sigASVs_vec_3_2 <- sigASVs_3_2 %>% pull(ASV)
pd_DESeq_3_2 <- prune_taxa(sigASVs_vec_3_2,pd_phyloseq)
sigASVs_3_2 <- tax_table(pd_DESeq_3_2) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_3_2) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
sigASVs_3_2 <- sigASVs_3_2 %>%
  filter(!(is.na(Genus) | grepl("^NA\\.", Genus)))
bar_3_2 <- ggplot(sigASVs_3_2) +
  geom_bar(aes(x=log2FoldChange, y=Genus, fill = factor(log2FoldChange < 0)), stat = "identity") +
  geom_errorbar(aes(xmin=log2FoldChange-lfcSE, xmax=log2FoldChange+lfcSE, y=Genus)) +
  scale_fill_manual(values = c("TRUE" = "indianred2", "FALSE" = "#9ACD32")) +
  theme(axis.text.x = element_text(angle=0, hjust=1, vjust=0.5)) +
  ggtitle("Entacapone vs PD-untreated") +
  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) +
  theme(legend.position = "none")

#4_2
sigASVs_4_2 <- res_4_2 %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row) %>%
  arrange(log2FoldChange) %>%
  group_by(sign(log2FoldChange)) %>%
  slice_head(n = 11)
sigASVs_vec_4_2 <- sigASVs_4_2 %>% pull(ASV)
pd_DESeq_4_2 <- prune_taxa(sigASVs_vec_4_2,pd_phyloseq)
sigASVs_4_2 <- tax_table(pd_DESeq_4_2) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_4_2) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
sigASVs_4_2 <- sigASVs_4_2 %>%
  filter(!(is.na(Genus) | grepl("^NA\\.", Genus)))
bar_4_2 <- ggplot(sigASVs_4_2) +
  geom_bar(aes(x=log2FoldChange, y=Genus, fill = factor(log2FoldChange < 0)), stat = "identity") +
  geom_errorbar(aes(xmin=log2FoldChange-lfcSE, xmax=log2FoldChange+lfcSE, y=Genus)) +
  scale_fill_manual(values = c("TRUE" = "indianred2", "FALSE" = "#569866")) +
  theme(axis.text.x = element_text(angle=0, hjust=1, vjust=0.5)) +
  ggtitle("Pramipexole vs PD-untreated") +
  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) +
  theme(legend.position = "none")

#5_2
sigASVs_5_2 <- res_5_2 %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row) %>%
  arrange(log2FoldChange) %>%
  group_by(sign(log2FoldChange)) %>%
  slice_head(n = 11)
sigASVs_vec_5_2 <- sigASVs_5_2 %>% pull(ASV)
pd_DESeq_5_2 <- prune_taxa(sigASVs_vec_5_2,pd_phyloseq)
sigASVs_5_2 <- tax_table(pd_DESeq_5_2) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_5_2) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
sigASVs_5_2 <- sigASVs_5_2 %>%
  filter(!(is.na(Genus) | grepl("^NA\\.", Genus)))
bar_5_2 <- ggplot(sigASVs_5_2) +
  geom_bar(aes(x=log2FoldChange, y=Genus, fill = factor(log2FoldChange < 0)), stat = "identity") +
  geom_errorbar(aes(xmin=log2FoldChange-lfcSE, xmax=log2FoldChange+lfcSE, y=Genus)) +
  scale_fill_manual(values = c("TRUE" = "indianred2", "FALSE" = "#104E8B")) +
  theme(axis.text.x = element_text(angle=0, hjust=1, vjust=0.5)) +
  ggtitle("Rasagiline vs PD-untreated") +
  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) +
  theme(legend.position = "none")

#6_2
sigASVs_6_2 <- res_6_2 %>%
  filter(padj < 0.01 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row) %>%
  arrange(log2FoldChange) %>%
  group_by(sign(log2FoldChange)) %>%
  slice_head(n = 11)
sigASVs_vec_6_2 <- sigASVs_6_2 %>% pull(ASV)
pd_DESeq_6_2 <- prune_taxa(sigASVs_vec_6_2,pd_phyloseq)
sigASVs_6_2 <- tax_table(pd_DESeq_6_2) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs_6_2) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))
sigASVs_6_2 <- sigASVs_6_2 %>%
  filter(!(is.na(Genus) | grepl("^NA\\.", Genus)))
bar_6_2 <- ggplot(sigASVs_6_2) +
  geom_bar(aes(x=log2FoldChange, y=Genus, fill = factor(log2FoldChange < 0)), stat = "identity") +
  geom_errorbar(aes(xmin=log2FoldChange-lfcSE, xmax=log2FoldChange+lfcSE, y=Genus)) +
  scale_fill_manual(values = c("TRUE" = "indianred2", "FALSE" = "steelblue2")) +
  theme(axis.text.x = element_text(angle=0, hjust=1, vjust=0.5)) +
  ggtitle("Amantadine vs PD-untreated") +
  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) +
  theme(legend.position = "none")


bar_3_1
bar_4_1
bar_5_1
bar_6_1
bar_3_2
bar_4_2
bar_5_2
bar_6_2
wrap_plots(bar_3_1, bar_4_1, bar_5_1, bar_6_1, bar_3_2, bar_4_2, bar_5_2, bar_6_2, ncol = 4)



ggsave(filename="bar_3_1.png",bar_3_1, width = 7, height = 10)
ggsave(filename="bar_4_1.png",bar_4_1, width = 7, height = 10)
ggsave(filename="bar_5_1.png",bar_5_1, width = 7, height = 10)
ggsave(filename="bar_6_1.png",bar_6_1, width = 7, height = 10)
ggsave(filename="bar_3_2.png",bar_3_2, width = 7, height = 10)
ggsave(filename="bar_4_2.png",bar_4_2, width = 7, height = 10)
ggsave(filename="bar_5_2.png",bar_5_2, width = 7, height = 10)
ggsave(filename="bar_6_2.png",bar_6_2, width = 7, height = 10)

ggsave(filename="vol_3_1.png",vol_3_1)
ggsave(filename="vol_4_1.png",vol_4_1)
ggsave(filename="vol_5_1.png",vol_5_1)
ggsave(filename="vol_6_1.png",vol_6_1)
ggsave(filename="vol_3_2.png",vol_3_2)
ggsave(filename="vol_4_2.png",vol_4_2)
ggsave(filename="vol_5_2.png",vol_5_2)
ggsave(filename="vol_6_2.png",vol_6_2)


write.csv(sigASVs_3_1, file = "ASVs_3_1.csv", row.names = FALSE)
write.csv(sigASVs_4_1, file = "ASVs_4_1.csv", row.names = FALSE)
write.csv(sigASVs_5_1, file = "ASVs_5_1.csv", row.names = FALSE)
write.csv(sigASVs_6_1, file = "ASVs_6_1.csv", row.names = FALSE)
write.csv(sigASVs_3_2, file = "ASVs_3_2.csv", row.names = FALSE)
write.csv(sigASVs_4_2, file = "ASVs_4_2.csv", row.names = FALSE)
write.csv(sigASVs_5_2, file = "ASVs_5_2.csv", row.names = FALSE)
write.csv(sigASVs_6_2, file = "ASVs_6_2.csv", row.names = FALSE)

####  If we want labels on the volcano plots  ####
vol_3_1 <- res_3_1 %>% 
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 2) %>%
  ggplot() + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj), col = significant)) +
  geom_text_repel(data = sigASVs_3_1, aes(x = log2FoldChange, y = -log10(padj), label = Genus), 
                  box.padding = 0.5, point.padding = 0.5) + 
  ggtitle("Group 3 vs 1") +
  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) + 
  theme(legend.position = "none")

vol_4_1 <- res_4_1 %>% 
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 2) %>%
  ggplot() + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj), col = significant)) +
  geom_text_repel(data = sigASVs_4_1, aes(x = log2FoldChange, y = -log10(padj), label = Genus), 
                  box.padding = 0.5, point.padding = 0.5) + 
  ggtitle("Group 4 vs 1") +
  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) + 
  theme(legend.position = "none")

vol_5_1 <- res_5_1 %>% 
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 2) %>%
  ggplot() + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj), col = significant)) +
  geom_text_repel(data = sigASVs_5_1, aes(x = log2FoldChange, y = -log10(padj), label = Genus), 
                  box.padding = 0.5, point.padding = 0.5) + 
  ggtitle("Group 5 vs 1") +
  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) + 
  theme(legend.position = "none")

vol_6_1 <- res_6_1 %>% 
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 2) %>%
  ggplot() + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj), col = significant)) +
  geom_text_repel(data = sigASVs_6_1, aes(x = log2FoldChange, y = -log10(padj), label = Genus), 
                  box.padding = 0.5, point.padding = 0.5) + 
  ggtitle("Group 6 vs 1") +
  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) + 
  theme(legend.position = "none")

vol_3_2 <- res_3_2 %>% 
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 2) %>%
  ggplot() + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj), col = significant)) +
  geom_text_repel(data = sigASVs_3_2, aes(x = log2FoldChange, y = -log10(padj), label = Genus), 
                  box.padding = 0.5, point.padding = 0.5) + 
  ggtitle("Group 3 vs 2") +
  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) + 
  theme(legend.position = "none")

vol_4_2 <- res_4_2 %>% 
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 2) %>%
  ggplot() + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj), col = significant)) +
  geom_text_repel(data = sigASVs_4_2, aes(x = log2FoldChange, y = -log10(padj), label = Genus), 
                  box.padding = 0.5, point.padding = 0.5) + 
  ggtitle("Group 4 vs 2") +
  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) + 
  theme(legend.position = "none")

vol_5_2 <- res_5_2 %>% 
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 2) %>%
  ggplot() + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj), col = significant)) +
  geom_text_repel(data = sigASVs_5_2, aes(x = log2FoldChange, y = -log10(padj), label = Genus), 
                  box.padding = 0.5, point.padding = 0.5) + 
  ggtitle("Group 5 vs 2") +
  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) + 
  theme(legend.position = "none")

vol_6_2 <- res_6_2 %>% 
  mutate(significant = padj < 0.01 & abs(log2FoldChange) > 2) %>%
  ggplot() + 
  geom_point(aes(x = log2FoldChange, y = -log10(padj), col = significant)) +
  geom_text_repel(data = sigASVs_6_2, aes(x = log2FoldChange, y = -log10(padj), label = Genus), 
                  box.padding = 0.5, point.padding = 0.5) + 
  ggtitle("Group 6 vs 2") +
  theme(plot.title = element_text(family = "Arial", size = 20, hjust = 0.5, face = "bold")) + 
  theme(legend.position = "none")

wrap_plots(vol_3_1, vol_4_1, vol_5_1, vol_6_1)
wrap_plots(vol_3_2, vol_4_2, vol_5_2, vol_6_2)
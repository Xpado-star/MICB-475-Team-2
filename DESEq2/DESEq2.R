library(tidyverse)
library(phyloseq)
library(DESeq2)

#### DESeq ####
# Convert phyloseq object to DESeq object
load(".../.../phyloseq/pd_phyloseq.RData")

pd_plus1 <- transform_sample_counts(pd_phyloseq, function(x) x+1)
pd_deseq <- phyloseq_to_deseq2(pd_plus1, ~`treatment`)
DESEQ_pd <- DESeq(pd_deseq)


table(pd_phyloseq@sam_data[["treatment"]])

res_1_3 <- results(DESEQ_pd, tidy=TRUE, contrast = c("treatment","1","3"))
res_1_4 <- results(DESEQ_pd, tidy=TRUE, contrast = c("treatment","1","4"))
res_1_5 <- results(DESEQ_pd, tidy=TRUE, contrast = c("treatment","1","5"))
res_1_6 <- results(DESEQ_pd, tidy=TRUE, contrast = c("treatment","1","6"))

res_2_3 <- results(DESEQ_pd, tidy=TRUE, contrast = c("treatment","2","3"))
res_2_4 <- results(DESEQ_pd, tidy=TRUE, contrast = c("treatment","2","4"))
res_2_5 <- results(DESEQ_pd, tidy=TRUE, contrast = c("treatment","2","5"))
res_2_6 <- results(DESEQ_pd, tidy=TRUE, contrast = c("treatment","2","6"))

View(res_1_3)
View(res_1_4)
View(res_1_5)
View(res_1_6)

View(res_2_3)
View(res_2_4)
View(res_2_5)
View(res_2_6)

ggplot(res_1_3) + geom_point(aes(x=log2FoldChange, y=-log10(padj)))
vol_plot <- res %>% mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

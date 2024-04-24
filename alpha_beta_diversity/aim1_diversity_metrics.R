#### ALPHA AND BETA DIVERSITY METRICS ####

### load packaages ###
library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(vegan)
library(dunn.test)
library(FSA)
library(ggstatsplot)
library(ggthemes)

### load phyloseq object ###
load("../phyloseq/pd_rare.RData")


pd_rare_combo <- pd_rare

# remove group 7 treatment (combo)
pd_rare <- subset_samples(pd_rare, treatment != 7)


#### ALPHA DIVERSITY ####

# view initial diversity metrics
plot_richness(pd_rare)

# plot as box plots
gg_richness <- plot_richness(pd_rare, x = "treatment") +
  xlab("Treatment Group") +
  geom_boxplot() + theme_bw()

# view graph
gg_richness

# save plot
ggsave(filename = "plot_richness_all.png"
       , gg_richness
       , height=6, width=8)

# plot Shannon only
gg_richness_shannon <- plot_richness(pd_rare, x = "treatment", measures = "Shannon") +
  xlab("Treatment Group") +
  geom_boxplot() + 
  theme_bw()
gg_richness_shannon

# save plot
ggsave(filename = "plot_richness_shannon.png"
       , gg_richness_shannon
       , height=6, width=8)

### ALPHA DIVERSITY STATISTICAL ANALYSIS ###
# create dataframe with Shannon's index
shannon_div <- estimate_richness(pd_rare, split = TRUE, measures='Shannon')
samp_dat_wdiv <- data.frame(sample_data(pd_rare), shannon_div)



#re-label treatment values for final plot
samp_dat_wdiv$treatment = dplyr::recode(samp_dat_wdiv$treatment,
                                  "1" = "Healthy",
                                  "2" = "PD-untreated",
                                  "3" = "Entacapone",
                                  "4" = "Pramipexole",
                                  "5" = "Rasagiline",
                                  "6" = "Amantadine")


# plot distributions to determine statistical test
shannon_distributions <- ggplot(samp_dat_wdiv, aes(x = Shannon)) +
  geom_histogram() +
  labs(title = "Distribution of Shannon Index by Treatment",
       x = "Shannon index",
       y = "Frequency") +
  facet_wrap(~treatment, scales = "free")

# view histograms
shannon_distributions

# save plots
ggsave(filename = "shannon_distributions_hist.png"
       , shannon_distributions
       , height=6, width=8)

# create log values and repeat
samp_dat_wdiv <- mutate(samp_dat_wdiv, log_Shannon = log(samp_dat_wdiv$Shannon))

shannon_distribution_log <- ggplot(samp_dat_wdiv, aes(x = log_Shannon)) +
  geom_histogram() +
  labs(title = "Distribution of Shannon Index by Treatment",
       x = "Shannon index",
       y = "Frequency") +
  facet_wrap(~treatment, scales = "free")

# view histogram
shannon_distribution_log


#save plots
ggsave(filename = "shannon_distributions_hist_log.png"
       , shannon_distribution_log
       , height=6, width=8)

# Kruskal-Wallis test
kruskal_result <- kruskal.test(Shannon ~ treatment, data = samp_dat_wdiv)
# Kruskal-Wallis chi-squared = 7.6084, df = 5, p-value = 0.1792
# P-value > 0.05 - fail to reject null

# Adjusted p-value
p.adjust(kruskal_result$p.value, method = "fdr")
# [1] 0.1791794

# Perform Dunn Test for post-hoc comparisons
dunn.test(samp_dat_wdiv$Shannon, g = samp_dat_wdiv$treatment, method = "bonferroni")
#none of the pairwise comparisons for any groups is significant


#View mean and confidence intervals for each treatment
aggregate(Shannon ~ treatment, data = samp_dat_wdiv, FUN = mean)

table(samp_dat_wdiv$treatment)

data_sub <- subset.data.frame(samp_dat_wdiv, treatment == 'Amantadine')
table(data_sub$treatment)
t.test(data_sub$Shannon)$conf.int

# color reference for plots
colors <- c("sienna1","indianred2", "#9ACD32", "#569866", "#104E8B", "steelblue2")

# create final box plot
shannons_final_comparison <- ggbetweenstats(
  data = samp_dat_wdiv,
  x = treatment,
  y = Shannon,
  ylab = "Shannon's diversity index",
  xlab = "Dopaminergic therapeutic",
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = TRUE,
  results.subtitle = FALSE
) + ggplot2::scale_color_manual(values = c("sienna1","indianred2", "#9ACD32", "#569866", "#104E8B", "steelblue2")) +
  theme_clean(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "none")
  
shannons_final_comparison


ggsave(filename = "shannon_alpha_with_p.png"
       , shannons_final_comparison
       , height=8, width=9)


# Faith's PD

# calculate Faith's PD
phylo_dist <- pd(t(otu_table(pd_rare)), phy_tree(pd_rare),
                 include.root=F) 

#add to data frame
samp_dat_wdiv$PD <- phylo_dist$PD

# plot distributions to determine statistical test
faith_distributions <- ggplot(samp_dat_wdiv, aes(x = PD)) +
  geom_histogram() +
  labs(title = "Distribution of Faith's PD by Treatment",
       x = "PD",
       y = "Frequency") +
  facet_wrap(~treatment, scales = "free")

# view histograms
faith_distributions
# non-normal distribution --> nonparametric test needed

# Kruskal-Wallis test
kruskal_result_faith <- kruskal.test(PD ~ treatment, data = samp_dat_wdiv)
p.adjust(kruskal_result_faith$p.value, method = "fdr")
# [1] 0.06368479
# result was not significant

# Perform Dunn Test for post-hoc comparisons
dunn.test(samp_dat_wdiv$PD, g = samp_dat_wdiv$treatment, method = "bonferroni")

# create final box plot
faiths_final_comparison <- ggbetweenstats(
  data = samp_dat_wdiv,
  x = treatment,
  y = PD,
  ylab = "Faith's phylogenetic diversity",
  xlab = "Dopaminergic therapeutic",
  type = "nonparametric", # ANOVA or Kruskal-Wallis
  plot.type = "box",
  pairwise.comparisons = TRUE,
  pairwise.display = "significant",
  centrality.plotting = TRUE,
  results.subtitle = FALSE
) + ggplot2::scale_color_manual(values = c("sienna1","indianred2", "#9ACD32", "#569866", "#104E8B", "steelblue2")) +
  theme_clean(base_size = 18) +
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "none")

faiths_final_comparison


ggsave(filename = "faith_alpha_final.png"
       , faiths_final_comparison
       , height=8, width=9)



### BETA DIVERSITY ###

# create ordination object
distance_method <- phyloseq::distance(pd_rare, method= 'wunifrac')
pcoa_wunifrac <- ordinate(pd_rare, method="PCoA", distance= distance_method)

# plot using ggplot
gg_pcoa <- plot_ordination(pd_rare, pcoa_wunifrac, color = "treatment") +
  labs(color = "Treatment") +
  stat_ellipse() +
  theme_clean() +
  scale_color_manual(values = c("sienna1","indianred2", "#9ACD32", "#569866", "#104E8B", "steelblue2"),
                     labels = c("Healthy", "PD-untreated", "Entacapone", "Pramipexole", "Rasagiline", "Amantadine")) +
  theme(panel.grid.major.x = element_line(color = "grey", size = 0.4, linetype = "dotted"))
gg_pcoa

?plot_ordination

#save plot
ggsave(filename = "beta_diversity_treatment_wunifrac.png"
       , gg_pcoa
       , height=8, width=12)


### BETA DIVERSITY STATISTICAL ANALYSIS ###

#weighted unifrac
dm_unifrac <- UniFrac(pd_rare, weighted=TRUE)
data_frame_permanova <- data.frame(sample_data(pd_rare))
perm_result <- adonis2(dm_unifrac ~ treatment, data=data_frame_permanova, permutations = 10000)
#             Df SumOfSqs      R2 F      Pr(>F)  
# treatment   5  0.03014 0.03166 1.7853 0.05459
perm_result
p.adjust(perm_result$`Pr(>F)`, method = "fdr")
# [1] 0.05459454         NA         NA



### BETA DIVERSITY PAIRWISE ###

#example - groups 1 and 6 --> edit the line below to view all other comparisons
pd_subset <- subset_samples(pd_rare, treatment == 1 | treatment == 2)

distance_method_sub <- phyloseq::distance(pd_subset, method= 'wunifrac')
pcoa_wunifrac_sub <- ordinate(pd_subset, method="PCoA", distance= distance_method_sub)
gg_pcoa_sub <- plot_ordination(pd_subset, pcoa_wunifrac_sub, color = "treatment") +
  labs(color = "Treatment") +
  stat_ellipse()
gg_pcoa_sub

# PERMANOVA
dm_unifrac_sub <- UniFrac(pd_subset, weighted=TRUE)
data_frame_permanova_sub <- data.frame(sample_data(pd_subset))
adonis2(dm_unifrac_sub ~ treatment, data=data_frame_permanova_sub, permutations = 10000)



### TAXONOMIC BAR PLOTS ###

# Plot bar plot of taxonomy
plot_bar(pd_rare, fill="Phylum")

# Convert to relative abundance
pd_RA <- transform_sample_counts(pd_rare, function(x) x/sum(x))

# To remove black bars, "glom" by phylum first
pd_phylum <- tax_glom(pd_RA, taxrank = "Phylum", NArm=FALSE)

gg_taxa <- plot_bar(pd_phylum, fill="Phylum") + 
  facet_wrap(~treatment, scales = "free_x") +
  theme_clean()

gg_taxa

ggsave("plot_taxonomy.png"
       , gg_taxa
       , height=8, width =12)


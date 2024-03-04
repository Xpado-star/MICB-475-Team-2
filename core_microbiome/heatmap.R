# Load data
library(microbiome)
load("~/Desktop/MICB-475-Team-2/phyloseq/pd_phyloseq.RData")

# Calculate compositional version of the data
# (relative abundances)
pdseq.rel <- microbiome::transform(pd_phyloseq, "compositional")

head(prevalence(pdseq.rel, detection = 1/100, sort = TRUE))
head(prevalence(pdseq.rel, detection = 1/100, sort = TRUE, count = TRUE))

# Convert to relative abundance
pd_RA <- transform_sample_counts(pdseq.rel, fun=function(x) x/sum(x))

# Filter dataset by drug use
pd_entacapone <- subset_samples(pd_RA, `entacapone`=="1")

# Core with compositionals:
library(RColorBrewer)
library(reshape)

prevalences <- seq(.05, 1, .05)

detections <- round(10^seq(log10(0.01), log10(.2), length = 9), 3)

# Also define gray color palette
gray <- gray(seq(0,1,length=5))

#Added pseq.rel, I thin... must be checked if it was in the the rednred version,; where it is initialized
#pseq.rel<- microbiome::transform(pseq, 'compositional')
#min-prevalence gets the 100th highest prevalence
p <- plot_core(pd_entacapone,
               plot.type = "heatmap", 
               colours = gray,
               prevalences = prevalences, 
               detections = detections, 
               min.prevalence = prevalence(pd_entacapone, sort = TRUE)[100]) +
  labs(x = "Detection Threshold\n(Relative Abundance (%))") +
  
  #Adjusts axis text size and legend bar height
  theme(axis.text.y= element_text(size=8, face="italic"),
        axis.text.x.bottom=element_text(size=8),
        axis.title = element_text(size=10),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10))

print(p)

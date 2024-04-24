### BACTERIAL ABUNDANCE REGRESSION ANALYSIS ###

#load libraries
library(tidyverse)
library(phyloseq)
library(ggplot2)
library(ggthemes)
library(ggpubr)

#### Bifidobacterium regression ####

#load metadata file
load("bifidobac_abs.RData")
assign("bifidobac_data", df_joined)

#remove NA values from UPDRS score
bifidobac_data <- bifidobac_data[!is.na(bifidobac_data$updrstotal),]
options(max.print = .Machine$integer.max)

#re-label Treatment values
bifidobac_data$Treatment = dplyr::recode(bifidobac_data$Treatment,
                                    "2" = "PD-untreated",
                                    "3" = "Entacapone",
                                    "4" = "Pramipexole",
                                    "5" = "Rasagiline",
                                    "6" = "Amantadine")


# create a scatter plot to view relationship
ggplot(bifidobac_data, aes(x=bifidobac_data$abs, 
                      y=bifidobac_data$updrstotal)) +
  geom_point() +
  theme_minimal() +
  xlab("Bifidobacterium abundance") +
  ylab("UPDRS score")

# some abundance values are >1000 and does not allow for a proper scatterplot.
# To fix this, we will log-transformation the abundance values

# create a new column of log(abundance + 1)
bifidobac_data$log_abs <- log(bifidobac_data$abs + 1)


# add colors based on Treatment type
colors <- c("indianred2", "#9ACD32", "#569866", "#104E8B", "steelblue2")

#re-plot graph
bifidobac_plot <- ggplot(bifidobac_data, aes(x=bifidobac_data$log_abs, 
                           y=bifidobac_data$updrstotal,
                           color = Treatment)) +
  geom_point() +
  theme_clean() +
  xlab("Log Bifidobacterium abundance") +
  ylab("UPDRS score") + 
  geom_smooth(aes(group = 1), method = "lm", color = "black", fill = "lightgrey") +
  scale_color_manual(values = colors) + 
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold")) +
  stat_cor(aes(group = 1), label.x = 0, label.y = 115) +
  stat_regline_equation(aes(group = 1), label.x = 0, label.y = 120)


# view plot
bifidobac_plot

#view plots based on treatment type
bifidobac_plot + facet_wrap(~Treatment, scales = "free_y")


#### ANALYSIS ####
#Use lm() for linear regression analysis
bifidobac_regression <- lm(updrstotal ~ log_abs, data = bifidobac_data)

# View regression equation
bifidobac_regression

# View p-value
summary(bifidobac_regression)
# 0.809

#residual plot
plot(residuals(bifidobac_regression) ~ log_abs, data = bifidobac_data)
abline(h=0)

plot(predict(bifidobac_regression), bifidobac_data$updrstotal, 
     xlab = "Predicted Values", 
     ylab = "Observed Values") 

bif_glm <- glm(updrstotal ~ log_abs, family = gaussian(link = "identity"), data = bifidobac_data)
summary(bif_glm)


                              #### PREVOTELLA REGRESSION ####

#load metadata file
load("prevotella_abs.RData")

#remove NA values from UPDRS score
prevotella_data <- prevotella_data[!is.na(prevotella_data$updrstotal),]

# create a scatter plot to view relationship
ggplot(prevotella_data, aes(x=prevotella_data$abs, 
                           y=prevotella_data$updrstotal)) +
  geom_point() +
  theme_minimal() +
  xlab("Prevotella abundance") +
  ylab("UPDRS score")

# create a new column of log(abundance + 1)
prevotella_data$log_abs <- log(prevotella_data$abs + 1)


#re-plot graph
prevotella_plot <- ggplot(prevotella_data, aes(x=prevotella_data$log_abs, 
                                             y=prevotella_data$updrstotal,
                                             color = Treatment)) +
  geom_point() +
  theme_clean() +
  xlab("Log Prevotella abundance") +
  ylab("UPDRS score") + 
  geom_smooth(aes(group = 1), method = "lm", color = "black", fill = "lightgrey") +
  scale_color_manual(values = colors) + 
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold")) +
  stat_cor(aes(group = 1), label.x = 0, label.y = 115) +
  stat_regline_equation(aes(group = 1), label.x = 0, label.y = 120)

#view plot
prevotella_plot

#view plots based on treatment type
prevotella_plot + facet_wrap(~Treatment, scales = "free_y")



#### LACHNOSPIRACEAE REGRESSION ####

#load metadata file
load("lachnos_abs.RData")

#remove NA values from UPDRS score
lachnos_data <- lachnos_data[!is.na(lachnos_data$updrstotal),]

# create a scatter plot to view relationship
ggplot(lachnos_data, aes(x=lachnos_data$abs, 
                            y=lachnos_data$updrstotal)) +
  geom_point() +
  theme_minimal() +
  xlab("Lachnospiraceae abundance") +
  ylab("UPDRS score")

# create a new column of log(abundance + 1)
lachnos_data$log_abs <- log(lachnos_data$abs + 1)


#re-plot graph
lachnos_plot <- ggplot(lachnos_data, aes(x=lachnos_data$log_abs, 
                                               y=lachnos_data$updrstotal,
                                               color = Treatment)) +
  geom_point() +
  theme_clean() +
  xlab("Log Lachnospiraceae abundance") +
  ylab("UPDRS score") + 
  geom_smooth(aes(group = 1), method = "lm", color = "black", fill = "lightgrey") +
  scale_color_manual(values = colors) + 
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold")) +
  stat_cor(aes(group = 1), label.x = 0, label.y = 110) +
  stat_regline_equation(aes(group = 1), label.x = 0, label.y = 120)

#view plot
lachnos_plot

#view plots based on treatment type
lachnos_plot + facet_wrap(~Treatment, scales = "free_y") +
  theme(strip.text.x = element_text(face = "bold"))


#### AKKERMANSIA REGRESSION ###

#load metadata file
load("akker_abs.RData")

#remove NA values from UPDRS score
akker_data <- akker_data[!is.na(akker_data$updrstotal),]

# create a scatter plot to view relationship
ggplot(akker_data, aes(x=akker_data$updrstotal, 
                         y=akker_data$log_abs)) +
  geom_point() +
  theme_minimal() +
  xlab("Lachnospiraceae abundance") +
  ylab("UPDRS score")

# create a new column of log(abundance + 1)
akker_data$log_abs <- log(akker_data$abs + 1)


#re-plot graph
akker_plot <- ggplot(akker_data, aes(x=akker_data$log_abs, 
                                         y=akker_data$updrstotal,
                                         color = Treatment)) +
  geom_point() +
  theme_clean() +
  xlab("Log Akkermansia abundance") +
  ylab("UPDRS score") + 
  geom_smooth(aes(group = 1), method = "lm", color = "black", fill = "lightgrey") +
  scale_color_manual(values = colors) + 
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold")) +
  stat_cor(aes(group = 1), label.x = 0, label.y = 115) +
  stat_regline_equation(aes(group = 1), label.x = 0, label.y = 120)

#view plot
akker_plot

#view plots based on treatment type
akker_plot + facet_wrap(~Treatment, scales = "free_y")



### Faecalibacterium regression ###

taxa = c("g__Faecalibacterium")
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

faecalibac_data <- df_joined

#remove NA values from UPDRS score
faecalibac_data <- faecalibac_data[!is.na(faecalibac_data$updrstotal),]

# create a new column of log(abundance + 1)
faecalibac_data$log_abs <- log(faecalibac_data$abs + 1)

# plot graph
faecalibac_plot <- ggplot(faecalibac_data, aes(x=faecalibac_data$abs, 
                                     y=faecalibac_data$updrstotal,
                                     color = Treatment)) +
  geom_point() +
  theme_clean() +
  xlab("Log Faecalibacterium abundance") +
  ylab("UPDRS score") + 
  geom_smooth(aes(group = 1), method = "lm", color = "black", fill = "lightgrey") +
  scale_color_manual(values = colors) + 
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold")) +
  stat_cor(aes(group = 1), label.x = 0, label.y = 115) +
  stat_regline_equation(aes(group = 1), label.x = 0, label.y = 120)

#view plot
faecalibac_plot

#view plots based on treatment type
faecalibac_plot + facet_wrap(~Treatment, scales = "free_y")

### Bacteroides regression ###

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

bacter_data <- df_joined

#remove NA values from UPDRS score
bacter_data <- bacter_data[!is.na(bacter_data$updrstotal),]

# create a new column of log(abundance + 1)
bacter_data$log_abs <- log(bacter_data$abs + 1)


#re-plot graph
bacter_plot <- ggplot(bacter_data, aes(x=bacter_data$abs, 
                                               y=bacter_data$updrstotal,
                                               color = Treatment)) +
  geom_point() +
  theme_clean() +
  xlab("Log Faecalibacterium abundance") +
  ylab("UPDRS score") + 
  geom_smooth(aes(group = 1), method = "lm", color = "black", fill = "lightgrey") +
  scale_color_manual(values = colors) + 
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold")) +
  stat_cor(aes(group = 1), label.x = 0, label.y = 115) +
  stat_regline_equation(aes(group = 1), label.x = 0, label.y = 120)

#view plot
bacter_plot

#view plots based on treatment type
bacter_plot + facet_wrap(~Treatment, scales = "free_y")

############ Adding columns to PD metadata tsv ###########

library(tidyverse)



#load in data

pd_data <- read_delim(file = "../../qiime_files/parkinsons_export/parkinsons_metadata.txt", delim = "\t")

summary(pd_data)

#convert all NAs to 0s
pd_data_zero <- pd_data %>% 
  mutate(entacapone = ifelse(is.na(entacapone), 0, entacapone),
                    pramipexole = ifelse(is.na(pramipexole), 0, pramipexole),
                    rasagiline = ifelse(is.na(rasagiline), 0, rasagiline),
                    amantadine = ifelse(is.na(amantadine), 0, amantadine))

#Add a new column called 'treatment' with the following:
# control = 1
# pd_untreated = 2
# pd_entac = 3 
# pd_prami = 4
# pd_rasag = 5
# pd_amant = 6
# pd_combo = 7


pd_data_treatment <- pd_data_zero %>% 
  mutate(treatment = ifelse(Disease == "Control", 1,
                            ifelse(Disease == "PD" & (entacapone + pramipexole + rasagiline + amantadine) == 0, 2,
                                   ifelse(Disease == "PD" & entacapone == 1 & (pramipexole + rasagiline + amantadine) == 0, 3,
                                          ifelse(Disease == "PD" & pramipexole == 1 & (entacapone + rasagiline + amantadine) == 0, 4,
                                                 ifelse(Disease == "PD" & rasagiline == 1& (entacapone + pramipexole + amantadine) == 0, 5,
                                                        ifelse(Disease == "PD" & amantadine == 1 & (entacapone + pramipexole + rasagiline) == 0, 6,
                                                               ifelse(Disease == "PD" & (entacapone + pramipexole + rasagiline + amantadine) > 1, 7, NA))))))))

pd_data_treatment <- arrange(pd_data_treatment, treatment)
pd_data_treatment$treatment <- as.factor(pd_data_treatment$treatment)
class(pd_data_treatment$treatment)

freq_table <- table(pd_data_treatment$treatment)
freq_table

save(pd_data_treatment, file = "pd_metadata_treatment.RData")
write.table(pd_data_treatment, file="pd_metadata_treatment.tsv", sep="\t", quote=FALSE, row.names = FALSE)

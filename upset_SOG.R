library(UpSetR)
library(tidyverse)
library(dplyr)
library(magrittr)
library(tidyr)
library(readr)
library(stringr)

#### LOAD DATA ####

# Sog vs Mývatn 
SOG_MYV = read.delim('../stacks_populations_result/iceland.recode.p.fst_MYV-SOG.tsv')
# Sog vs Svínavatn 
SOG_SVI = read.delim('../stacks_populations_result/iceland.recode.p.fst_SOG-SVI.tsv')
# Sog vs Hestvatn 
SOG_HES = read.delim('../stacks_populations_result/iceland.recode.p.fst_HES-SOG.tsv')
# Sog vs Öxara
SOG_OXA = read.delim('../stacks_populations_result/iceland.recode.p.fst_OXA-SOG.tsv')
# Sog vs STF 
SOG_STF = read.delim('../stacks_populations_result/iceland.recode.p.fst_SOG-STF.tsv')



#### MERGE DATA ####

one_data <- rbind(SOG_MYV,SOG_SVI,SOG_HES,SOG_OXA,SOG_STF)



one_data$Location_Merge <- paste(one_data$Chr,one_data$BP,sep="_")
one_data$Population_Merge <- paste(one_data$Pop.1.ID,one_data$Pop.2.ID,sep="_")

one_data$SIGNIFICANT <- 0
# change Fst value
significant_value <- 0.5
one_data$SIGNIFICANT[one_data$AMOVA.Fst>significant_value] <- 1




#### SUBSET DATA #####


one_data_subset <- data.frame(Location = one_data$Location_Merge,
                              FST = one_data$SIGNIFICANT,
                              Populations = one_data$Population_Merge)

# remove second duplicate
one_data_subset_duplicates <- distinct(one_data_subset, Location, Populations, .keep_all= TRUE)

# spurja AP
one_data_subset_duplicates_wide <- one_data_subset_duplicates %>%
  pivot_wider(names_from = Populations, values_from = FST)


one_data_subset_wide_final <- one_data_subset_duplicates_wide[complete.cases(one_data_subset_duplicates_wide), ]


one_data_subset_wide_final <- as.data.frame(one_data_subset_wide_final)
one_data_subset_wide_final$Location <- factor(one_data_subset_wide_final$Location)
one_data_subset_wide_final$MYV_SOG <- as.integer(one_data_subset_wide_final$MYV_SOG)
one_data_subset_wide_final[-1] <- sapply(one_data_subset_wide_final[-1], as.integer)



str(one_data_subset_wide_final)




#### LOAD DATA #######


png('../comparisons_SOG/test_upset.png',width = 1000, height = 600,units = "px")
upset(one_data_subset_wide_final,
      #sets =names(color),
      #sets.bar.color = color,
      order.by = "freq",
      text.scale = c(2.5, # size of y-axis title
                     1.5, # size of values on y-axis
                     2, # size of x-axis title (set size)
                     1.2, # Size of values above set size
                     1, # Size of text for group
                     2), # Size of numbers on bars
      keep.order = TRUE)
dev.off()

library(UpSetR)
library(tidyverse)
library(dplyr)
library(magrittr)
library(tidyr)
library(readr)
library(stringr)

#### FUNCTIONS ####


filter_data <- function(merged_data,significant_value){
  merged_data$Location_Merge <- paste(merged_data$Chr,merged_data$BP,sep="_")
  merged_data$Population_Merge <- paste(merged_data$Pop.1.ID,merged_data$Pop.2.ID,sep="_")
  
  merged_data$SIGNIFICANT <- 0
  # change Fst value
  merged_data$SIGNIFICANT[merged_data$AMOVA.Fst>significant_value] <- 1
  
  # SUBSET DATA 
  
  one_data_subset <- data.frame(Location = merged_data$Location_Merge,
                                FST = merged_data$SIGNIFICANT,
                                Populations = merged_data$Population_Merge)
  
  one_data_subset_duplicates <- distinct(one_data_subset, Location, Populations, .keep_all= TRUE)
  
  
  one_data_subset_duplicates_wide <- one_data_subset_duplicates %>%
    pivot_wider(names_from = Populations, values_from = FST)
  
  # remove rows with NA values
  one_data_subset_duplicates_wide_final <- one_data_subset_duplicates_wide[complete.cases(one_data_subset_duplicates_wide), ]
  
  one_data_subset_duplicates_wide_final <- as.data.frame(one_data_subset_duplicates_wide_final)
  
  one_data_subset_duplicates_wide_final[-1] <- sapply(one_data_subset_duplicates_wide_final[-1], as.integer)
  
  one_data_subset_duplicates_wide_final$Location <- factor(one_data_subset_duplicates_wide_final$Location)

  return(one_data_subset_duplicates_wide_final)
  
}


#### LOAD DATA ####

# SOG: 

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


# SCO:

# Sco vs Mývatn 
#SCO_MYV = read.delim('stacks_populations_result/iceland.recode.p.fst_MYV-SCO.tsv')
# Sco vs Svínavatn 
#SCO_SVI = read.delim('stacks_populations_result/iceland.recode.p.fst_SCO-SVI.tsv')
# Sco vs Hestvatn 
#SCO_HES = read.delim('stacks_populations_result/iceland.recode.p.fst_HES-SCO.tsv')
# Sco vs Öxara
#SCO_OXA = read.delim('stacks_populations_result/iceland.recode.p.fst_OXA-SCO.tsv')
# Sco vs STF 
#SCO_STF = read.delim('stacks_populations_result/iceland.recode.p.fst_SCO-STF.tsv')



# SOG ONLY

merge_all_SOG <- rbind(SOG_MYV,SOG_SVI,SOG_HES,SOG_OXA,SOG_STF)


# Fst 0.5
## DEFAULT MYND
# SOG comparison folder
merge_all_SOG_filtered_05 <-  filter_data(merge_all_SOG,0.5)

pdf('../comparisons_SOG/upset_SOG_figure_fst05.pdf',width = 15, height=10, onefile=FALSE)
upset(merge_all_SOG_filtered_05,
      #sets =names(color),
      #sets.bar.color = color,
      order.by = "freq",
      nintersects = 40,
      set_size.show = TRUE,
      set_size.numbers_size = 15,
      set_size.scale_max = 280,
      text.scale = c(2.5, # size of y-axis title
                     1.5, # size of values on y-axis
                     2, # size of x-axis title (set size)
                     1.2, # Size of values above set size
                     2, # Size of text for group
                     2), # Size of numbers on bars
      keep.order = TRUE)
dev.off()





# SCO ONLY

merge_all_SCO <- rbind(SCO_HES,SCO_MYV,SCO_OXA,SCO_STF,SCO_SVI)


# Fst 0.5
## DEFAULT MYND
merge_all_SCO_filtered_05 <-  filter_data(merge_all_SCO,0.5)

png('plots/230215_new/upset_SCO_figure_fst05.png',width = 1000, height = 600,units = "px")
upset(merge_all_SCO_filtered_05,
      #sets =names(color),
      #sets.bar.color = color,
      order.by = "freq",
      nintersects = 40,
      set_size.show = TRUE,
      set_size.numbers_size = 15,
      set_size.scale_max = 280,
      text.scale = c(2.5, # size of y-axis title
                     1.5, # size of values on y-axis
                     2, # size of x-axis title (set size)
                     1.2, # Size of values above set size
                     1, # Size of text for group
                     2), # Size of numbers on bars
      keep.order = TRUE)
dev.off()


# ALL

merge_all <- rbind(SOG_MYV,SOG_SVI,SOG_HES,SOG_OXA,SOG_STF,SCO_MYV,SCO_SVI,SCO_HES,SCO_OXA,SCO_STF)


# Fst 0.5
merge_all_filtered_05 <-  filter_data(merge_all,0.5)

png('plots/230215_new/upset_ALL_figure_fst05.png',width = 2000, height = 600,units = "px")
upset(merge_all_filtered_05,
      #sets =names(color),
      #sets.bar.color = color,
      order.by = "freq",
      nintersects = 100,
      nsets = 10,
      set_size.show = TRUE,
      set_size.numbers_size = 10,
      set_size.scale_max = 220,
      text.scale = c(2.5, # size of y-axis title
                     1.5, # size of values on y-axis
                     2, # size of x-axis title (set size)
                     1.2, # Size of values above set size
                     1, # Size of text for group
                     2), # Size of numbers on bars
      keep.order = TRUE)
dev.off()



myList <- list()


myList[[1]] <- upset(merge_all_filtered_05,
              #sets =names(color),
              #sets.bar.color = color,
              order.by = "freq",
              nintersects = 100,
              nsets = 10,
              set_size.show = TRUE,
              set_size.numbers_size = 10,
              set_size.scale_max = 220,
              text.scale = c(2.5, # size of y-axis title
                             1.5, # size of values on y-axis
                             2, # size of x-axis title (set size)
                             1.2, # Size of values above set size
                             1, # Size of text for group
                             2), # Size of numbers on bars
              keep.order = TRUE)


pdf("plots/230215_new/test.pdf",16,12)
myList
dev.off()

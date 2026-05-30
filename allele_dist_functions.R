library(vcfR)
library(sqldf)
library(dplyr)
library(tidyr)
library(ggpubr)
library(UpSetR)




#### LOAD VCF ####

# read the vcf file 
vcfData <- read.vcfR("../marcos_brown_trout_pairwise_comparisons/iceland.recode.vcf", verbose = FALSE )

# convert to tidy
vcfData_tidy <- vcfR2tidy(vcfData)

# more like the meta-data, information about each SNP 
vcfData_tidy_fix <- vcfData_tidy$fix

# all the vcf data - one line for each sample on each SNP
vcfData_tidt_gt <- vcfData_tidy$gt

# match CHROM to gt data
# needed the chromosome name to match - used ChromKey for that
vcfData_tidt_gt$CHROM <- vcfData_tidy_fix$CHROM[match(vcfData_tidt_gt$ChromKey,vcfData_tidy_fix$ChromKey)]

# change dots to - this was more to have the data the same (some of my data had dots where some had underscores)
# vcfData_tidt_gt$CHROM <- str_replace_all(vcfData_tidt_gt$CHROM, c("\\." = "_"))

# create new id that is a compination of the chromosome and location of the SNP
vcfData_tidt_gt$ID <- paste(vcfData_tidt_gt$CHROM,vcfData_tidt_gt$POS,sep = "_")

# subset the table only with necessary information
# Here I subtracted the base genotypes, you could also ectract the 0/1 values, then just use the column gt_GT
gt_ID_table <- data.frame(IND=vcfData_tidt_gt$Indiv,ID=vcfData_tidt_gt$ID,ALLELE=vcfData_tidt_gt$gt_GT_alleles)


# extract only the populations we are using: 
# MYV (Myv), SVI (Svi), HES (Hes), OXA (Oxa), STF (Stf), SOG (Sog), SCO (Slp)

gt_ID_table$short_pop <- substr(gt_ID_table$IND, 0, 3)


gt_ID_table$Population <- "NO"
gt_ID_table$Population[gt_ID_table$short_pop=="Myv"] <- "MYV"
gt_ID_table$Population[gt_ID_table$short_pop=="Svi"] <- "SVI"
gt_ID_table$Population[gt_ID_table$short_pop=="Hes"] <- "HES"
gt_ID_table$Population[gt_ID_table$short_pop=="Oxa"] <- "OXA"
gt_ID_table$Population[gt_ID_table$short_pop=="Stf"] <- "STF"
gt_ID_table$Population[gt_ID_table$short_pop=="Sog"] <- "SOG"
gt_ID_table$Population[gt_ID_table$short_pop=="Slp"] <- "SLP"

df_subpopulations <- gt_ID_table[gt_ID_table$Population!="NO",]

df_all_filtered <- df_subpopulations[df_subpopulations$ALLELE!=".",]



#### FUNCTION TO GENERATE FIGURES ####

# allele distribution figure

allele_figure <- function(SNP,allele_df){
  # extract only rows with this SNP-ID
  allele_df_SNP <- allele_df[allele_df$ID %in% SNP, ]
  
  allele_df_SNP_subset <- data.frame(IND=allele_df_SNP$IND,
                                     POPULATION=allele_df_SNP$Population,
                                     ID=allele_df_SNP$ID,
                                     ALLELE=allele_df_SNP$ALLELE)
  # remove duplicated rows if any
  allele_df_SNP_subset_duplicates <- distinct(allele_df_SNP_subset, IND, ID, .keep_all= TRUE)
  
  # count allele-instances
  allele_total <- sqldf('SELECT Population, ALLELE, COUNT(ALLELE)
                         FROM allele_df_SNP_subset_duplicates
                         GROUP BY Population, ALLELE')
  
  colnames(allele_total) <- c("Population","Allele","count")
  
  allele_total_final <- allele_total %>% group_by(Population) %>% mutate(perAllele = count/sum(count)*100)
  
  allele_total_final$Population <- as.factor(allele_total_final$Population)
  allele_total_final$Population <- factor(allele_total_final$Population, levels = c("SOG","SLP", "HES","SVI","MYV","STF","OXA"))
  
  allele_distribution_figure <- allele_total_final %>%
    ggplot(aes(x=Population, 
               y=perAllele, 
               fill=Allele)) +
    geom_bar(stat="identity",colour="black") +
    geom_text(aes(label = paste(format(round(perAllele,1), nsmall = 1), "% (",count,")",sep="")), 
              color = "black", size = 5, position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = c("#E0FFFF","#87CEEB","#1E90FF")) +
    labs(x="Population",y="Percentage",fill="Alleles")+
    ggtitle(SNP)+
    theme_bw() +
    theme(plot.title = element_text(vjust = 0.5, size=22, face="bold"),
          axis.title.x = element_text(size=20, face="bold"),
          axis.title.y = element_text(size=20, face="bold"),
          legend.title = element_text(size=15, face="bold"),
          legend.text = element_text(size  = 15),
          axis.text.x = element_text(size  = 20),
          axis.text.y = element_text(size  = 20))
  
  return(allele_distribution_figure)
  
}

# filter data to generate upset-plot
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


#### LOAD FST ESTIMATED DATA ####


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




# merge all SNP data
merge_all <- rbind(SOG_MYV,SOG_SVI,SOG_HES,SOG_OXA,SOG_STF)

merge_all_05 <-  filter_data(merge_all,0.5)




#### FIGURES FOR PATTERNS 3 OR MORE ####

# default upset graph
upset_graph <- upset(merge_all_05,
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
                                    1), # Size of numbers on bars
                     keep.order = TRUE)

myList <- list()
i <- 1
myList[[i]] <- upset_graph


# first
upset_3_6 <- merge_all_05$Location[merge_all_05$MYV_SOG==1 &  
                                     merge_all_05$OXA_SOG==1 &
                                     merge_all_05$SOG_SVI==0 &  
                                     merge_all_05$HES_SOG==0 &
                                     merge_all_05$SOG_STF==1]
# NC_042958.1_32852905 -
# NC_042959.1_36387368 - 
# NC_042964.1_5260841 -
# NC_042970.1_51130646 -
# NC_042974.1_39567488 -
# NC_042992.1_2352301 -


SNP_list <- upset_3_6
i <- i+1

for (SNP in SNP_list) {
  SNP_frequency_figure <- allele_figure(SNP,df_all_filtered)
  myList[[i]] <- SNP_frequency_figure
  i <- i+1
}


myList[[i]] <- upset_graph

# second
upset_3_3 <- merge_all_05$Location[merge_all_05$MYV_SOG==0 &  
                                     merge_all_05$OXA_SOG==1 &
                                     merge_all_05$SOG_SVI==1 &  
                                     merge_all_05$HES_SOG==0 &
                                     merge_all_05$SOG_STF==1]

# NC_042969.1_64385521
# NC_042983.1_22238557
# NC_042990.1_23677206

SNP_list <- upset_3_3
i <- i+1

for (SNP in SNP_list) {
  SNP_frequency_figure <- allele_figure(SNP,df_all_filtered)
  myList[[i]] <- SNP_frequency_figure
  i <- i+1
}


myList[[i]] <- upset_graph

# third
upset_2_5 <- merge_all_05$Location[merge_all_05$MYV_SOG==1 &  
                                     merge_all_05$OXA_SOG==1 &
                                     merge_all_05$SOG_SVI==1 &  
                                     merge_all_05$HES_SOG==1 &
                                     merge_all_05$SOG_STF==1]
# NC_042966.1_28586393
# NC_042981.1_25831499

SNP_list <- upset_2_5
i <- i+1

for (SNP in SNP_list) {
  SNP_frequency_figure <- allele_figure(SNP,df_all_filtered)
  myList[[i]] <- SNP_frequency_figure
  i <- i+1
}

myList[[i]] <- upset_graph



upset_131 <- merge_all_05$Location[merge_all_05$MYV_SOG==1 &  
                                     merge_all_05$OXA_SOG==0 &
                                     merge_all_05$SOG_SVI==0 &  
                                     merge_all_05$HES_SOG==1 &
                                     merge_all_05$SOG_STF==1]
# NC_042961.1_36658707
SNP_list <- upset_131
i <- i+1

for (SNP in SNP_list) {
  SNP_frequency_figure <- allele_figure(SNP,df_all_filtered)
  myList[[i]] <- SNP_frequency_figure
  i <- i+1
}


myList[[i]] <- upset_graph

upset_132 <- merge_all_05$Location[merge_all_05$MYV_SOG==1 &  
                                     merge_all_05$OXA_SOG==1 &
                                     merge_all_05$SOG_SVI==1 &  
                                     merge_all_05$HES_SOG==0 &
                                     merge_all_05$SOG_STF==0]
# NC_042982.1_20624823

SNP_list <- upset_132
i <- i+1

for (SNP in SNP_list) {
  SNP_frequency_figure <- allele_figure(SNP,df_all_filtered)
  myList[[i]] <- SNP_frequency_figure
  i <- i+1
}


myList[[i]] <- upset_graph


upset_133 <- merge_all_05$Location[merge_all_05$MYV_SOG==0 &  
                                     merge_all_05$OXA_SOG==1 &
                                     merge_all_05$SOG_SVI==0 &  
                                     merge_all_05$HES_SOG==1 &
                                     merge_all_05$SOG_STF==1]
# NC_042987.1_12358347
SNP_list <- upset_133
i <- i+1

for (SNP in SNP_list) {
  SNP_frequency_figure <- allele_figure(SNP,df_all_filtered)
  myList[[i]] <- SNP_frequency_figure
  i <- i+1
}


myList[[i]] <- upset_graph

upset_1_4 <- merge_all_05$Location[merge_all_05$MYV_SOG==1 &  
                                     merge_all_05$OXA_SOG==1 &
                                     merge_all_05$SOG_SVI==0 &  
                                     merge_all_05$HES_SOG==1 &
                                     merge_all_05$SOG_STF==1]
# NC_042963.1_22010091
SNP_list <- upset_1_4
i <- i+1

for (SNP in SNP_list) {
  SNP_frequency_figure <- allele_figure(SNP,df_all_filtered)
  myList[[i]] <- SNP_frequency_figure
  i <- i+1
}

pdf("../comparisons_SOG/patterns_upset_3ormore.pdf",15,10)
myList
dev.off()





#### FIGURE FOR ARTICLE ####


allele_figure_multiple <- function(SNP,allele_df,position_text){
  # extract only rows with this SNP-ID
  allele_df_SNP <- allele_df[allele_df$ID %in% SNP, ]
  
  allele_df_SNP_subset <- data.frame(IND=allele_df_SNP$IND,
                                     POPULATION=allele_df_SNP$Population,
                                     ID=allele_df_SNP$ID,
                                     ALLELE=allele_df_SNP$ALLELE)
  # remove duplicated rows if any
  allele_df_SNP_subset_duplicates <- distinct(allele_df_SNP_subset, IND, ID, .keep_all= TRUE)
  
  # count allele-instances
  allele_total <- sqldf('SELECT Population, ALLELE, COUNT(ALLELE)
                         FROM allele_df_SNP_subset_duplicates
                         GROUP BY Population, ALLELE')
  
  colnames(allele_total) <- c("Population","Allele","count")
  
  allele_total_final <- allele_total %>% group_by(Population) %>% mutate(perAllele = count/sum(count)*100)
  
  allele_total_final$Population <- as.factor(allele_total_final$Population)
  allele_total_final$Population <- factor(allele_total_final$Population, levels = c("SOG","SLP", "HES","SVI","MYV","STF","OXA"))
  
  f <- allele_total_final %>%
    ggplot(aes(x=Population, 
               y=perAllele, 
               fill=Allele)) +
    geom_bar(stat="identity",colour="black") +
    geom_text(aes(label = paste("(",count,")",sep="")), 
              color = "black", size = 2, position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = c("#E0FFFF","#87CEEB","#1E90FF")) +
    labs(x="Population",y="Percentage",fill="Alleles")+
    ggtitle(SNP)+
    theme_bw() +
    theme(plot.title = element_text(vjust = 0.5, size=10, face="bold"),
          axis.title.x = element_text(size=10, face="bold"),
          axis.title.y = element_text(size=10, face="bold"),
          legend.title = element_text(size=10, face="bold"),
          legend.text = element_text(size  = 10),
          axis.text.x = element_text(size  = 10),
          axis.text.y = element_text(size  = 10))
  
  if (position_text=="00") {
    ## figure number == 00: no y-axis title, no x-axis title
    f <- f + labs(x= "", y= "")
  }
  else if (position_text=="10") {
    ## figure number == 10: y-axis title, no x-axis title 
    f <- f + labs(x= "", y= "Percentage")
  }
  else if (position_text=="01") {
    ## figure number == 01: no y-axis title, x-axis title
    f <- f + labs(x= "Population", y= "")
  }
  else if (position_text=="11") {
    ## figure number == 11: y-axis title, x-axis title
    f <- f + labs(x= "Population", y= "Percentage")
  }
  
  
  return(f)
  
}


# Choose SNPs

# Pattern: SNPs in all 5 comparisons, then 4, then 3, then 2, and then only in 1 comparison

# 5 comparisons - use both:
five_comparisons <- merge_all_05$Location[merge_all_05$MYV_SOG==1 &  
                                            merge_all_05$OXA_SOG==1 &
                                            merge_all_05$SOG_SVI==1 &  
                                            merge_all_05$HES_SOG==1 &
                                            merge_all_05$SOG_STF==1]

# 4 comparisons - only use one:
four_comparisons <- merge_all_05$Location[merge_all_05$MYV_SOG==1 &  
                                            merge_all_05$OXA_SOG==1 &
                                            merge_all_05$SOG_SVI==0 &  
                                            merge_all_05$HES_SOG==1 &
                                            merge_all_05$SOG_STF==1]

# 3 comparisons - use three, one (first) from each case:

three_comparisons_1 <- merge_all_05$Location[merge_all_05$MYV_SOG==1 &  
                                               merge_all_05$OXA_SOG==1 &
                                               merge_all_05$SOG_SVI==0 &  
                                               merge_all_05$HES_SOG==0 &
                                               merge_all_05$SOG_STF==1]

three_comparisons_2 <- merge_all_05$Location[merge_all_05$MYV_SOG==0 &  
                                               merge_all_05$OXA_SOG==1 &
                                               merge_all_05$SOG_SVI==1 &  
                                               merge_all_05$HES_SOG==0 &
                                               merge_all_05$SOG_STF==1]

three_comparisons_3 <- merge_all_05$Location[merge_all_05$MYV_SOG==1 &  
                                               merge_all_05$OXA_SOG==0 &
                                               merge_all_05$SOG_SVI==0 &  
                                               merge_all_05$HES_SOG==1 &
                                               merge_all_05$SOG_STF==1]


# 2 comparisons - use two, one (first) from both case:


two_comparisons_1 <- merge_all_05$Location[merge_all_05$MYV_SOG==1 &  
                                             merge_all_05$OXA_SOG==1 &
                                             merge_all_05$SOG_SVI==0 &  
                                             merge_all_05$HES_SOG==0 &
                                             merge_all_05$SOG_STF==0]

two_comparisons_2 <- merge_all_05$Location[merge_all_05$MYV_SOG==0 &  
                                               merge_all_05$OXA_SOG==1 &
                                               merge_all_05$SOG_SVI==0 &  
                                               merge_all_05$HES_SOG==0 &
                                               merge_all_05$SOG_STF==1]


# 2 comparisons - use two, one (first) from both case:


one_comparisons_1 <- merge_all_05$Location[merge_all_05$MYV_SOG==0 &  
                                             merge_all_05$OXA_SOG==1 &
                                             merge_all_05$SOG_SVI==0 &  
                                             merge_all_05$HES_SOG==0 &
                                             merge_all_05$SOG_STF==0]

one_comparisons_2 <- merge_all_05$Location[merge_all_05$MYV_SOG==0 &  
                                             merge_all_05$OXA_SOG==0 &
                                             merge_all_05$SOG_SVI==0 &  
                                             merge_all_05$HES_SOG==0 &
                                             merge_all_05$SOG_STF==1]





# SOG comparisons
allele_distributions <- ggarrange(allele_figure_multiple(five_comparisons[1],df_all_filtered,"10"),allele_figure_multiple(five_comparisons[2],df_all_filtered,"00"),
                                   allele_figure_multiple(four_comparisons[1],df_all_filtered,"10"),allele_figure_multiple(three_comparisons_1[1],df_all_filtered,"00"),
                                   allele_figure_multiple(three_comparisons_2[1],df_all_filtered,"10"),allele_figure_multiple(three_comparisons_3[1],df_all_filtered,"00"),
                                   allele_figure_multiple(two_comparisons_1[1],df_all_filtered,"10"),allele_figure_multiple(two_comparisons_2[1],df_all_filtered,"00"),
                                   allele_figure_multiple(one_comparisons_1[1],df_all_filtered,"11"),allele_figure_multiple(one_comparisons_2[1],df_all_filtered,"01"),
                                   ncol=2,nrow=5, common.legend = FALSE,labels="AUTO")

# SOG folder
ggsave("../comparisons_SOG/allele_distributions_allexapmles.pdf", width = 15, height = 10)

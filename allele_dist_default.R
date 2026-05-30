library(vcfR)
library(sqldf)
library(dplyr)
library(tidyr)




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

# create LG+BP ID for each SNP
merge_all$Location_Merge <- paste(merge_all$Chr,merge_all$BP,sep="_")

# merge the two compared populations
merge_all$Population_Merge <- paste(merge_all$Pop.1.ID,merge_all$Pop.2.ID,sep="_")



merge_all$sign05 <- 0

# change Fst value
merge_all$sign05[merge_all$AMOVA.Fst>0.5] <- 1

# subset dataframe
significant_SNPs_fst05_subset <- data.frame(Location = merge_all$Location_Merge,
                                            FST = merge_all$sign05,
                                            Populations = merge_all$Population_Merge)


# remove duplication
significant_SNPs_fst05_subset_dup <- distinct(significant_SNPs_fst05_subset, Location, Populations, .keep_all= TRUE)

# change format
significant_SNPs_fst05_subset_dup_wide <- significant_SNPs_fst05_subset_dup %>%
  pivot_wider(names_from = Populations, values_from = FST)

# remove rows with missing values (NA)
significant_SNPs_fst05_subset_dup_wide_na <- significant_SNPs_fst05_subset_dup_wide[complete.cases(significant_SNPs_fst05_subset_dup_wide), ]

# only keep rows with at least one case where a SNP was sifnificant
significant_SNPs_fst05_subset_dup_wide_na_final <- significant_SNPs_fst05_subset_dup_wide_na[rowSums(significant_SNPs_fst05_subset_dup_wide_na[c(2:6)])>0,]


# all figures



SNP_list <- significant_SNPs_fst05_subset_dup_wide_na_final$Location

SNP <- SNP_list[1]

df_subpopulations_filtered_oneSNP <- df_all_filtered[df_all_filtered$ID %in% SNP, ]
  
df_INDonly <- data.frame(IND=df_subpopulations_filtered_oneSNP$IND,
                         POPULATION=df_subpopulations_filtered_oneSNP$Population,
                         ID=df_subpopulations_filtered_oneSNP$ID,
                         ALLELE=df_subpopulations_filtered_oneSNP$ALLELE)
  
df_INDonly_duplicates <- distinct(df_INDonly, IND, ID, .keep_all= TRUE)
  
total <- sqldf('SELECT Population, ALLELE, COUNT(ALLELE)
                FROM df_INDonly_duplicates
                GROUP BY Population, ALLELE')
  
colnames(total) <- c("Population","Allele","count")
  
total_final <- total %>% group_by(Population) %>% mutate(perAllele = count/sum(count)*100)
  
total_final$Population <- as.factor(total_final$Population)
total_final$Population <- factor(total_final$Population, levels = c("SOG","SLP", "HES","SVI","MYV","STF","OXA"))
  
freq_allele_label <- total_final %>%
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
  
freq_allele_label

  
ggsave("../comparisons_SOG/allele_plot_onetest.pdf", width = 10, height = 10)

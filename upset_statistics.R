#### LOAD DATA ####
library(dplyr)
library(tidyr)

SOG_MYV = read.delim('../stacks_populations_result/iceland.recode.p.fst_MYV-SOG.tsv')
SOG_SVI = read.delim('../stacks_populations_result/iceland.recode.p.fst_SOG-SVI.tsv')
SOG_HES = read.delim('../stacks_populations_result/iceland.recode.p.fst_HES-SOG.tsv')
SOG_OXA = read.delim('../stacks_populations_result/iceland.recode.p.fst_OXA-SOG.tsv')
SOG_STF = read.delim('../stacks_populations_result/iceland.recode.p.fst_SOG-STF.tsv')


merge_all_SOG <- rbind(SOG_MYV,SOG_SVI,SOG_HES,SOG_OXA,SOG_STF)


merge_all_SOG$Location_Merge <- paste(merge_all_SOG$Chr,merge_all_SOG$BP,sep="_")
merge_all_SOG$Population_Merge <- paste(merge_all_SOG$Pop.1.ID,merge_all_SOG$Pop.2.ID,sep="_")

merge_all_SOG$SIGNIFICANT <- 0
# change Fst value
merge_all_SOG$SIGNIFICANT[merge_all_SOG$AMOVA.Fst>0.5] <- 1


# SUBSET DATA 

one_data_subset <- data.frame(Location = merge_all_SOG$Location_Merge,
                              FST = merge_all_SOG$SIGNIFICANT,
                              Populations = merge_all_SOG$Population_Merge)

one_data_subset_duplicates <- distinct(one_data_subset, Location, Populations, .keep_all= TRUE)


one_data_subset_duplicates_wide <- one_data_subset_duplicates %>%
  pivot_wider(names_from = Populations, values_from = FST)

# remove rows with NA values
one_data_subset_duplicates_wide_final <- one_data_subset_duplicates_wide[complete.cases(one_data_subset_duplicates_wide), ]



#adapt data to correct dataformat for the other package

# choose the correct vectors

SOG_MYV_sig_vec <- one_data_subset_duplicates_wide_final$Location[one_data_subset_duplicates_wide_final$MYV_SOG==1]
SOG_SVI_sig_vec <- one_data_subset_duplicates_wide_final$Location[one_data_subset_duplicates_wide_final$SOG_SVI==1]
SOG_HES_sig_vec <- one_data_subset_duplicates_wide_final$Location[one_data_subset_duplicates_wide_final$HES_SOG==1]
SOG_OXA_sig_vec <- one_data_subset_duplicates_wide_final$Location[one_data_subset_duplicates_wide_final$OXA_SOG==1]
SOG_STF_sig_vec <- one_data_subset_duplicates_wide_final$Location[one_data_subset_duplicates_wide_final$SOG_STF==1]

# create a list

SOG_comparisons_list <- list()

SOG_comparisons_list[["MYV_SOG"]] <- SOG_MYV_sig_vec
SOG_comparisons_list[["SOG_SVI"]] <- SOG_SVI_sig_vec
SOG_comparisons_list[["HES_SOG"]] <- SOG_HES_sig_vec
SOG_comparisons_list[["OXA_SOG"]] <- SOG_OXA_sig_vec
SOG_comparisons_list[["SOG_STF"]] <- SOG_STF_sig_vec


# number of SNPs

total <- length(unique(one_data_subset_duplicates_wide_final$Location))



#### RUN PACKAGE ####

library("SuperExactTest")



str(SOG_comparisons_list)

(length.gene.sets=sapply(SOG_comparisons_list,length))

(num.expcted.overlap=total*do.call(prod,as.list(length.gene.sets/total)))

(p=sapply(0:101,function(i) dpsets(i, length.gene.sets, n=total)))

common.genes=intersect(SOG_comparisons_list[[1]], SOG_comparisons_list[[2]], SOG_comparisons_list[[3]],
                       SOG_comparisons_list[[4]], SOG_comparisons_list[[5]])
(num.observed.overlap=length(common.genes))

(FE=num.observed.overlap/num.expcted.overlap)

dpsets(num.observed.overlap, length.gene.sets, n=total)

cpsets(num.observed.overlap-1, length.gene.sets, n=total, lower.tail=FALSE)

fit=MSET(SOG_comparisons_list, n=total, lower.tail=FALSE)
fit$FE
fit$p.value



res=supertest(SOG_comparisons_list, n=total)

plot(res, sort.by="size", margin=c(2,2,2,2), color.scale.pos=c(0.85,1), legend.pos=c(0.9,0.15))




plot(res, Layout="landscape", degree=2:4, sort.by="size", margin=c(0.5,5,1,2))


pdf('../comparisons_SOG/upset_SOG_figure_fst05_statistics.pdf',width = 15, height=10, onefile=FALSE)
plot(res, Layout="landscape", degree=2:4, sort.by="size", margin=c(0.5,5,1,2))
dev.off()

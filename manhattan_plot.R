#### PACKAGES ####

library(HistogramTools)
library(ggplot2)
library(qqman)
library(magrittr)
library(dplyr)
library(ggrepel)
library(readr)



#### LOAD REFERENCE GENOME ####

salmo_genome <- read.csv(file="marcos_brown_trout_pairwise_comparisons/Salmo_trutta_chr_lengths.csv", stringsAsFactors=FALSE, sep="")

# tmp <- Stchr$CumLen
# names(tmp) <- Stchr$Chr_ID

salmo_genome <- salmo_genome %>% group_by(LG) %>% mutate(Middle=min(CumLen)+(LGlength)/2)

LG_color_df <- data.frame(LG_unique=salmo_genome$LG)
LG_color_df <- LG_color_df %>% distinct()
LG_color_df$pattern <- rep(c(TRUE,FALSE), length.out=nrow(LG_color_df))

salmo_genome <- left_join(salmo_genome, LG_color_df, by=c('LG'='LG_unique'))



#### LOAD DATA ####


fst_data = read.delim('marcos_brown_trout_pairwise_comparisons/iceland.habitat.clean.recode.p.fst_river-lake.tsv')

# ok$CumPos <- ok$BP + tmp[ok$Chr]
fst_data_join <- left_join(fst_data, salmo_genome, by=c('Chr'='Chr_ID'))
fst_data_join <- as.data.frame(fst_data_join)
fst_data_join$realBP <- fst_data_join$BP + fst_data_join$CumLen


# tmp <- Stchr$LG
# names(tmp) <- Stchr$Chr_ID
# ok$Chr <- tmp[ok$Chr]


plot(fst_data_join$realBP,fst_data_join$AMOVA.Fst, ylim=c(0,1), pch=20, col=fst_data_join$LG)



#### SIGNIFICANT SNPS - HISTOGRAM ####

Histogram_SNPS <- hist(fst_data_join$AMOVA.Fst)

ApproxQuantile(Histogram_SNPS, 0.995)

# significant value
hist_values <- unname(ApproxQuantile(Histogram_SNPS, 0.995))

plot(fst_data_join$realBP, fst_data_join$AMOVA.Fst, pch=20, cex=1.5, 
     col=ifelse(fst_data_join$AMOVA.Fst>hist_values, "red", "black"), 
     main="River vs lake", xlab="Position in Brown trout Genome (BP)", ylab="Fst values")


fst_data_join$pattern[fst_data_join$AMOVA.Fst>hist_values] <- "Significant"




#### PLOT ####

# color palette

color_p <- c("#999999","#98cce5","red")
names(color_p) <- c("TRUE","FALSE","Significant")


ggplot(fst_data_join, aes(x=realBP, y=AMOVA.Fst)) +
  geom_point(aes(color=as.factor(pattern)), alpha=1, size=2.0) +
  scale_color_manual(values = color_p) +
  
  # custom X axis:
  scale_x_continuous(label = unique(fst_data_join$LG), breaks= unique(fst_data_join$Middle)) +
  scale_y_continuous(limits=c(0,1), expand = c(0.01, 0.01) ) +    # remove space between plot area and x axis
  
  #Add labels
  labs(x= "Chromosome", y= "Fst") +
  
  # Custom the theme:
  theme_bw() +
  theme( legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(vjust = 0.5, size=15, face="bold"),
    axis.text.x =
      element_text(size  = 8,
                   angle = 90,
                   hjust = 1,
                   vjust = 0.5),
    axis.text.y =
      element_text(size  = 8,
                   angle = NULL,
                   hjust = 1,
                   vjust = NULL))









#### DUPLICATED IDs ####

fst_data_join$Location_Merge <- paste(fst_data_join$Chr,fst_data_join$BP,sep="_")

fst_data_join$Location_Merge[duplicated(fst_data_join$Location_Merge)]

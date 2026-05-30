#### PACKAGES ####

library(HistogramTools)
library(ggplot2)
library(qqman)
library(magrittr)
library(dplyr)
library(ggrepel)
library(readr)
library(ggpubr)



#### LOAD REFERENCE GENOME ####

salmo_genome <- read.csv(file="../marcos_brown_trout_pairwise_comparisons/Salmo_trutta_chr_lengths.csv", stringsAsFactors=FALSE, sep="")

# calculate middle of chromosome
salmo_genome <- salmo_genome %>% group_by(LG) %>% mutate(Middle=min(CumLen)+(LGlength)/2)

# group chromosomes - used later for the shade
LG_color_df <- data.frame(LG_unique=salmo_genome$LG)
LG_color_df <- LG_color_df %>% distinct()
LG_color_df$pattern <- rep(c(TRUE,FALSE), length.out=nrow(LG_color_df))

# merge
salmo_genome <- left_join(salmo_genome, LG_color_df, by=c('LG'='LG_unique'))


#### LOAD DATA ####

# test data:

# Sog vs Mývatn 
# also SOG_MYV
SNP_data = read.delim('stacks_populations_result/iceland.recode.p.fst_MYV-SOG.tsv')
SCO_MYV = read.delim('stacks_populations_result/iceland.recode.p.fst_MYV-SCO.tsv')


fst_data_join <- left_join(SNP_data, salmo_genome, by=c('Chr'='Chr_ID'))
fst_data_join <- as.data.frame(fst_data_join)
fst_data_join$realBP <- fst_data_join$BP + fst_data_join$CumLen

# group significant SNPs - can be changed
# changeable significant value
significant_value <- 0.5
fst_data_join$pattern[fst_data_join$AMOVA.Fst>significant_value] <- "Significant"


color_p <- c("#999999","#98cce5","red")
names(color_p) <- c("TRUE","FALSE","Significant")


df <- fst_data_join
f <- ggplot(df, aes(x=realBP, y=AMOVA.Fst)) +
  geom_point(aes(color=as.factor(pattern)), alpha=1, size=1.5) +
  scale_color_manual(values = color_p) +
  
  # custom X axis:
  scale_x_continuous(label = unique(df$LG), breaks= unique(df$Middle)) +
  scale_y_continuous(limits=c(0,1), expand = c(0.01, 0.01) ) +    # remove space between plot area and x axis
  
  #Add labels
  labs(x= "", y= "") +
  ggtitle("title") + 
  
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
f

f <- f + labs(x= "chromosomes", y= "Fst") 

f


#### SCATTER ####

SOG_MYV$Location_Merge <- paste(SOG_MYV$Chr,SOG_MYV$BP,sep="_")
SOG_MYV$SOG_sig <- "NOT"
SOG_MYV$SOG_sig[SOG_MYV$AMOVA.Fst>0.5] <- "SIG"

SOG_sub <- data.frame(Location = SOG_MYV$Location_Merge,
                      SOG_FST = SOG_MYV$AMOVA.Fst,
                      SOG_SIG = SOG_MYV$SOG_sig)


SCO_MYV$Location_Merge <- paste(SCO_MYV$Chr,SCO_MYV$BP,sep="_")
SCO_MYV$SCO_sig <- "NOT"
SCO_MYV$SCO_sig[SCO_MYV$AMOVA.Fst>0.5] <- "SIG"

SCO_sub <- data.frame(Location = SCO_MYV$Location_Merge,
                      SCO_FST = SCO_MYV$AMOVA.Fst,
                      SCO_SIG = SCO_MYV$SCO_sig)


MYV_merge <- merge(SOG_sub, SCO_sub, by.x = "Location", by.y = "Location")


MYV_merge$pattern <- "NO"

MYV_merge$pattern[MYV_merge$SOG_SIG=="SIG" & MYV_merge$SCO_SIG=="SIG"] <- "Both"
MYV_merge$pattern[MYV_merge$SOG_SIG=="SIG" & MYV_merge$SCO_SIG=="NOT"] <- "SOGonly"
MYV_merge$pattern[MYV_merge$SOG_SIG=="NOT" & MYV_merge$SCO_SIG=="SIG"] <- "SCOonly"

ggplot(MYV_merge, aes(x=SOG_FST, y=SCO_FST, color=pattern)) + 
  geom_point() +
  ggtitle("MYV") +
  theme_bw()



ggsave("plots/230206_new/scatter_MYV_fst05.png", width = 15, height = 10)

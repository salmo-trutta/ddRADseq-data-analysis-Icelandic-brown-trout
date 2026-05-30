
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


# function to analyse the data and make figures

analyse_data <- function(salmo_reference,SNP_data,comparison_string,significant_value,figure_number) {
  # join SNP and reference to calculate "real" location 
  fst_data_join <- left_join(SNP_data, salmo_reference, by=c('Chr'='Chr_ID'))
  fst_data_join <- as.data.frame(fst_data_join)
  fst_data_join$realBP <- fst_data_join$BP + fst_data_join$CumLen
  
  # group significant SNPs - can be changed
  # changeable significant value
  fst_data_join$pattern[fst_data_join$AMOVA.Fst>significant_value] <- "Significant"
  
  manhattan_figure <- create_figure(fst_data_join,comparison_string,figure_number)
  
  return(manhattan_figure)
  
}
## figure number == 00: no y-axis title, no x-axis title
## figure number == 10: y-axis title, no x-axis title 
## figure number == 01: no y-axis title, x-axis title
## figure number == 11: y-axis title, x-axis title

create_figure <- function(df,string_title,f_number_string) {
  color_p <- c("#c7c7c7","#595959","red")
  names(color_p) <- c("TRUE","FALSE","Significant")
  
  f <- ggplot(df, aes(x=realBP, y=AMOVA.Fst)) +
    geom_point(aes(color=as.factor(pattern)), alpha=1, size=1.5) +
    scale_color_manual(values = color_p) +
    
    # custom X axis:
    scale_x_continuous(label = unique(df$LG), breaks= unique(df$Middle)) +
    scale_y_continuous(limits=c(0,1), expand = c(0.01, 0.01) ) +    # remove space between plot area and x axis
    
    #Add labels
    # labs(x= "Chromosome", y= "Fst") +
    ggtitle(string_title) + 
    
    # Custom the theme:
    theme_bw() +
    theme( legend.position="none",
           panel.border = element_blank(),
           panel.grid.major.x = element_blank(),
           panel.grid.minor.x = element_blank(),
           plot.title = element_text(vjust = 0.5, size=15, face="bold"),
           axis.title.x = element_text(size=17, face="bold"),
           axis.title.y = element_text(size=17, face="bold"),
           axis.text.x =
             element_text(size  = 8,
                          angle = 90,
                          hjust = 1,
                          vjust = 0.5),
           axis.text.y =
             element_text(size  = 10,
                          angle = NULL,
                          hjust = 1,
                          vjust = NULL))
  
  if (f_number_string=="00") {
    ## figure number == 00: no y-axis title, no x-axis title
    f <- f + labs(x= "", y= "")
  }
  else if (f_number_string=="10") {
    ## figure number == 10: y-axis title, no x-axis title 
    f <- f + labs(x= "", y= "Fst")
  }
  else if (f_number_string=="01") {
    ## figure number == 01: no y-axis title, x-axis title
    f <- f + labs(x= "Chromosome", y= "")
  }
  else if (f_number_string=="11") {
    ## figure number == 11: y-axis title, x-axis title
    f <- f + labs(x= "Chromosome", y= "Fst")
  }
  
  
  
  return(f)
  
}




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


# Slp vs Mývatn 
# SLP_MYV = read.delim('stacks_populations_result/iceland.recode.p.fst_MYV-SCO.tsv')
# Slp vs Svínavatn 
# SLP_SVI = read.delim('stacks_populations_result/iceland.recode.p.fst_SCO-SVI.tsv')
# Slp vs Hestvatn 
# SLP_HES = read.delim('stacks_populations_result/iceland.recode.p.fst_HES-SCO.tsv')
# Slp vs Öxara
# SLP_OXA = read.delim('stacks_populations_result/iceland.recode.p.fst_OXA-SCO.tsv')
# Slp vs STF 
# SLP_STF = read.delim('stacks_populations_result/iceland.recode.p.fst_SCO-STF.tsv')



#### FIGURES ####


# New figures from meeting 230201


# SOG comparisons
SOG_comparisons_fst05 <- ggarrange(analyse_data(salmo_genome,SOG_HES,"SOG vs. HES",0.5,"10"),
                                   analyse_data(salmo_genome,SOG_SVI,"SOG vs. SVI",0.5,"10"),
                                   analyse_data(salmo_genome,SOG_MYV,"SOG vs. MYV",0.5,"10"),
                                   analyse_data(salmo_genome,SOG_STF,"SOG vs. STF",0.5,"10"),
                                   analyse_data(salmo_genome,SOG_OXA,"SOG vs. OXA",0.5,"11"),
                                   ncol=1,nrow=5, common.legend = FALSE,labels="AUTO")

# SOG folder
ggsave("../comparisons_SOG/manhattan_pairwisecomparisons_SOG_fst05.pdf", width = 15, height = 10)



# ALL comparisons
ALL_comparisons_fst05 <- ggarrange(analyse_data(salmo_genome,SOG_HES,"SOG vs. HES",0.5,"10"),analyse_data(salmo_genome,SCO_HES,"SCO vs. HES",0.5,"00"),
                                   analyse_data(salmo_genome,SOG_SVI,"SOG vs. SVI",0.5,"10"),analyse_data(salmo_genome,SCO_SVI,"SCO vs. SVI",0.5,"00"),
                                   analyse_data(salmo_genome,SOG_MYV,"SOG vs. MYV",0.5,"10"),analyse_data(salmo_genome,SCO_MYV,"SCO vs. MYV",0.5,"00"),
                                   analyse_data(salmo_genome,SOG_STF,"SOG vs. STF",0.5,"10"),analyse_data(salmo_genome,SCO_STF,"SCO vs. STF",0.5,"00"),
                                   analyse_data(salmo_genome,SOG_OXA,"SOG vs. OXA",0.5,"11"),analyse_data(salmo_genome,SCO_OXA,"SCO vs. OXA",0.5,"01"),
                                   ncol=2,nrow=5, common.legend = FALSE)

#ggsave("plots/230206_new/manhattan_pairwisecomparisons_ALL_fst05.png", width = 20, height = 18)

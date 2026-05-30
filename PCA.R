library(vcfR)
library(adegenet)

vcf <- read.vcfR("iceland.recode.vcf")

x <- vcfR2genind(vcf)

# Dealing with NA values
y <- tab(x, NA.method="mean")

# Performing a PCA
y.pca = dudi.pca(y, scannf=FALSE, nf=20)

library(factoextra)
library(ggpubr)
library(ggplot2)


fviz_eig(y.pca)

fviz_eig(y.pca, addlabels = TRUE) + theme(text = element_text(family="Times New Roman")) + labs(title = NULL)

pop <- read.table(file="iceland.tsv", header=F)
sampling_sites <- as.factor(pop$V2)

str(sampling_sites)

## with dim1 vs dim2
A <-fviz_pca_ind(y.pca,
             habillage = sampling_sites,
             axes = 1:2,
             geom = "point", 
             repel = FALSE,    
             addEllipses = F,
             label="none",
             pointshape=19,
             pointsize=3.5,
             palette=c("#ff4c00", "#ff5d19", "#009700", "#580058", "#6533cb", "#4ca64c", "#8b0000", "#125699", "#004c4c", "#c1adea",
                       "#004d4d", "#0483ff", "#1e90ff", "#00cccc", "#006767", "#720072", "#ff69b4", "#961919", "#7447d1", "#a500a5", 
                       "#ff6f32", "#ff814c", "#007f7f", "#003100", "#ff9366", "#be00be", "#666666", "#00b2b2", "#009a9a", "#7f0000", 
                       "#a23232", "#ffa57f", "#d800d8", "#007e00", "#006ad1", "#6f0000", "#b3b300", "#00b100", "#4b0082", "#329932")) +
  coord_cartesian(ylim=c(-60, 80)) + 
  scale_y_continuous(breaks = c(-60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80)) + 
  coord_cartesian(xlim=c(-50, 60)) + 
  scale_x_continuous(breaks = c(-50, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60)) + theme(legend.position="none") +
  theme(plot.title=element_blank()) +
  labs(title = NULL, x = "PC1 (12%)", y = "PC2 (8.8%)") # + theme(legend.title = element_text(size=11, face="bold"), legend.text = element_text(size=9)) + 
  #guides(col = guide_legend(nrow = 42)) + labs(color = "Sampling site", face="bold")

A

## with dim3 vs dim4
B <-fviz_pca_ind(y.pca,
                 habillage = sampling_sites,
                 axes = 3:4,
                 geom = "point", 
                 repel = FALSE,    
                 addEllipses = F,
                 label="none",
                 pointshape=19,
                 pointsize=3.5,
                 palette=c("#ff4c00", "#ff5d19", "#009700", "#580058", "#6533cb", "#4ca64c", "#8b0000", "#125699", "#004c4c", "#c1adea",
                           "#004d4d", "#0483ff", "#1e90ff", "#00cccc", "#006767", "#720072", "#ff69b4", "#961919", "#7447d1", "#a500a5", 
                           "#ff6f32", "#ff814c", "#007f7f", "#003100", "#ff9366", "#be00be", "#666666", "#00b2b2", "#009a9a", "#7f0000", 
                           "#a23232", "#ffa57f", "#d800d8", "#007e00", "#006ad1", "#6f0000", "#b3b300", "#00b100", "#4b0082", "#329932")) + 
  coord_cartesian(ylim=c(-60, 80)) + 
  scale_y_continuous(breaks = c(-60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80)) + 
  coord_cartesian(xlim=c(-50, 60)) + 
  scale_x_continuous(breaks = c(-50, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60)) + theme(legend.position="none") +
  theme(plot.title=element_blank()) +
  labs(title = NULL, x = "PC3 (4.5%)", y = "PC4 (3.1%)") #+ theme(legend.title = element_text(size=11, face="bold"), legend.text = element_text(size=9)) + 
  #guides(col = guide_legend(nrow = 42)) + labs(color = "Sampling site", face="bold")

B

## with dim5 vs dim6
C <-fviz_pca_ind(y.pca,
                 habillage = sampling_sites,
                 axes = 5:6,
                 geom = "point", 
                 repel = FALSE,    
                 addEllipses = F,
                 label="none",
                 pointshape=19,
                 pointsize=3.5,
                 palette=c("#ff4c00", "#ff5d19", "#009700", "#580058", "#6533cb", "#4ca64c", "#8b0000", "#125699", "#004c4c", "#c1adea",
                           "#004d4d", "#0483ff", "#1e90ff", "#00cccc", "#006767", "#720072", "#ff69b4", "#961919", "#7447d1", "#a500a5", 
                           "#ff6f32", "#ff814c", "#007f7f", "#003100", "#ff9366", "#be00be", "#666666", "#00b2b2", "#009a9a", "#7f0000", 
                           "#a23232", "#ffa57f", "#d800d8", "#007e00", "#006ad1", "#6f0000", "#b3b300", "#00b100", "#4b0082", "#329932")) + 
  coord_cartesian(ylim=c(-60, 80)) + 
  scale_y_continuous(breaks = c(-60, -50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80)) + 
  coord_cartesian(xlim=c(-50, 60)) + 
  scale_x_continuous(breaks = c(-50, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60)) + theme(legend.position="none") +
  theme(plot.title=element_blank())+
  labs(title = NULL, x = "PC5 (2.5%)", y = "PC6 (2.1%)") #+ theme(legend.title = element_text(size=11, face="bold"), legend.text = element_text(size=9)) + 
  #guides(col = guide_legend(nrow = 42)) + labs(color = "Sampling site", face="bold")

C

library(ggplot2)
library(ggpubr)
library(dplyr)

ggarrange(A, B, C, labels = c("A", "B", "C"), font.label = list(size = 14), ncol = 1, nrow=3, common.legend = T, legend="none")

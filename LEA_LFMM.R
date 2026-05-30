##USING THE PACKAGE LEA##
library(LEA)
library(vcfR)

vcf2geno("iceland.temp.recode.vcf", "temp.geno")

pc = pca("temp.geno", scale = F)

tw = tracy.widom(pc)
tw$pvalues[1:5]
plot(tw$percentage)

project = snmf("temp.geno",
               K = 7:9,
               entropy = TRUE,
               CPU = 4,
               alpha = 100,
               seed = 125,
               repetitions = 2500,
               project = "new")

plot(project, col = "blue", pch = 19, cex = 1.2)


project<-load.snmfProject("temp.snmfProject")

best = which.min(cross.entropy(project, K = 7))

impute(project, "temp.geno", method = 'mode', K = 7, run = best)

obj.lfmm <- lfmm("temp.lfmm_imputed.lfmm", "temp.env", K = 7, rep = 100, project = "new", all = F, seed = 123456, 
                 missing.data = F, CPU = 6,  iterations = 20000, burnin = 10000)


obj.lfmm<-load.lfmmProject("temp.lfmmProject")
#Record z-scores from the 100 runs in the zs matrix
zs = z.scores(obj.lfmm, K = 7)

#Combine z-scores using the median
zs.median = apply(zs, MARGIN = 1, median)

#Compute the GIF
lambda = median(zs.median^2)/qchisq(0.5, df = 1)
lambda

# lambda = 0.6724973

# compute adjusted p-values from the combined z-scores
adj.p.values = pchisq(zs.median^2/lambda, df = 1, lower = FALSE)


#histogram of p-values
hist(adj.p.values, col = "red")

adj.p.values

# compute adjusted p-values from the combined z-scores
adj.p.values = pchisq(zs.median^2/0.45, df = 1, lower = FALSE)

#histogram of p-values
hist(adj.p.values, col = "green", main=NULL, xlab = "Adjusted p values")

plot(-log10(adj.p.values), pch = 19, col = "red", cex = .7)

## FDR control: Benjamini-Hochberg at level q
## L = number of loci
L = 5614
#fdr level q
q = 0.0001
w = which(sort(adj.p.values) < q * (1:L)/L)
candidates.bh = order(adj.p.values)[w]

write.csv(candidates.bh, "candidates")

## FDR control: Storey's q-values 
library(qvalue)
plot(qvalue(adj.p.values))
candidates.qv = which(qvalue(adj.p.values, fdr = 0.0001)$signif)
candidates.qv

plot(-log10(adj.p.values), main="Manhattan plot", xlab = "Locus", cex = .7, col = "grey")

points(candidates.qv, -log10(adj.p.values)[candidates.qv], pch = 19, cex = .7, col = "red")

### making a Manhattan plot

Stchr <- read.csv(file="Salmo_trutta_chr_lengths.csv", stringsAsFactors=FALSE, sep="")
tmp <- Stchr$CumLen
names(tmp) <- Stchr$Chr_ID

# making file with position and ID of SNP's
# grep -v "^##" pop_step03.recode.vcf | cut -f1-3 > icelandVCF.tsv

# Loading content of tsv file where Fst values are located
ok = read.delim('icelandVCF.tsv')

ok$CumPos <- ok$POS + tmp[ok$CHROM]
tmp <- Stchr$LG
names(tmp) <- Stchr$Chr_ID
ok$CHROM <- tmp[ok$CHROM]

BP <- ok$CumPos
P <- -log10(adj.p.values)
CHR <- ok$CHROM
RBP <- ok$POS

library(dplyr)
library(ggplot2)

my.df <- data.frame(BP, CHR, P, RBP)

axisdf = my.df %>% group_by(CHR) %>% summarize(center=( max(BP) + min(BP) ) / 2 )

# identifying SNP's
selected = -log10(adj.p.values)[candidates.qv]
mysnps <- my.df[my.df$P %in% selected,] 


ggplot(my.df, aes(x=BP, y=P)) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=1, size=2.0) +
  scale_color_manual(values = rep(c("#999999", "#98cce5"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center) +
  scale_y_continuous(expand = c(0.01, 0.01) ) +     # remove space between plot area and x axis

  # Add highlighted points
  geom_point(data=mysnps, aes(x=BP,y=P), color= "red", size=2.0) +
  
    #Add labels
    labs(x= "Chromosome", y= "-log10 (p values)") +
  
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
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
                     vjust = NULL)
    )


write.csv(mysnps, "mysnps.csv")


############################################################################
### Manhattan plot with FDR > 0.5
############################################################################
candidates.qv = which(qvalue(adj.p.values, fdr = 0.05)$signif)
candidates.qv

plot(-log10(adj.p.values), main="Manhattan plot", xlab = "Locus", cex = .7, col = "grey")

points(candidates.qv, -log10(adj.p.values)[candidates.qv], pch = 19, cex = .7, col = "red")


### making a Manhattan plot

Stchr <- read.csv(file="Salmo_trutta_chr_lengths.csv", stringsAsFactors=FALSE, sep="")
tmp <- Stchr$CumLen
names(tmp) <- Stchr$Chr_ID

# making file with position and ID of SNP's
# grep -v "^##" pop_step03.recode.vcf | cut -f1-3 > icelandVCF.tsv

# Loading content of tsv file where Fst values are located
ok = read.delim('icelandVCF.tsv')

ok$CumPos <- ok$POS + tmp[ok$CHROM]
tmp <- Stchr$LG
names(tmp) <- Stchr$Chr_ID
ok$CHROM <- tmp[ok$CHROM]

BP <- ok$CumPos
P <- -log10(adj.p.values)
CHR <- ok$CHROM
RBP <- ok$POS

library(dplyr)
library(ggplot2)

my.df <- data.frame(BP, CHR, P, RBP)

axisdf = my.df %>% group_by(CHR) %>% summarize(center=( max(BP) + min(BP) ) / 2 )

# identifying my SNP's
selected = -log10(adj.p.values)[candidates.qv]
mysnps <- my.df[my.df$P %in% selected,] 


ggplot(my.df, aes(x=BP, y=P)) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=1, size=2.0) +
  scale_color_manual(values = rep(c("#999999", "#98cce5"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center) +
  scale_y_continuous(expand = c(0.01, 0.01) ) +     # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=mysnps, aes(x=BP,y=P), color= "red", size=2.0) +
  
  #Add labels
  labs(x= "Chromosome", y= "-log10 (p values)") +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
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
                   vjust = NULL)
  )


write.csv(mysnps, "mysnps_fdr0.05.csv")

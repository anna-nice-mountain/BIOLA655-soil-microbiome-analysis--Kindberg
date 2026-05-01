# Microbial Community Analysis in R
# Author: Anna Kindberg
# Description: Performs statistical analysis for beta diversity, beta dispersion, PCoA ordination, and alpha diversity
# Input: Bray-Curtis distance matrix, shannon alpha diversity scores + metadata
# Output: PCoA plot and statistical test results

# R version: 4.5.2

#Download needed packages
library(ggplot2)
library(vegan)
library(car)

#
####make metadata
#young <- c("BFY7_R1", "BFY8_R1", "BFY9_R3", "HR1A_R1")
#mature <- c("FP5D_R1", "GSM1_R1", "GSM2_R1", "GSM3_R3", "GSM4_R1", "GSM5_R1", "MDM1_R1", "TRM1_R2", "UP4A_R1", "UP4B_R1", "UP4D_R1", "WCM4_R1")
samples <- rownames(beta_mat)
maturity <- ifelse(samples %in% c("BFY7_R1", "BFY8_R1", "BFY9_R3", "HR1A_R1"),
                   "young",
                   "mature")
meta <- data.frame(
  sample=samples,
  maturity = maturity)

rownames(meta) <-meta$sample

#
#####beta diversity analysis
beta <- read.delim("data_distance-matrix.tsv", row.names = 1, check.names = FALSE)
#convert to matrix
beta_mat <-as.matrix(beta)
#convert to distance object
beta_dist <- as.dist(beta_mat)

#PERMANOVA: p<0.05 = microbial community differences by maturity level could not have arisen by chance. 
adonis2(beta_dist ~  maturity, 
        data = meta,
        permutations = 999)

#beta dispersion test
disp <- betadisper(beta_dist, meta$maturity)
anova(disp) 

#PCoA plot
pcoa <- cmdscale(beta_dist, k = 2, eig = TRUE)
plot(pcoa$points[,1], pcoa$points[,2])
str(pcoa)

pcoa_df <- as.data.frame(pcoa$points)
pcoa_df$sample <- rownames(pcoa_df)
pcoa_df <- merge(pcoa_df, meta, by = "sample")

ggplot(pcoa_df, aes(x = V1, y = V2, color = maturity)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(x = "PCoA1", y = "PCoA2", title = "PCoA of Beta Diversity by Maturity Group")

#
######Alpha diversity analysis

alpha <- read.delim("data_alpha-diversity.tsv", row.names = 1, check.names = FALSE)
alpha$shannon_entropy <- alpha[,1]
meta$alpha <- alpha[rownames(meta), 1]

#Levenes Test 
leveneTest(alpha ~ maturity, data = meta) #test for homogeneity of variance

#Wilcoxon rank sum test for significance
# p<0.05 = Young and Mature forest soil communities differ significantly by alpha diversity
wilcox.test(Shannon ~ Maturity, data = meta_sub) 

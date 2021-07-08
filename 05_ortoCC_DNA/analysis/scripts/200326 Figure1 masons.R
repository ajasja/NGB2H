################################################# Standard libraries + viridis to make things pretty #################################################
library(viridis)
library(cowplot)
library(wesanderson)
library(GGally)
library(MASS)
library(dplyr)
library(tidyr)
library(extrafont)
################################################# Read in the data and mash it together #################################################
setwd("../data/")
load("masons.Rda")

map_mason1 <- read.table("mason-c1.x-y.txt", header = TRUE) 
map_mason1 <- filter(map_mason1, !is.na(X_peptide), !is.na(Y_peptide))%>%
  separate(X_peptide, into = c("X_peptide","X_codon","xjunk"), sep = "-") %>%
  separate(Y_peptide, into = c("Y_peptide","Y_codon","yjunk"), sep = "-") %>%
  select(-xjunk, -yjunk)

map_mason2 <- read.table("mason-c2.x-y.txt", header = TRUE) 
map_mason2 <- filter(map_mason2, !is.na(X_peptide), !is.na(Y_peptide))%>%
  separate(X_peptide, into = c("X_peptide","X_codon","xjunk"), sep = "-") %>%
  separate(Y_peptide, into = c("Y_peptide","Y_codon","yjunk"), sep = "-") %>%
  select(-xjunk, -yjunk)

mason_1D1 <- read.table("mason-1_DNA-1.sccount")
mason_1D2 <- read.table("mason-1_DNA-2.sccount")
mason_1R1 <- read.table("mason-1_RNA-1.sccount")
mason_1R2 <- read.table("mason-1_RNA-2.sccount")

mason_2D1 <- read.table("mason-2_DNA-1.sccount")
mason_2D2 <- read.table("mason-2_DNA-2.sccount")
mason_2R1 <- read.table("mason-2_RNA-1.sccount")
mason_2R2 <- read.table("mason-2_RNA-2.sccount")

colnames(mason_1D1) <- c("Barcode", "DNA_1_counts")
colnames(mason_1D2) <- c("Barcode", "DNA_2_counts")
colnames(mason_1R1) <- c("Barcode", "RNA_1_counts")
colnames(mason_1R2) <- c("Barcode", "RNA_2_counts")

colnames(mason_2D1) <- c("Barcode", "DNA_1_counts")
colnames(mason_2D2) <- c("Barcode", "DNA_2_counts")
colnames(mason_2R1) <- c("Barcode", "RNA_1_counts")
colnames(mason_2R2) <- c("Barcode", "RNA_2_counts")

mason1 <- full_join(mason_1D1, mason_1D2, by = "Barcode")
mason1 <- full_join(mason1, mason_1R1, by = "Barcode")
mason1 <- full_join(mason1, mason_1R2, by = "Barcode")
mason1[is.na(mason1)] <- 0

mason2 <- full_join(mason_2D1, mason_2D2, by = "Barcode")
mason2 <- full_join(mason2, mason_2R1, by = "Barcode")
mason2 <- full_join(mason2, mason_2R2, by = "Barcode")
mason2[is.na(mason2)] <- 0

mason_mapped1 <- inner_join(mason1, map_mason1, by = "Barcode")
mason_mapped2 <- inner_join(mason2, map_mason2, by = "Barcode")

med_mason1 <- group_by(filter(mason_mapped1, DNA_1_counts > 10, DNA_2_counts > 10), X_peptide, Y_peptide) %>%
  summarise(RD1 = median(RNA_1_counts/DNA_1_counts), RD2 = median(RNA_2_counts/DNA_2_counts), medRD_mas1 = median((RNA_1_counts + RNA_2_counts)/(DNA_1_counts + DNA_2_counts)), bcnum =n(), MAD = median( ((RNA_1_counts + RNA_2_counts)/(DNA_1_counts + DNA_2_counts)) - median((RNA_1_counts + RNA_2_counts)/(DNA_1_counts + DNA_2_counts))) )

med_mason2 <- group_by(filter(mason_mapped2, DNA_1_counts > 10, DNA_2_counts > 10), X_peptide, Y_peptide) %>%
  summarise(RD1 = median(RNA_1_counts/DNA_1_counts), RD2 = median(RNA_2_counts/DNA_2_counts), medRD_mas2 = median((RNA_1_counts + RNA_2_counts)/(DNA_1_counts + DNA_2_counts)), bcnum =n())

medboth <- inner_join(med_mason1, med_mason2, by = c("X_peptide", "Y_peptide"))

get_density <- function(x,y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

medboth <- data.frame(medboth, "density" = NA)
medboth$density <- get_density(log(medboth$medRD_mas1), log(medboth$medRD_mas2), n = 70)

med_mason1 <- data.frame(med_mason1, "density" = NA)
med_mason2 <- data.frame(med_mason2, "density" = NA)
med_mason1$density <- get_density(log(med_mason1$RD1), log(med_mason1$RD2), n = 70)
med_mason2$density <- get_density(log(med_mason2$RD1), log(med_mason2$RD2), n = 70)


################################################################################## Figure 1C bioreps
mas1repcor <- data.frame(x = -2, y = 0, label = paste("italic(r)==",round(cor(log(med_mason1$RD1),log(med_mason1$RD2)), digits = 3)))
g.medmason1biorep <- ggplot(med_mason1, aes(log(RD1), log(RD2), color = density)) + geom_point(size = 4, alpha = 0.5)  + scale_color_viridis(name = "Density") + scale_x_continuous(limits = c(-3.3, 1.7)) + scale_y_continuous(limits = c(-3.3, 1.7)) +
  labs(x = "Interaction score biological replicate 1", y = "Interaction score biological replicate 2", title  = "190611 Mason-1 biological replicates") + 
  geom_text(data = mas1repcor, aes(x=x,y=y,label=label), parse = T, color = "black", size =6, family = "Myriad Web Pro") + theme(axis.text = element_text(size = 14, family = "Myriad Web Pro"),text = element_text(size = 18, family = "Myriad Web Pro"))
g.medmason1biorep
ggsave(filename = "../figures/figure 1/190624 Mason-1 bioreps.pdf", plot = g.medmason1biorep, scale = 0.6)

################################################################################## Figure 1D mason barcode replicates
bccor <- data.frame(x = -2.5, y = 0, label = paste("italic(r)==",round(cor(log(medboth$medRD_mas1),log(medboth$medRD_mas2)), digits = 3)))
g.bccor <- ggplot(medboth, aes(log(medRD_mas1), log(medRD_mas2), color = density)) + geom_point(size = 4, alpha = 0.7) + scale_x_continuous(limits = c(-3.3, 2), breaks = c(-3, -2, -1, 0, 1, 2)) + 
  scale_y_continuous(limits = c(-3.3, 2), breaks = c(-3, -2, -1, 0, 1, 2)) + scale_color_viridis(name = "Density") +
  labs(x = "Interaction score barcoding 1", y = "Interaction score barcoding 2", title  = "190611 Mason comparisons between different barcoding schemes") +
  geom_text(data = bccor, aes(x=x,y=y,label=label), parse = T, color = "black", size =6, family = "Myriad Web Pro") + 
  theme(axis.text = element_text(size = 14, family = "Myriad Web Pro"),text = element_text(size = 18, family = "Myriad Web Pro"))
g.bccor

ggsave(filename = "../figures/figure 1/190611 Mason barcoding replicates.pdf", plot = g.bccor, scale = 0.6)

################################################################################## Figure 1E Codon correlation
med_mas1_codon <- group_by(filter(mason_mapped1, DNA_1_counts > 10, DNA_2_counts > 10), X_peptide, Y_peptide, X_codon, Y_codon) %>%
  summarise(medRD = median((RNA_1_counts + RNA_2_counts)/(DNA_1_counts + DNA_2_counts)), bcnum = n())

med_C1C6 <- filter(med_mas1_codon, (X_codon == 'C2' & Y_codon == 'C3') | (X_codon == 'C1' & Y_codon == 'C1')) 
med_C1C6 <- ungroup(med_C1C6)
med_C1C6 <-  mutate(med_C1C6, "codons" = paste(X_codon, Y_codon, sep = "")) %>%
  select(-X_codon, -Y_codon, -bcnum)
med_C1C6 <- spread(med_C1C6, codons, medRD) %>%
  filter(!is.na(C1C1), !is.na(C2C3))
med_C1C6 <- data.frame(med_C1C6, "density" = NA)
med_C1C6$density <- get_density(log(med_C1C6$C1C1), log(med_C1C6$C2C3), n = 70)

medcodcorr <- data.frame(x = -2.5, y = 0, label = paste("italic(r)==",round(cor(log(med_C1C6$C1C1),log(med_C1C6$C2C3), use = "complete.obs"), digits = 3)))
g.medmedcodon <- ggplot(med_C1C6, aes(log(C1C1), log(C2C3), color = density)) + geom_point(size = 4, alpha = 0.7) + scale_x_continuous(limits = c(-4.1, 1.9), breaks = c(-4, -3, -2, -1, 0, 1)) + 
  scale_y_continuous(limits = c(-4.1, 1.9), breaks = c(-4, -3, -2, -1, 0, 1)) +labs(x = "Interaction score codon usage 1", y = "Interaction score codon usage 6", title = "190702 Mason codon usage replicate") +
  geom_text(data = medcodcorr, aes(x=x,y=y,label=label),parse = TRUE, size = 6, color = "black") + 
  theme(text = element_text(family= "Myriad Web Pro", size = 18), axis.text = element_text(size = 14, family = "Myriad Web Pro")) + scale_color_viridis(name = "Density")
g.medmedcodon

ggsave(filename = "../figures/figure 1/190704 Mason Codon Usage replicates.pdf", plot = g.medmedcodon, scale  = 0.6)


################################################################################## Figure 1F mason reciprical correlations
med_mason1$X_peptide <- factor(med_mason1$X_peptide, levels = c("O1","O2","O3","O4","O5","O6","O7","O8","O9","O10","O11","O12","O13","O14","O15","O16"))
med_mason1$Y_peptide <- factor(med_mason1$Y_peptide, levels = c("O1","O2","O3","O4","O5","O6","O7","O8","O9","O10","O11","O12","O13","O14","O15","O16"))
mason1rev <- data.frame(med_mason1, altor = 0, altbcs = 0)

for (row in 1:nrow(mason1rev)){
  if(as.integer(mason1rev$X_peptide[row]) - as.integer(mason1rev$Y_peptide[row]) > 0) { 
    a <- filter(mason1rev, X_peptide == mason1rev[row, 2], Y_peptide == mason1rev[row, 1])
    mason1rev$altor[row] <- a$medRD_mas1
  }
}

scoredmason1rev <- filter(mason1rev, altor != 0)
scoredmason1rev <- mutate(scoredmason1rev, "density" = NA)
scoredmason1rev$density <- get_density(log(scoredmason1rev$medRD_mas1), log(scoredmason1rev$altor), n = 70)

corrmas1recip <- data.frame(x = -2, y = 1, label = paste("italic(r)==", round(cor(log(scoredmason1rev$medRD_mas1), log(scoredmason1rev$altor)), digits = 3)))
g.mason1recip <- ggplot(scoredmason1rev, aes(log(medRD_mas1), log(altor), color = density)) + geom_point(alpha = 0.7, size =4) + geom_text(data = corrmas1recip, aes(x=x, y=y,label=label), parse = TRUE, size  = 6, color = "black") +
  labs(x =  "Interaction score", y = "Reciprocal orientation interaction score", title = "190611 Mason-1 reciprocal orientations median RNA/DNA") +
  theme(axis.text = element_text(size = 14, family = "Myriad Web Pro"),text = element_text(size = 18, family = "Myriad Web Pro")) + geom_abline(slope = 1, color = "blue") + scale_color_viridis(name = "Density") +
  scale_x_continuous(limits = c(-3.2, 1.4)) + scale_y_continuous(limits = c(-3.2, 1.4))
g.mason1recip

ggsave(filename = "../figures/figure 1/190624 Mason-1 recip orientations.pdf", plot = g.mason1recip, scale = 0.6)


################################################################################## Figure 1G mason data as Tms
mas_tmscut <- mason_tms

mas_tmscut$X_peptide <- factor(mas_tmscut$X_peptide, levels = c("O1","O2","O3","O4","O5","O6","O7","O8","O9","O10","O11","O12","O13","O14","O15","O16"))
mas_tmscut$Y_peptide <- factor(mas_tmscut$Y_peptide, levels = c("O1","O2","O3","O4","O5","O6","O7","O8","O9","O10","O11","O12","O13","O14","O15","O16"))

a <- c()
for (row in 1:nrow(mas_tmscut)){
  if(as.integer(mas_tmscut$X_peptide[row]) - as.integer(mas_tmscut$Y_peptide[row]) > 0) { 
    a <- c(a,row)
  }
}
mas_tmscut <- mas_tmscut[-a, ]
mas_tmscut$Tm[mas_tmscut$Tm < 40] <- 40

g.mas_tmgrid <- ggplot(mas_tmscut, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill=Tm)) + scale_fill_viridis(name = "Previously\npublished Tm") +
  labs(x = "X protein", y = "Y protein", title = "190716 Mason Tms > 40 grid graph") + theme(text = element_text(family = "Myriad Web Pro", size = 14), axis.title = element_text(size = 18), axis.text.x = element_text(angle = 90))
g.mas_tmgrid

ggsave(filename = "../figures/figure 1/200326 Mason Tms greater than 40.pdf", plot = g.mas_tmgrid, scale = 0.6)

################################################################################## Figure 1H tile graph
med_mason1$X_peptide <- factor(med_mason1$X_peptide, levels = c("O1","O2","O3","O4","O5","O6","O7","O8","O9","O10","O11","O12","O13","O14","O15","O16"))
med_mason1$Y_peptide <- factor(med_mason1$Y_peptide, levels = c("O1","O2","O3","O4","O5","O6","O7","O8","O9","O10","O11","O12","O13","O14","O15","O16"))

g.masongrid <- ggplot(med_mason1, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill = log(medRD_mas1))) + scale_fill_viridis(name = "Interaction\nscore") +
  labs(x = "X protein", y = "Y protein", title = "190613 Mason-1 median RNA/DNA") + theme(text = element_text(family = "Myriad Web Pro", size = 14), axis.title = element_text(size = 18), axis.text.x = element_text(angle = 90))
g.masongrid

ggsave(filename = "../figures/figure 1/190613 mason 1 tile graph.pdf", plot = g.masongrid, scale =0.6)
################################################################################ Figure 1I
corrlist <- data_frame("corr" = 1,"i" = 0, "rep"=0)
for(j in 1:5){
  for(i in seq(10,160, 10)){
    down_masD1 <- sample_n(mason_1D1, size = sum(mason_1D1$DNA_1_counts)/i, replace = T, weight = DNA_1_counts)
    down_masD1 <- group_by(down_masD1, Barcode) %>%
      summarise(DNA_1_counts = n())
    down_masD2 <- sample_n(mason_1D2, size = sum(mason_1D2$DNA_2_counts)/i, replace = T, weight = DNA_2_counts)
    down_masD2 <- group_by(down_masD2, Barcode) %>%
      summarise(DNA_2_counts = n())
    down_masR1 <- sample_n(mason_1R1, size = sum(mason_1R1$RNA_1_counts)/i, replace = T, weight = RNA_1_counts)
    down_masR1 <- group_by(down_masR1, Barcode) %>%
      summarise(RNA_1_counts = n())
    down_masR2 <- sample_n(mason_1R2, size = sum(mason_1R2$RNA_2_counts)/i, replace = T, weight = RNA_2_counts)
    down_masR2 <- group_by(down_masR2, Barcode) %>%
      summarise(RNA_2_counts = n())
    
    downsamp <- full_join(down_masD1, down_masD2, by = "Barcode")
    downsamp <- full_join(downsamp, down_masR1, by = "Barcode")
    downsamp <- full_join(downsamp, down_masR2, by = "Barcode")
    downsamp[is.na(downsamp)] <- 0
    downsamp <- inner_join(downsamp,map_mason1, by = "Barcode")
    
    med_downsamp <- group_by(filter(downsamp, DNA_1_counts > 10, DNA_2_counts > 10), X_peptide, Y_peptide) %>%
      summarise(RD1 = median(RNA_1_counts/DNA_1_counts), RD2 = median(RNA_2_counts/DNA_2_counts), medRD_down = median((RNA_1_counts + RNA_2_counts)/(DNA_1_counts + DNA_2_counts)), bcnum =n())
    med_downsamp[med_downsamp == 0] <- 0.005
    
    med_downsamp <- inner_join(med_downsamp, med_mason1, by = c("X_peptide", "Y_peptide"))
    corrdownsamp <- cor(log(med_downsamp$medRD_down),log(med_downsamp$medRD_mas1))
    corradd <- data.frame("corr" = corrdownsamp, "i" = i, "rep" =j )
    corrlist <- bind_rows(corrlist, corradd) 
  }
}
sumcorrlist <- group_by(corrlist, i) %>%
  summarise(meancorr = mean(corr), med_corr = median(corr), sdcorr = sd(corr))

sdcorr <- aes(ymin = meancorr -sdcorr, ymax = meancorr + sdcorr)
g.sumcorrlist <- ggplot(filter(sumcorrlist, i <160), aes(i, meancorr)) + geom_bar(stat = "identity", fill ="#00A08A") + geom_errorbar(sdcorr) + #scale_fill_manual(values = c("#00A08A")) +
  labs(x = "Fold sub-sampled", y = "Correlation with complete data", title = "191220 Mason-1 median RNA/DNA subsampling correlation with full data set") +
  theme(axis.text = element_text(size = 14, family = "Myriad Web Pro"),text = element_text(size = 18, family = "Myriad Web Pro"))
g.sumcorrlist


g.sumcorrlist <- ggplot(filter(corrlist, i < 160), aes(as.factor(i), corr)) + geom_boxplot() +labs(x = "Fraction of total data used", y = "Interaction score\ncorrelation with all data", title = "191220 Mason-1 median RNA/DNA subsampling correlation with full data set") +
  theme(axis.text = element_text(size = 14, family = "Myriad Web Pro"),text = element_text(size = 18, family = "Myriad Web Pro")) + scale_y_continuous(limits = c(0, 1)) +
  scale_x_discrete(labels = c("1","1/20","1/40","1/60","1/80","1/100","1/120","1/140"), breaks = c("0", "20", "40", "60", "80", "100", "120", "140"))
g.sumcorrlist

ggsave(filename = "../figures/figure 1/191220 Mason downsamples.pdf", plot = g.sumcorrlist, scale =0.6)



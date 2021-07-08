################################################# Standard libraries + viridis to make things pretty + extrafont for comic sans #################################################
library(viridis)
library(cowplot)
library(wesanderson)
library(GGally)
library(MASS)
library(extrafont)
library(dplyr)
library(tidyr)
################################################# Read in the data and mash it together #################################################
setwd("../data")
load("masons.Rda")
load("R18k_with_aliases.Rda")

#map_18k <- read.table("R18k-c1.x-y.txt", header = TRUE) 
#map_18k <- filter(map_18k, !is.na(X_peptide), !is.na(Y_peptide))%>%
#  separate(X_peptide, into = c("X_peptide","junkx"), sep = "\\+") %>%
#  separate(Y_peptide, into = c("junky", "Y_peptide"), sep = "\\+") %>%
#  separate(Y_peptide, into = c("Y_peptide", "Y_codon"), sep = -3) %>%
#  select(-junkx, -junky)

R18k_D1 <- read.table("R18k-pDNA-1.sccount")
R18k_D2 <- read.table("R18k-pDNA-2.sccount")
R18k_R1 <- read.table("R18k-cDNA-1.sccount")
R18k_R2 <- read.table("R18k-cDNA-2.sccount")

colnames(R18k_D1) <- c("Barcode", "DNA_1_counts")
colnames(R18k_D2) <- c("Barcode", "DNA_2_counts")
colnames(R18k_R1) <- c("Barcode", "RNA_1_counts")
colnames(R18k_R2) <- c("Barcode", "RNA_2_counts")

R18k <- full_join(R18k_D1, R18k_D2, by = "Barcode")
R18k <- full_join(R18k, R18k_R1, by = "Barcode")
R18k <- full_join(R18k, R18k_R2, by = "Barcode")
R18k[is.na(R18k)] <- 0

R18k_map <- inner_join(R18k, map_18k, by = "Barcode")

med18k <- group_by(filter(R18k_map, DNA_1_counts > 10, DNA_2_counts > 10), X_peptide, Y_peptide) %>%
  summarise(medRD1 = median(RNA_1_counts/DNA_1_counts), medRD2 = median(RNA_2_counts/DNA_2_counts), medRDboth = median((RNA_1_counts + RNA_2_counts)/(DNA_1_counts + DNA_2_counts)), bcnum = n(),
           CIhigh = median((RNA_1_counts + RNA_2_counts)/(DNA_1_counts + DNA_2_counts)) + (1.96 * sd((RNA_1_counts + RNA_2_counts)/(DNA_1_counts + DNA_2_counts))/sqrt(n())), 
           CIlow = median((RNA_1_counts + RNA_2_counts)/(DNA_1_counts + DNA_2_counts)) - (1.96 * sd((RNA_1_counts + RNA_2_counts)/(DNA_1_counts + DNA_2_counts))/sqrt(n())))
#med18k <- separate(med18k, X_peptide, into = c("X_peptide","X_group"), sep = -3)
#med18k <- separate(med18k, Y_peptide, into = c("Y_peptide","Y_group"), sep = -3)

get_density <- function(x,y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

med18k$density[med18k$medRD1 != 0 & med18k$medRD2 != 0] <- get_density(log(filter(med18k, medRD1 != 0,medRD2 != 0)$medRD1),log(filter(med18k, medRD1 != 0,medRD2 != 0)$medRD2), n = 1000)
med18k$density[is.na(med18k$density)] <- 0


med18k <-  mutate(med18k, "interact" = paste(X_peptide, Y_peptide, sep = ""), "revint" = paste(Y_peptide, X_peptide, sep = ""))
med18k <- data.frame(med18k, classified = 0, altor = 0, altbcs = 0, altCIhi = 0, altCIlow = 0)

for (row in 1:nrow(med18k)) {
  a <- med18k$interact[row]
  b <-  filter(med18k, revint == a)$medRDboth 
  c <-  filter(med18k, revint == a)$bcnum
  if (length(b) != 0) {
    med18k$altor[row] <- b
    med18k$altbcs[row] <- c
  }
} 

med18k$CIhigh[is.na(med18k$CIhigh)] <- 0
med18k$CIlow[is.na(med18k$CIlow)] <- 0


for (row in 1:nrow(med18k)) {
  a <- med18k$interact[row]
  revintrow <-  filter(med18k, revint == a)
  if (length(revintrow$X_peptide) != 0) {
    med18k$altCIhi[row] <- revintrow$CIhigh
    med18k$altCIlow[row] <- revintrow$CIlow
  }
} 
med18k <- mutate(med18k, medavg = (medRDboth*bcnum + altor*altbcs)/(bcnum + altbcs))

med18k <- mutate(med18k, CIhighavg = (CIhigh*bcnum + altCIhi*altbcs)/(bcnum + altbcs))
med18k <- mutate(med18k, CIlowavg = (CIlow*bcnum + altCIlow*altbcs)/(bcnum + altbcs))
med18k <- filter(med18k, !(bcnum == 1 & altbcs == 1) & !(bcnum == 1 & altbcs == 0))

med18ksnone <- filter(med18k, altbcs == 0, bcnum != 1)
colnames(med18ksnone) <- c("Y_peptide", "X_peptide","medRD1","medRD2","medRDboth","altbcs", 'CIhigh', 'CIlow',"density","revint","interact","classified","altor","bcnum", "altCIhi", "altCIlow","medavg", 'CIhighavg','CIlowavg')
med18ksnone <- select(med18ksnone, X_peptide, Y_peptide, medRD1, medRD2, medRDboth, bcnum, CIhigh, CIlow, density,interact, revint, classified, altor, altbcs, altCIhi, altCIlow, medavg, CIhighavg, CIlowavg)
med18k <- bind_rows(med18k, med18ksnone)

small18k <- select(med18k, X_peptide, Y_peptide, CIhighavg, CIlowavg)
write.csv(small18k, file = 'R18kCIavgs.csv', quote = FALSE, row.names = FALSE)

maxset <- read.csv(file = "210603_R18k_largest_orthogonal_subset_gap0.csv", header = T, quote = "\'")
maxset <- gather(maxset,key = "peptidenum",value = "pepname", -name, -classifier, -minpos,-maxneg, -lnogap)
maxset$pepname <- gsub("\\[|\\]", "", maxset$pepname)
maxset <- filter(maxset, pepname != "", pepname != " ")
maxset$name <- factor(maxset$name)
maxset$pepname <- trimws(maxset$pepname)
summaxset <- group_by(maxset, name, classifier) %>%
  summarise(numpeps = n(), ogap = mean(lnogap), minpos = mean(minpos), maxneg = mean(maxneg))


#####################################################################################################################################################
#####################################                                  Figure 4B                           ##################################### 
#####################################################################################################################################################

numints <- data.frame()
for (i in levels(maxset$name)){
  for (j in unique(maxset$classifier)){
    peps <- filter(maxset, name == i, classifier == j)$pepname
    interactions <- filter(med18k, X_peptide %in% peps, Y_peptide %in% peps)
    interactions <- mutate(interactions, "classy" = 0)
    interactions$classy[interactions$medavg > j & (interactions$X_peptide != interactions$Y_peptide)] <- 0.5
    interactions$classy[interactions$medavg > j & (interactions$X_peptide == interactions$Y_peptide)] <- 1
    row <- data.frame(i, j, sum(interactions$classy))
    numints <- bind_rows(numints, row)
  } 
}
colnames(numints) <- c("name","classifier","orthints")

numints <- left_join(summaxset, numints, by = c("name", "classifier"))
sumnumints <- group_by(numints, name) %>%
  summarise(maxos = max(orthints))
numints2 <- left_join(numints, sumnumints, by = "name")
numints2 <- numints2 %>% filter(maxos == orthints) %>% group_by(name) %>% filter(ogap == max(ogap)) %>% filter(classifier == last(classifier))


g.sumnumints <- ggplot(sumnumints, aes(maxos)) + geom_histogram(fill ="#00A08A", bins = 13, color = "grey70") + #scale_x_continuous(limits = c(0,25)) + labs(x = "Number of orthogonal interactions", y = "Count", title = "200310 R18k Number of orthogonal interactions per set") +
  theme(text = element_text(family = "Myriad Web Pro"))
g.sumnumints

ggsave(file = "../figures/figure 4/200529 Figure 4B number of ortho interactions.pdf", scale = 0.6)


#####################################################################################################################################################
#####################################                                  Figure 4C                                ##################################### 
#####################################################################################################################################################


listofmostpeps <- filter(maxset, name == filter(numints2, orthints == max(numints2$orthints))$name, classifier == filter(numints2, orthints == max(numints2$orthints))$classifier)$pepname
mostints <- filter(med18k, X_peptide %in% listofmostpeps, Y_peptide %in% listofmostpeps)
mostints <- separate(mostints, X_peptide, into = c("X_peptide","X_junk"), sep = "-") %>%
  separate(Y_peptide, into = c("Y_peptide", "Y_junk"), sep = "-")

mostints$X_peptide <- factor(mostints$X_peptide, levels = c("4H3697", "4H1142", "4H1579","4H1864", "4H1597","4H3410","4H171","4H3016","4H485","4H2725",
                                                            "4H3395", "4H3120","4H536", "4H891", "4H1385", "4H1626", "4H1842", "4H2162","4H2253","4H2274",
                                                            "4H2330", "4H2698","4H310","4H3165","4H629","4H797", "4H909",  "4H94", "4H930", "4H3590","4H2438"))

mostints$Y_peptide <- factor(mostints$Y_peptide, levels =  c("4H3697", "4H1142", "4H1579","4H1864", "4H1597","4H3410","4H171","4H3016","4H485","4H2725",
                                                             "4H3395", "4H3120","4H536", "4H891", "4H1385", "4H1626", "4H1842", "4H2162","4H2253","4H2274",
                                                             "4H2330", "4H2698","4H310","4H3165","4H629","4H797", "4H909",  "4H94", "4H930", "4H3590","4H2438"))

mostints$medavg[mostints$medavg < 1.64] <- 1.64
g.mostints <- ggplot(mostints, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill = log(medavg))) + geom_tile(data = filter(mostints, medavg > 4.9), color = "grey70", alpha = 0, size = 1.5)+ scale_fill_viridis(name = "Interaction\nscore", labels = c('<0.5','1.0','1.5','2.0','2.5')) + theme(axis.text.x = element_text(angle = 90), text = element_text(family = "Myriad Web Pro")) +
  labs(x = "X protein", y = "Y protein", title = "200310 R18k Largest set orthogonal set") 
g.mostints

ggsave(g.mostints, filename = "../figures/Figure 4/200310 Figure 4C R18k largest orthogonal set.pdf", scale = 0.6)



#####################################################################################################################################################
#####################################                                  Figure 4D                                ##################################### 
#####################################################################################################################################################
m18chop <- separate(med18k, X_peptide, into = c("X_peptide","xjunk"), by ="-") %>%
  separate(Y_peptide, into = c("Y_peptide","yjunk"), by = "-")
masontm <- left_join(mason_tms, m18chop)

sumsummaxset <- group_by(numints, orthints) %>%
  filter(ogap == max(ogap)) %>%
  select(-classifier) %>%
  distinct()

sumsummaxset$orthints <- factor(sumsummaxset$orthints, levels = c("0","1", "2", "3", "4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"))

preds <-  lm(masontm$Tm ~ log(masontm$medavg))

#g.orthogapsbysetsize <- ggplot(filter(sumsummaxset, !is.na(orthints), orthints != 0, orthints != 1) , aes(orthints, log(minpos-maxneg))) + geom_bar(stat = "identity", fill = "#00A08A", color = "grey70") + labs(x = "Number of orthogonal proteins per set", y = "Orthogonality gap", title = "200304 R18k Maximum Size of orthogonality gap by number of peptides present") + 
#  annotate("text", x= " ", y = c(seq(0,3,1)), label = as.character(round(seq(0, 3, 1)*preds$coefficients["log(masontm$medavg)"] + preds$coefficients['(Intercept)'], 1))) + scale_x_discrete(limits = c(levels(sumsummaxset$orthints),'_' ,' ', "  ")) + annotate("line", x = "_", y =c(0,3)) + 
#  annotate("text", x ="  ", y = 0.5, label = "Estimated Tm", angle = 90, size = 5 ) + theme(text = element_text(family = "Myriad Web Pro"))
#g.orthogapsbysetsize

g.orthogapsbysetsize <- ggplot(filter(sumsummaxset, !is.na(orthints), orthints != 0, orthints != 1) , aes(orthints, ogap)) + geom_bar(stat = "identity", fill = "#00A08A", color = "grey70") + labs(x = "Number of orthogonal proteins per set", y = "Orthogonality gap", title = "200304 R18k Maximum Size of orthogonality gap by number of peptides present") + 
    scale_x_discrete(limits = c(levels(sumsummaxset$orthints)))  + theme(text = element_text(family = "Myriad Web Pro"))
g.orthogapsbysetsize

ggsave(g.orthogapsbysetsize, filename = "../figures/figure 4/200310 Figure 4D R18k maxset size by orthogap ints2.pdf", scale = 0.6)


ggplot(masontm, aes(log(medavg), Tm)) + geom_point(alpha = 0.4) + geom_abline(intercept = 34.06, slope = 11.43)

#####################################################################################################################################################
#####################################                                  Figure 4e                                ##################################### 
#####################################################################################################################################################
modelscore <-read.table(file = "combined_scoring.csv", sep = ",", colClasses = c("NULL","character","character", rep("NULL", 7), "character",rep("NULL", 4)), header = T)
colnames(modelscore) <- c("X_peptide", "Y_peptide","D5_b_nc")
trim18k <- separate(med18k, X_peptide,into = c("X_peptide","xjunk"),  sep = "-") %>%
  separate(Y_peptide, into = c("Y_peptide", "yjunk"), sep = "-")
modelscore <- inner_join(trim18k, modelscore, by = c("X_peptide", "Y_peptide"))
modelscore$D5_b_nc <- as.numeric(modelscore$D5_b_nc)

get_density <- function(x,y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

modelscore2 <- data.frame(modelscore, "density" = "")
modelscore2$density <- get_density(log(modelscore2$medavg), modelscore2$D5_b_nc, n = 100)

icipacor <- data.frame(y = 0-2, x = 50, label =paste("italic(R^2)==", round(cor(log(modelscore$medavg),modelscore$D5_b_nc)^2, digits = 3)))
g.icipa <- ggplot(modelscore2, aes( D5_b_nc,log(medavg))) + geom_point(aes(color = density), alpha = 0.5) +  labs(y = "Interaction score", x = "iCipa predicted Tm", title = "200310 R18k iCipa predicted vs medRD") +
  geom_smooth(method = "lm", se = FALSE, color = "Black")  + theme(text = element_text(family = "Myriad Web Pro")) + geom_text(data = icipacor, aes(x=x,y=y,label=label), parse = T) + scale_color_viridis(name = "Density") +
  scale_y_continuous(limits = c(-3,2.8))
g.icipa

ggsave(g.icipa, filename = "../figures/figure 4/200310 Figure 4E R18k medavg icipa corrleation.pdf", scale = 0.6)


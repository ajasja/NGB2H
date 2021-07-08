################################################# Standard libraries + viridis to make things pretty + extrafont for comic sans #################################################
library(viridis)
library(cowplot)
library(wesanderson)
library(extrafont)
library(GGally)
library(MASS)
library(dplyr)
library(tidyr)
################################################# Read in the data and mash it together #################################################
setwd("../data/")


map_R8000 <- read.table("R8000-c1.x-y.txt", header = TRUE) 
map_R8000 <- separate(map_R8000, X_peptide, into = c("X_peptide","junkx"), sep = "\\+") %>%
  separate(Y_peptide, into = c("junky", "Y_peptide"), sep = "\\+") %>%
  separate(Y_peptide, into = c("Y_peptide", "Y_codon"), sep = "_") %>%
  filter(!is.na(X_peptide), !is.na(Y_peptide)) %>%
  select(-junkx, -junky)


R8000_D1 <- read.table("R8000-pDNA-1.sccount")
R8000_D2 <- read.table("R8000-pDNA-2.sccount")
R8000_R1 <- read.table("R8000-cDNA-1.sccount")
R8000_R2 <- read.table("R8000-cDNA-2.sccount")

colnames(R8000_D1) <- c("Barcode", "DNA_1_counts")
colnames(R8000_D2) <- c("Barcode", "DNA_2_counts")
colnames(R8000_R1) <- c("Barcode", "RNA_1_counts")
colnames(R8000_R2) <- c("Barcode", "RNA_2_counts")

R8000 <- full_join(R8000_D1, R8000_D2, by = "Barcode")
R8000 <- full_join(R8000, R8000_R1, by = "Barcode")
R8000 <- full_join(R8000, R8000_R2, by = "Barcode")
R8000[is.na(R8000)] <- 0

R8000_map <- inner_join(R8000, map_R8000, by = "Barcode")

medR8000 <- group_by(filter(R8000_map, DNA_1_counts > 10, DNA_2_counts > 10), X_peptide, Y_peptide) %>%
  summarise(medRD1 = median(RNA_1_counts/DNA_1_counts), medRD2 = median(RNA_2_counts/DNA_2_counts), medRDboth = median((RNA_1_counts + RNA_2_counts)/(DNA_1_counts + DNA_2_counts)), bcnum = n(), sumD1 = sum(DNA_1_counts), sumD2 = sum(DNA_2_counts), sumDboth = sum(DNA_1_counts + DNA_2_counts))
medR8000 <- separate(medR8000, X_peptide, into  = c("X_peptide", "X_group"), sep = "-") %>%
  separate(Y_peptide, into = c("Y_peptide","Y_group"), sep = "-") %>%
  filter(X_group == Y_group | (is.na(X_group) & is.na(Y_group)))

#take the data on largest orthogonal sets and get it into tidy form
maxset <- read.csv(file = "200305_R8000_max_ortho_sets.csv", header = T, quote = "\'")
maxset <- gather(maxset,key = "peptidenum",value = "pepname", -name, -classifier, -orthogap)
maxset$pepname <- gsub("\\[|\\]", "", maxset$pepname)
maxset <- filter(maxset, pepname != "")
maxset <- separate(maxset, name, into = c("name", "group"), sep = -3) %>%
  separate(pepname, sep = "-", into = c("pepname","pepgroup"))
maxset <- separate(maxset, group, into =c("junk","group"), sep = "-") %>%
  select(-junk)
maxset$name <- factor(maxset$name)
maxset$pepname <- trimws(maxset$pepname)

R8000maxs <- full_join(medR8000, maxset, by = c("X_group" = "group", "X_peptide" = "pepname"))
R8000maxs <- inner_join(R8000maxs, maxset, by = c("Y_group" = "group", "Y_peptide" = "pepname", "name", "classifier", "orthogap", "pepgroup"))
R8000maxs <- filter(R8000maxs, !is.na(name))
R8000maxs <- data.frame(R8000maxs, "classy" = 0)
R8000maxs$classy[R8000maxs$medRDboth > R8000maxs$classifier] <- 1
###########################################################################################
############################## Figure 2B R8000 set with most interacting orthogonal interactions
###########################################################################################
orthos <- filter(R8000maxs, name == "4h-1or2N-not_first_only_B07_bc-8.05_nc-7.05-all.00.fasta", X_group == "bN", classifier == 2.3)
orthos2 <- filter(R8000maxs, name == "4h-1or2N-not_first_only_B07_bc-8.05_nc-7.05-all.00.fasta", X_group == "bN", classifier == 2.3)
colnames(orthos2) <- c("Y_peptide","X_group","X_peptide","Y_group","medRD1","medRD2","medRDboth","bcnum","sumD1","sumD2","sumDboth", "name","classifier","orthogap","peptidenum.x","pepgroup","peptidenum.y", "classy") 
orthos2  <- select(orthos2, "X_peptide","X_group","Y_peptide","Y_group","medRD1","medRD2","medRDboth","bcnum","sumD1","sumD2","sumDboth", "name","classifier","orthogap","peptidenum.x","pepgroup","peptidenum.y", "classy") 
orthos <- bind_rows(orthos, orthos2)
orthos <- orthos[!duplicated(orthos),]

orthos$X_peptide <- factor(orthos$X_peptide, levels = c( "4H235","4H2760","4H264","4H2931","4H166","4H3297","4H2585","4H349","4H2162","4H2253","4H2444","4H3364","4H56","4H730","4H796","4H932"))
orthos$Y_peptide <- factor(orthos$Y_peptide, levels = c("4H235","4H2760","4H264","4H2931","4H166","4H3297","4H2585","4H349","4H2162","4H2253","4H2444","4H3364","4H56","4H730","4H796","4H932"))

onesideorthos <- data.frame()
for (i in 1:nrow(orthos)) {
  if(as.integer(orthos$X_peptide[i]) - as.integer(orthos$Y_peptide[i]) <= 0) {
    onesideorthos <- bind_rows(onesideorthos, orthos[i,])
  }
}

onesideorthos$medRDboth[onesideorthos$medRDboth < 1] <- 1

g.orthos <- ggplot(onesideorthos, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill =log(medRDboth))) + theme(axis.text.x = element_text(angle = 90), text = element_text(family = "Myriad Web Pro", size = 14), axis.title = element_text(size = 18)) +
  scale_fill_viridis(name = "Interaction\nScore") + labs(x = "X protein", y = "Y protein", title = "190617 R8000 largest number of orthogonal interactions") + geom_tile(data = filter(onesideorthos, classy == 1), alpha =0, size = 2.0, color = "grey")
g.orthos

ggsave(filename = "../figures/figure 2/191108 R8000 most orthogonal interactions.pdf", plot = g.orthos, scale = 0.6)

#################################################################################
############################## Figure 2C set with most R8000 orthogonal peptides
###########################################################################################
biggroup <- filter(R8000maxs, name == "4h-1or2N-not_first_only_B07_bc-8.05_nc-7.05-hetero-ex.01.00.fasta", X_group == "bN", classifier == 2.7)
biggroup2 <- filter(R8000maxs, name == "4h-1or2N-not_first_only_B07_bc-8.05_nc-7.05-hetero-ex.01.00.fasta", X_group == "bN", classifier == 2.7)
colnames(biggroup2) <- c("Y_peptide","X_group","X_peptide","Y_group","medRD1","medRD2","medRDboth","bcnum","sumD1","sumD2","sumDboth", "name","classifier","orthogap","peptidenum.x","pepgroup","peptidenum.y", "classy") 
biggroup2  <- select(biggroup2, "X_peptide","X_group","Y_peptide","Y_group","medRD1","medRD2","medRDboth","bcnum","sumD1","sumD2","sumDboth",  "name","classifier","orthogap","peptidenum.x","pepgroup","peptidenum.y", "classy") 
biggroup <- bind_rows(biggroup, biggroup2)
biggroup <- biggroup[!duplicated(biggroup),]

biggroup$X_peptide <- factor(biggroup$X_peptide, levels = c("4H2503","4H2223","4H182","4H3057","4H2839","4H319","4H2588","4H348","4H2915","4H40","4H2717","4H477","4H2692","4H492","4H2134","4H3277","4H741"))
biggroup$Y_peptide <- factor(biggroup$Y_peptide, levels = c("4H2503","4H2223","4H182","4H3057","4H2839","4H319","4H2588","4H348","4H2915","4H40","4H2717","4H477","4H2692","4H492","4H2134","4H3277","4H741"))

onesidebiggroup <- data.frame()
for (i in 1:nrow(biggroup)) {
  if(as.integer(biggroup$X_peptide[i]) - as.integer(biggroup$Y_peptide[i]) <= 0) {
    onesidebiggroup <- bind_rows(onesidebiggroup, biggroup[i,])
  }
}

onesidebiggroup$medRDboth[onesidebiggroup$medRDboth < 1] <- 1
  
g.biggroup <- ggplot(onesidebiggroup, aes(X_peptide, Y_peptide)) + geom_tile(aes(fill =log(medRDboth))) + theme(axis.text.x = element_text(angle = 90), text = element_text(family = "Myriad Web Pro", size = 14), axis.title = element_text(size = 18)) +
  scale_fill_viridis(name = "Interaction\nscore") + labs(x = "X protein", y = "Y protein", title = "190617 R8000 largest number of peptides") + geom_tile(data = filter(onesidebiggroup, classy == 1), alpha =0, size = 2.0, color = "grey")
g.biggroup
ggsave(filename = "../figures/figure 2/191108 R8000 largest set of peptides medRD.pdf", plot = g.biggroup, scale = 0.6)


#################################################################################
############################## Figure 2D most R8000 orthogonal peptides per set
###########################################################################################
sumR8000maxs <- group_by(R8000maxs, name, X_group, classifier) %>%
  summarise(numorthopeps = n(), numorthoints = sum(classy))
sumsumR8000maxs <- group_by(sumR8000maxs, name, X_group) %>%
  summarise(maxints = max(numorthoints), maxnumpeps = max(numorthopeps))

g.sumR8000maxs <- ggplot(sumsumR8000maxs, aes(maxints, fill =  X_group)) +geom_histogram(bins = 11, show.legend = T, size = 1.5) + geom_vline(xintercept = 7, linetype = "dashed") + 
  labs(x = "Number of interactions", y = "Number of sets", title = "190402 R8000 Orthogonal sets interaction counts") + scale_x_continuous(breaks = c(2, 4 ,6, 8, 10, 12)) +
  scale_y_continuous(breaks = c(0, 2, 4, 6,8, 10)) + theme(strip.background = element_rect("white")) + scale_fill_manual(values = c("#FF0000","#F2AD00","#00A08A","#F98400", "#5BBCD6"), name = "Backbone") + 
  theme(axis.text = element_text(size = 14, family = "Myriad Web Pro"),text = element_text(size = 18, family = "Myriad Web Pro"))
g.sumR8000maxs


ggsave(g.sumR8000maxs, file = "../figures/figure 2/1901108 R8000 orthogonal set sizes.pdf",scale = 0.6)

colorCount <- length(unique(sumsumR8000maxs$name))
fastapal <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(17)


g.sumR0000sets <- ggplot(sumsumR8000maxs, aes(maxints, fill = name)) + geom_histogram(bins = 11) + scale_fill_manual(values = fastapal, labels = NULL, name = "Sets with\nsame interfaces") +
  theme(axis.text = element_text(size = 14, family = "Myriad Web Pro"),text = element_text(size = 18, family = "Myriad Web Pro")) +scale_y_continuous(breaks = c(0, 2, 4, 6,8, 10)) + 
  labs(x = "Number of interactions", y = "Number of sets", title = "200622 R8000 Orthogonal sets by backbones") + scale_x_continuous(breaks = c(2, 4 ,6, 8, 10, 12))
g.sumR0000sets

ggsave(g.sumR0000sets, file = "../figures/Figure S13 200622 R8000 orthogonal set sizes.pdf",scale = 0.6)

###########################################################################################
###################################### Figure 2E Interactions per orthogonal interaction
###########################################################################################
fastas <- list.files(path = "../designs/170317 R8000 all sets !OUT/", pattern = ".fasta", full.names = TRUE)
allsets <- data.frame("name" = character(),"group" = character())
for (i in 1:length(fastas)){
  names <- read.table(file = fastas[i])
  odds <- seq(1,nrow(names),2)
  newset = data.frame("name" = names[odds,], "group" = fastas[i])
  allsets <- bind_rows(allsets, newset)
}
allsets$group <- gsub("../designs/170317 R8000 all sets !OUT//","",allsets$group)
sumallsets <- group_by(allsets, group) %>%
  summarise(numtotalpeps = n(), numtotalints = sum(1:n()))

sumR8000maxs <- left_join(sumR8000maxs, sumallsets, by = c("name" = "group"))
sumsumR8000maxs <- group_by(sumR8000maxs, name, X_group) %>%
  summarise(maxints = max(numorthoints), maxnumpeps = max(numorthopeps), ratio = max(numorthopeps)/max(numtotalints), numints = median(numtotalints))

fit <- lm(sumsumR8000maxs$maxnumpeps ~ sumsumR8000maxs$numints)
fit2 <-  lm(sumsumR8000maxs$maxnumpeps ~ poly(sumsumR8000maxs$numints, 2)) 

ratcor <- data.frame(x = 350, y = 150, label = paste("italic(rho)==",round(cor(sumsumR8000maxs$numints, sumsumR8000maxs$maxnumpeps, use = "complete.obs", method = 'spearman'), digits = 2)))
g.ratio <- ggplot(sumsumR8000maxs, aes(numints, maxnumpeps)) + geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = F) + labs(x = "Number of possible interactions per set", y = "Number of total orthogonal interactions", title = "191108 R8000 Number of interactions vs ratio of orthogonal interactions") +
  theme(text = element_text(family = "Myriad Web Pro", size = 14), axis.title = element_text(size = 18)) + geom_point(alpha = 0.5, size =4, position=position_jitter(h=5,w=5)) + geom_text(data = ratcor, aes(x =x, y=y, label = label), parse = T, size = 5) +
  theme_cowplot()
g.ratio

ggsave(filename = "../figures/figure 2/210513 R8000 numtot vs numortho interactions2.pdf", g.ratio, scale = 0.6)

###################################################################################












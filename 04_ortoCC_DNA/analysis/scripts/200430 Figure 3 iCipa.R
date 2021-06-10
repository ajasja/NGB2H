library(MASS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(GGally)
library(extrafont)
library(viridis)
############################################################################
# Read in the data and mash it together
##########################################################################
setwd("../data")

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

medR8000 <- group_by(filter(R8000_map, DNA_1_counts > 15, DNA_2_counts > 15), X_peptide, Y_peptide) %>%
  summarise(medRD1 = median(RNA_1_counts/DNA_1_counts), medRD2 = median(RNA_2_counts/DNA_2_counts), medRDboth = median((RNA_1_counts + RNA_2_counts)/(DNA_1_counts + DNA_2_counts)), bcnum = n(), sumD1 = sum(DNA_1_counts), sumD2 = sum(DNA_2_counts), sumDboth = sum(DNA_1_counts + DNA_2_counts))
medR8000 <- separate(medR8000, X_peptide, into = c("X_peptide","X_group"), sep = "-") %>%
  separate(Y_peptide, into = c("Y_peptide","Y_group"), sep = "-")
medR8000 <-  filter(medR8000, X_group == "bA", Y_group == "bA")


R8bcipa <- read.csv("../model_data/190408_R8000_bcipa.scoring", header = F)
R8complete <- read.csv("../model_data/190408_R8000_complete.scoring", header = F)
R8fong <- read.csv("../model_data/190408_R8000_fong_svm.scoring", header = F)
R8rfe <- read.csv("../model_data/190408_R8000_rfe.scoring", header = F)
R8vinson <- read.csv("../model_data/190408_R8000_vinson_ce.scoring", header = F)

colnames(R8bcipa) <- c("X_peptide", "Y_peptide", "bcipa")
colnames(R8complete) <- c("X_peptide", "Y_peptide", "complete")
colnames(R8fong) <- c("X_peptide", "Y_peptide", "fong")
colnames(R8rfe) <- c("X_peptide", "Y_peptide", "rfe")
colnames(R8vinson) <- c("X_peptide", "Y_peptide", "vinson")

modelsR8 <- inner_join(R8bcipa, R8complete)
modelsR8 <- inner_join(modelsR8, R8fong)
modelsR8 <- inner_join(modelsR8, R8rfe)
modelsR8 <- inner_join(modelsR8, R8vinson)
modelsR8 <- separate(modelsR8, X_peptide, into = c("X_peptide", "X_group"), sep = "-") %>%
  separate(Y_peptide, into = c("Y_peptide", "Y_group"), sep  = "-")
modelsR8 <- filter(modelsR8, X_group =="bA", Y_group =="bA")


get_density <- function(x,y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


medR8000model <- inner_join(medR8000, modelsR8)
densities <- data.frame(medR8000model[,1:4],"bcipa" = NA, "complete" = NA, "fong" = NA, "rfe" = NA, "vinson" = NA)
densities$bcipa <- get_density(log(medR8000model$medRDboth), medR8000model$bcipa, n = 400)
densities$complete <- get_density(log(medR8000model$medRDboth), medR8000model$complete, n=400)
densities$fong <- get_density(log(medR8000model$medRDboth), medR8000model$fong, n=400)
densities$rfe <- get_density(log(medR8000model$medRDboth), medR8000model$rfe, n=400)
densities$vinson <- get_density(log(medR8000model$medRDboth), medR8000model$vinson, n=400)
densities <- mutate(densities, "bcipa" = (bcipa - min(bcipa))/(max(bcipa) - min(bcipa)), "complete" = (complete - min(complete))/(max(complete) - min(complete)), "fong" = (fong - min(fong))/(max(fong) - min(fong)),
                    "rfe" = (rfe - min(rfe))/(max(rfe) - min(rfe)), "vinson" = (vinson - min(vinson))/(max(vinson) - min(vinson)))
 

medR8000model <- dplyr::select(medR8000model, -medRD1, -medRD2, -bcnum,-sumD1, -sumD2, -sumDboth)

medR8000model <- gather(medR8000model, key = "func", value = "score", -X_peptide, -Y_peptide, -medRDboth, -X_group, -Y_group)
densities <- gather(densities, key = "func", value = "density", -X_peptide, -Y_peptide, -X_group, -Y_group)
medR8000model <- left_join(medR8000model, densities, by = c("X_peptide","X_group", "Y_peptide", "Y_group", "func"))


############################################################################
########################## Figure 5A correlation of old models and bA data
############################################################################
modelcor <- data.frame(y = c(-2.2,-2.2,-2.2, -2.2, -2.2), x=c(-15,-6,-15,-6,0),func =c("bcipa","complete","fong","rfe","vinson") ,label =c(paste("italic(R^2)==",round(cor(log(filter(medR8000model,func == "bcipa")$medRDboth), filter(medR8000model,func == "bcipa")$score)^2, digits = 3)),
                                                                           paste("italic(R^2)==",round(cor(log(filter(medR8000model,func == "complete")$medRDboth), filter(medR8000model,func == "complete")$score)^2, digits = 3)),
                                                                           paste("italic(R^2)==",round(cor(log(filter(medR8000model,func == "fong")$medRDboth), filter(medR8000model,func == "fong")$score)^2, digits = 3)),
                                                                           paste("italic(R^2)==",round(cor(log(filter(medR8000model,func == "rfe")$medRDboth), filter(medR8000model,func == "rfe")$score)^2, digits = 3)),
                                                                           paste("italic(R^2)==",round(cor(log(filter(medR8000model,func == "vinson")$medRDboth), filter(medR8000model,func == "vinson")$score)^2, digits = 3))))

pfunc <- c("bcipa" = "bCipa", "complete" = "Popatov", "fong" = "Fong", "rfe" = "Popatov\nLite", "vinson" = "Vinson")
g.modelR8000 <- ggplot(medR8000model, aes(score,log(medRDboth),  color = density)) + geom_point(alpha = 0.3)  + geom_smooth(method = 'lm', se = FALSE, color = 'black') +geom_text(data= modelcor, aes(x=x,y=y,label=label),parse = T, size = 5, color = 'black') + facet_wrap(~func, scales = "free", ncol = 5, labeller = labeller("func" = pfunc)) +
  labs(x = "Interaction score", y= "Score (AU)", title = "190723 R8000 median RNA/DNA vs Scoring functions") + theme(strip.background = element_rect(fill = "white"), text = element_text(family = "Myriad Web Pro", size = 14), axis.title = element_text(size = 18)) +
  scale_color_viridis(name = "Density") + scale_y_continuous(limits = c(-2.7, 2.7))
g.modelR8000 

ggsave(filename = "../figures/figure 3/191105 R8000 medRD vs models2.pdf", plot = g.modelR8000, width = 15, height = 5)


############################################################################
########################## Figure 5B Correlation of new scoring functions with bA data
############################################################################
load(file = "icipa_candidatescores.Rda")
allscores <- mutate(allscores, "lgmedRD" = log(medRDboth))
allscore <- dplyr::select(allscores, lgmedRD, contains("b_"), ends_with("b"))

allscore <- dplyr::select(allscores, medRDboth, contains("b_"), ends_with("b"))
allscore <- mutate(allscore, "lgmedRDboth" = log(medRDboth))
samplesallscore = data_frame("name" = character(), "iteration" = numeric(), "score" = numeric())
for (i in 1:100){
sampallscore <- sample_frac(allscore, 0.1, replace = F)
sampallscorecor <- cor(sampallscore)
sampallscoforbind <- data.frame("name" = rownames(sampallscorecor), "iteration" = i, "score" = (data.frame(sampallscorecor)$lgmedRDboth)^2)
samplesallscore <- bind_rows(samplesallscore, sampallscoforbind)
}
samplesallscore <- data.frame(samplesallscore, "D_type" = 0)
samplesallscore$D_type[grepl("D5",samplesallscore$name)] <- 1
samplesallscore$name <-factor(samplesallscore$name, levels = c("D0_b", "D5_b","D0_b_cv","D5_b_cv","D0_b_nc","D5_b_nc", "lgmedRDboth", "medRDboth"))
sigs <- pairwise.t.test(samplesallscore$score,samplesallscore$name, p.adjust.method = 'BY')
usedsigs  <- formatC(c(sigs$p.value[1,1], sigs$p.value[3,3], sigs$p.value[5,5], sigs$p.value[5,4], sigs$p.value[5,2]), format = "e", digits = 2) 
               
g.allscore <- ggplot(filter(samplesallscore, name != "medRDboth", name != "lgmedRDboth"), aes(name, score, color = as.factor(D_type))) + geom_boxplot() + scale_color_manual(values = c("#FF0000","#F2AD00"), name = "Heptad Shifting", labels = c("No","Yes")) + 
  geom_jitter(alpha = 0.5, width = 0.2) + labs(x = "Model", y = "R^2 of interation score", title = "191105 R8000 medRD correlation with models") + scale_x_discrete(labels = c("Basic Model","","+A-position stacking terms","","+A-position N-terminal terms","")) +
  theme(text = element_text(family = "Myriad Web Pro", size = 14), axis.title = element_text(size = 18)) + annotate("text", x=c("D0_b","D0_b_cv","D0_b_nc", "D5_b_nc", "D5_b_nc"), y = c(0.43, 0.43, 0.43, 0.47, 0.5), label = usedsigs)
g.allscore

ggsave(filename = "../figures/figure 3/191105 R8000 icipa models correlation.pdf", plot = g.allscore, scale = 0.6)


############################################################################
########################## Figure 5C Correlation of icipa scoring function with bA data
############################################################################

mylm <- lm(log(allscores.medRDboth)~allscores.D5_b_nc, data.frame(allscores$medRDboth, allscores$D5_b_nc))

allscores <- data.frame(allscores, density = 0)
allscores$density <- get_density(log(allscores$medRDboth), allscores$D5_b_nc, n = 100)
corrD5bnc <- data.frame(y=1, x=60, label=paste("italic(r^2)==", round(cor(log(allscores$medRDboth), allscores$D5_b_nc)^2, digits = 3)))
g.D5_b_nc <- ggplot(allscores, aes( D5_b_nc, log(medRDboth),  color = density)) + geom_point(alpha = 0.4, size = 2) + geom_smooth(method = "lm",formula = y~x, se = FALSE, color = 'black') + geom_text(data = corrD5bnc, aes(x=x,y=y,label=label),parse=TRUE, size = 5, color = 'black') +
  labs(x = "Interaction score", y="Predicted value of interaction (D5bnc)", title = "190726 R8000-bAs predicted values") + theme(text = element_text(family = "Myriad Web Pro", size = 14), axis.title = element_text(size = 18)) +
  scale_color_viridis(name = "Density") + scale_y_continuous(limits =c(-2.7, 2.7)) + scale_x_continuous(limits =c(15,70))
g.D5_b_nc

ggsave(filename = "../figures/figure 3/191105 R8000 predicted by D5Bnc.pdf", plot= g.D5_b_nc, scale = 0.6)

############################################################################
########################## Figure 5D weights for icipa scoring function
############################################################################
weights <- read.csv(file = "DNA-ALL-basic-rep-nter_core-Ridge-WbnRD10.features.csv")
weights <- filter(weights, N_iter ==  5)
weights <- separate(weights, feature, into = c("start", "residues"), sep = -3)
weights$residues <- factor(weights$residues, levels =  c("_NN","_II","_IN","_EK","_KK","_EE","_LL","ion"))

g.weights <- ggplot(weights, aes(residues, coef, fill = start)) + geom_bar(stat = "identity", position = "dodge") + scale_fill_manual(name = "Type",labels = c("Position A", "Electrostatic","Interface\nRepulsion","N-terminal\nPosition A"), values = c("#FF0000","#F2AD00", "#00A08A","#5BBCD6")) +
  scale_x_discrete(labels = c("_NN" = "NN","_II" = "II","_IN" = "IN", "_EK" = "EK","_KK" = "KK","_EE" = "EE","_LL" = "LL", "ion" ="Overall charge")) + labs(x = "Residues", y = "Weight", title = "190724 iCipa (DNA-ALL-basicL-rep-nter_core-Ridge-WbnRD10) Weights") + 
  theme(text = element_text(family = "Myriad Web Pro", size = 14), axis.title = element_text(size = 18))
g.weights

ggsave(filename = "../figures/figure 3/1901105 icipa weights.pdf", plot= g.weights, scale = 0.6)

############################################################################
########################## Figure 5E correlation of icipa and others with mason Tm data
############################################################################
###########################################################
mason_icipa <- read.csv(file = "mason_5_DNA-ALL-basicL-rep-nter_cor-Ridge_WbnRD10.csv")
mason_icipa <- separate(mason_icipa, ID1, into = c("X_peptide","X_group"), sep = "-") %>%
  separate(ID2, into = c("Y_peptide","Y_group"), sep = "-")
load("masons.Rda")
masons <- inner_join(mason_icipa, mason_tms, by = c("X_peptide", "Y_peptide"))
medR8000model2 <- inner_join(medR8000, modelsR8) 
medR8000model2 <- dplyr::select(medR8000model2, -X_group, -Y_group, -medRD1, -medRD2, -bcnum,-sumD1, -sumD2, -sumDboth)

masons <- inner_join(masons, medR8000model2)
masons <- mutate(masons, medRDboth = log(medRDboth))
masons <- dplyr::select(masons, -score_bc7, -RD_Tm, -Tm.x, -norm_sd_RD, -cv_RD)
thinmason <- gather(masons, key = "sfunction", value = "fscore", -X_peptide, -Y_peptide, -X_group, -Y_group,  -Tm.y)
icipa <- filter(thinmason, sfunction == "score")
popatov <- filter(thinmason, sfunction == "complete")
bcipa <- filter(thinmason, sfunction == "bcipa")
icipa <- mutate(icipa, score = (fscore - min(fscore))/(max(fscore) - min(fscore)))
popatov <- mutate(popatov, score = 1 - (fscore - min(fscore))/(max(fscore) - min(fscore)))
bcipa <- mutate(bcipa, score = 1 - (fscore - min(fscore))/(max(fscore) - min(fscore)))


corrmasonmall <- data.frame(x = c(0.25, .85, 0.05), y = c(65, 0, 45),sfunction = c("score","complete", "bcipa") ,label = c(paste("italic(r^2)==", round(cor(icipa$score, icipa$Tm.y)^2, digits = 3)), 
                                                                                                                           paste("italic(r^2)==", round(cor(popatov$score, popatov$Tm.y)^2, digits = 3)),
                                                                                                                           paste("italic(r^2)==", round(cor(bcipa$score, bcipa$Tm.y)^2, digits = 3))))
g.masonsmall <- ggplot(filter(thinmason, sfunction =="score" | sfunction == "complete" | sfunction == "bcipa"), aes(y =Tm.y, color = sfunction)) +
  geom_jitter(data = icipa, aes(x = score), alpha = 0.8, size = 4, width = 0.02) + geom_smooth(data = icipa, aes(x = score), method = "lm", se = F, show.legend = F) +
  geom_jitter(data = popatov, aes(x = score), alpha = 0.8, size = 4, width = 0.02) + geom_smooth(data = popatov, aes(x = score), method = "lm", se = F) + 
  geom_jitter(data = bcipa, aes(x = score), alpha = 0.8, size = 4, width = 0.02)  +geom_smooth(data = bcipa, aes(x = score), method = "lm", se = F) + 
  geom_text(data = corrmasonmall, aes(x=x,y=y,label=label), parse = TRUE, show.legend = FALSE, size = 5) +
  scale_color_manual(name = "Scoring\nfunction",labels = c("bCipa", "Popatov", "iCipa"), values = c("#FF0000","#F2AD00","#00A08A")) + labs(x = "Normalized score", y = "Previously reported Tm", title = "190729 scoring functions vs mason Tm") +
  theme(text = element_text(family = "Myriad Web Pro", size = 14), axis.title = element_text(size = 18))
g.masonsmall

ggsave(filename = "../figures/figure 3/191105 Mason vs scoring funcs.pdf", plot = g.masonsmall, scale = 0.6)










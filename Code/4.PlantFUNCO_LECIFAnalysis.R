###################################
#   TAIOTACTICAL-LECIF analyses   #
###################################

### 0. Load libs and pkgs

library(tidyverse)
library(hrbrthemes)
library(viridis)
library(ggplot2)
library(RColorBrewer)
library(philentropy)
library(proxyC)
library(cowplot)
library(patternplot)
library(LOLA)
library(GenomicRanges)
library(GenomicFeatures)
library(simpleCache)
library(circlize)

###################################################################
########## 1. LECIF score Chromatin States-based analysis ######### 
###################################################################

### 1. LECIF score over 16 chromatin states for 3 species and 6 LECIF comparisons

LECIF_AtOs <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/At-Os_LECIFscore_hiHMM.bedgraph", header = F, sep = "\t")
LECIF_AtZm <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea//output/At-Zm_LECIFscore_hiHMM.bedgraph", header = F, sep = "\t")
LECIF_OsAt <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/Os-At_LECIFscore_hiHMM.bedgraph", header = F, sep = "\t")
LECIF_OsZm <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/Os-Zm_LECIFscore_hiHMM.bedgraph", header = F, sep = "\t")
LECIF_ZmAt <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/Zm-At_LECIFscore_hiHMM.bedgraph", header = F, sep = "\t")
LECIF_ZmOs <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/Zm-Os_LECIFscore_hiHMM.bedgraph", header = F, sep = "\t")

### 1.2 Simmilarity of Chromatin States (hiHMM) per LECIF score percentile ranks bins
### We got only 1 epigenome because of inter-species chromatin states with multiple conditions
### So instead of frequency vector of chromatin states for each epigenome and pearson -->
### binary vectors with each element of the vectors constituying a region, 1 if that region is annotated by the Chromatin state, 0 if not
### Jaccard/SMC or other simmilarity metric between binary vectors so if the 1s are in the same position or less changes are needed to convert 
### 0s to 1s --> higher simmilarity

###############################################
########## 1.2.1 Arabidopsis - Oryza ######### 
##############################################

OsAt_regions <- read.delim("H:/Taiotactical/Taiotactical_LECIF/example/check_at-os/all.m", header = F, sep = "\t") 
AtOs_regions <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/At-Os_LECIFscore.bedgraph", header = F, sep = "\t") 
AtOs_bed <- AtOs_regions[,c(1,2,3)]
OsAt_bed <- OsAt_regions[,c(1,2,3)]
OsAt_bed <- OsAt_bed[grep("Chr", OsAt_bed$V1),]
OsAt_bed$V3 <- as.numeric(OsAt_bed$V3) - 1
OsAt_bed$V2 <- as.numeric(OsAt_bed$V2)
nrow(AtOs_bed) == nrow(OsAt_bed) # check the same number of at aligning regions to os
OsAt_bed$V3 <- OsAt_bed$V3 + (AtOs_bed$V3 - AtOs_bed$V2)
length(which(((AtOs_bed$V3 - AtOs_bed$V2) == (OsAt_bed$V3) -  (OsAt_bed$V2)) == TRUE)) # check same expansion of the first aligning base
write.table(OsAt_bed,"H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/Os_AtAligningRegions.bed", quote = F, sep = "\t", col.names = F, row.names = F)
# samtools sort -k1,1 -k2,2n Os_AtAligningRegions.bed > Os_AtAligningRegions.sorted.bed
# bedtools intersect -a Os_AtAligningRegions.sorted.bed -b hihmm.model1.K15.Os.ReMapped.ReNamed.Reduced.Reordered.bed -wa -wb | cut -f1,2,3,7 > Os_AtAligningRegions.hiHMM.sorted.bed
AtOs_hiHMMInfo <- LECIF_AtOs[,c(1,2,3,5)] 
OsAt_hiHMMInfo <- read.table("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/Os_AtAligningRegions.hiHMM.sorted.bed", header = F, sep = "\t")
length(percent_rank(x = AtOs_regions$V4)) == nrow(AtOs.df.F) # check
AtOs.df.F <- data.frame(Chr.At = AtOs_bed$V1, Start.At = AtOs_bed$V2, End.At = AtOs_bed$V3,
                        Chr.Os = OsAt_bed$V1, Start.Os = OsAt_bed$V2, End.Os = OsAt_bed$V3,
                        LECIF.Score = AtOs_regions$V4, Percentile.Rank.Score = percent_rank(x = AtOs_regions$V4), Percentile.Rank.Bin = NA,
                        CS1.At = 0, CS2.At = 0, CS3.At = 0, CS4.At = 0, CS5.At = 0, CS6.At = 0,
                        CS7.At = 0, CS8.At = 0, CS9.At = 0, CS10.At = 0, CS11.At = 0, CS12.At = 0,
                        CS13.At = 0, CS14.At = 0, CS15.At = 0, CS16.At = 0,
                        CS1.Os = 0, CS2.Os = 0, CS3.Os = 0, CS4.Os = 0, CS5.Os = 0, CS6.Os = 0,
                        CS7.Os = 0, CS8.Os = 0, CS9.Os = 0, CS10.Os = 0, CS11.Os = 0, CS12.Os = 0,
                        CS13.Os = 0, CS14.Os = 0, CS15.Os = 0, CS16.Os = 0,
                        Bivalent.At = 0, Active.At = 0, Divergent.At = 0, Heterochromatin.At = 0, Quies.At = 0,
                        Bivalent.Os = 0, Active.Os = 0, Divergent.Os = 0, Heterochromatin.Os = 0, Quies.Os = 0)

AtOs.df.F$Percentile.Rank.Bin[which(AtOs.df.F$Percentile.Rank.Score < 0.1)] <- "First" # 98927 regions
AtOs.df.F$Percentile.Rank.Bin[which(AtOs.df.F$Percentile.Rank.Score >= 0.1 & AtOs.df.F$Percentile.Rank.Score < 0.2)] <- "Second" # 98893
AtOs.df.F$Percentile.Rank.Bin[which(AtOs.df.F$Percentile.Rank.Score >= 0.2 & AtOs.df.F$Percentile.Rank.Score < 0.3)] <- "Thirdth" # 99076
AtOs.df.F$Percentile.Rank.Bin[which(AtOs.df.F$Percentile.Rank.Score >= 0.3 & AtOs.df.F$Percentile.Rank.Score < 0.4)] <- "Fourth" # 101338
AtOs.df.F$Percentile.Rank.Bin[which(AtOs.df.F$Percentile.Rank.Score >= 0.4 & AtOs.df.F$Percentile.Rank.Score < 0.5)] <- "Fifth" # 96302
AtOs.df.F$Percentile.Rank.Bin[which(AtOs.df.F$Percentile.Rank.Score >= 0.5 & AtOs.df.F$Percentile.Rank.Score < 0.6)] <- "Sixth" # 98879
AtOs.df.F$Percentile.Rank.Bin[which(AtOs.df.F$Percentile.Rank.Score >= 0.6 & AtOs.df.F$Percentile.Rank.Score < 0.7)] <- "Seventh" # 100699
AtOs.df.F$Percentile.Rank.Bin[which(AtOs.df.F$Percentile.Rank.Score >= 0.7 & AtOs.df.F$Percentile.Rank.Score < 0.8)] <- "Eigth" # 104676
AtOs.df.F$Percentile.Rank.Bin[which(AtOs.df.F$Percentile.Rank.Score >= 0.8 & AtOs.df.F$Percentile.Rank.Score < 0.9)] <- "Nineth" # 91322
AtOs.df.F$Percentile.Rank.Bin[which(AtOs.df.F$Percentile.Rank.Score >= 0.9 & AtOs.df.F$Percentile.Rank.Score < 1)] <- "Tenth" # 98900
(98927 + 98893 + 99076 + 101338 + 96302 + 98879 + 100699 + 104676 + 91322 + 98900) == nrow(AtOs.df.F) # check

for(i in 1:nrow(AtOs.df.F)){
  cat(paste0(i / nrow(AtOs.df.F) * 100, ' % completed \t'))
  AtCSinfo <- AtOs_hiHMMInfo$V5[which(AtOs_hiHMMInfo$V1 == AtOs.df.F$Chr.At[i] & AtOs_hiHMMInfo$V2 == AtOs.df.F$Start.At[i] & AtOs_hiHMMInfo$V3 == AtOs.df.F$End.At[i])]
  OsCSinfo <- OsAt_hiHMMInfo$V4[which(OsAt_hiHMMInfo$V1 == AtOs.df.F$Chr.Os[i] & OsAt_hiHMMInfo$V2 == AtOs.df.F$Start.Os[i] & OsAt_hiHMMInfo$V3 == AtOs.df.F$End.Os[i])]
  AtOs.df.F[i,which(c(paste0("CS",1:16)) %in% AtCSinfo)+9] <- 1
  AtOs.df.F[i,which(c(paste0("CS",1:16)) %in% OsCSinfo)+25] <- 1
  if(length(which(c("CS1", "CS2", "CS3", "CS4", "CS5") %in% AtCSinfo))>0){
    AtOs.df.F$Bivalent.At[i] <- 1
  }
  if(length(which(c("CS1", "CS2", "CS3", "CS4", "CS5") %in% OsCSinfo))>0){
    AtOs.df.F$Bivalent.Os[i] <- 1
  }
  if(length(which(c("CS6", "CS7", "CS8", "CS9") %in% AtCSinfo))>0){
    AtOs.df.F$Active.At[i] <- 1
  }
  if(length(which(c("CS6", "CS7", "CS8", "CS9") %in% OsCSinfo))>0){
    AtOs.df.F$Active.Os[i] <- 1
  }
  if(length(which(c("CS10") %in% AtCSinfo))>0){
    AtOs.df.F$Divergent.At[i] <- 1
  }
  if(length(which(c("CS10") %in% OsCSinfo))>0){
    AtOs.df.F$Divergent.Os[i] <- 1
  }
  if(length(which(c("CS11", "CS12", "CS13") %in% AtCSinfo))>0){
    AtOs.df.F$Heterochromatin.At[i] <- 1
  }
  if(length(which(c("CS11", "CS12", "CS13") %in% OsCSinfo))>0){
    AtOs.df.F$Heterochromatin.Os[i] <- 1
  }
  if(length(which(c("CS14", "CS15", "CS16") %in% AtCSinfo))>0){
    AtOs.df.F$Quies.At[i] <- 1
  }
  if(length(which(c("CS14", "CS15", "CS16") %in% OsCSinfo))>0){
    AtOs.df.F$Quies.Os[i] <- 1
  }
}

#write.table(AtOs.df.F, "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/AtOs_CSsimmilarity.DF.txt", col.names = T, row.names = F, sep = "\t", quote = F)
#AtOs.df.F. <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/AtOs_CSsimmilarity.DF.txt",header = T, row.names = F, sep = "\t", quote = F)
# loop to compute simmilarity between this huge binary vectos || metric performing well in boolean data but not too influences by super-sparsed data (low coverage from some states)
CS.Simmilarity.LECIF.AtOs <- data.frame(Bin = rep(c("First", "Second", "Thirdth", "Fourth", "Fifth", "Sixth", "Seventh", "Eigth", "Nineth", "Tenth"),46),
                                        Metric = c(rep(c("Jaccard"),10),rep(c("SMC"),10),rep(c("Sokal"),10), rep(c("Tanimoto"),10),
                                                   rep(c("Dice"),10), rep(c("Hanman"),10), rep(c("Ochiai"),10), rep(c("Sokal2"),10),
                                                   rep(c("Phi"),10), rep(c("S2"),10), 
                                                   rep(c("Eucledian"),10), rep(c("Sorensen"),10), rep(c("Canberra"),10), rep(c("Wavehedges"),10), 
                                                   rep(c("Tanimoto2"),10), rep(c("Fidelity"),10), rep(c("SquaredChord"),10), rep(c("SquaredChi"),10), 
                                                   rep(c("AdditiveSym"),10), rep(c("Topsoe"),10), rep(c("Kumarjhonson"),10), rep(c("Manhattan"),10),
                                                   rep(c("Lorentzian"),10), rep(c("Czekanowski"),10), rep(c("Hassebrook"),10),  rep(c("SquaredEucledian"),10),
                                                   rep(c("ProbSymm"),10), rep(c("Kullbackleiber"),10), rep(c("Jensenshannon"),10), rep(c("Minkowski"),10), 
                                                   rep(c("Soergel"),10), rep(c("Innerproduct"),10), rep(c("Hellinger"),10), rep(c("Divergence"),10), 
                                                   rep(c("Jeffreys"),10), rep(c("Jensendifference"),10), rep(c("Chebysev"),10), rep(c("Kulczynskid"),10), 
                                                   rep(c("Harmonicmean"),10), rep(c("Matusita"),10), rep(c("Neyman"),10), rep(c("Clark"),10), 
                                                   rep(c("Kdivergence"),10), rep(c("Taneja"),10), 
                                                   rep(c("Maximum"),10), rep(c("Hamming"),10)),
                                        CS1 = NA, CS2 = NA, CS3 = NA, CS4 = NA,
                                        CS5 = NA, CS6 = NA, CS7 = NA, CS8 = NA,
                                        CS9 = NA, CS10 = NA, CS11 = NA, CS12 = NA,
                                        CS13 = NA, CS14 = NA, CS15 = NA, CS16 = NA,
                                        Bivalent = NA, Active = NA, Divergent = NA, Heterochromatin = NA, Quies = NA)
## Reorder prev dataframe to enter the loop
AtOs.df.F <- AtOs.df.F[,c(1:25, 42:46, 26:41, 47:51)]
## Loop to fill the distance/simmilarity metrics between CS (states and groups) splited by LECIF percentile bins
for(i in 1:21){# 16 CS and 5 CS functional groups
  cat(paste0(i / 21 * 100, ' % CS groups started \t'))
  for(j in 1:10){ # lecif score percentile bins
    cat(paste0(j / 10 * 100, ' % LECIF bins started \t'))
    CS.Simmilarity.LECIF.AtOs[j,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 1)
    CS.Simmilarity.LECIF.AtOs[j+10,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 2) 
    CS.Simmilarity.LECIF.AtOs[j+20,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 3)  
    CS.Simmilarity.LECIF.AtOs[j+30,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 4) 
    CS.Simmilarity.LECIF.AtOs[j+40,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 5) 
    CS.Simmilarity.LECIF.AtOs[j+50,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 6) 
    CS.Simmilarity.LECIF.AtOs[j+60,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 7) 
    CS.Simmilarity.LECIF.AtOs[j+70,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 8) 
    CS.Simmilarity.LECIF.AtOs[j+80,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 9) 
    CS.Simmilarity.LECIF.AtOs[j+90,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 10) 
    CS.Simmilarity.LECIF.AtOs[j+100,2+i] <- philentropy::euclidean(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F) 
    CS.Simmilarity.LECIF.AtOs[j+110,2+i] <- philentropy::sorensen(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+120,2+i] <- philentropy::canberra(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+130,2+i] <- philentropy::wave_hedges(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+140,2+i] <- philentropy::tanimoto(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+150,2+i] <- philentropy::fidelity(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+160,2+i] <- philentropy::squared_chord(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+170,2+i] <- philentropy::squared_chi_sq(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+180,2+i] <- philentropy::additive_symm_chi_sq(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+190,2+i] <- philentropy::topsoe(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F,unit = "log2")
    CS.Simmilarity.LECIF.AtOs[j+200,2+i] <- philentropy::kumar_johnson(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.AtOs[j+210,2+i] <- philentropy::manhattan(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+220,2+i] <- philentropy::lorentzian(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.AtOs[j+230,2+i] <- philentropy::czekanowski(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+240,2+i] <- philentropy::kumar_hassebrook(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+250,2+i] <- philentropy::squared_euclidean(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+260,2+i] <- philentropy::prob_symm_chi_sq(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+270,2+i] <- philentropy::kullback_leibler_distance(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.AtOs[j+280,2+i] <- philentropy::jensen_shannon(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.AtOs[j+290,2+i] <- philentropy::minkowski(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, n = 3)
    CS.Simmilarity.LECIF.AtOs[j+300,2+i] <- philentropy::soergel(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+310,2+i] <- philentropy::inner_product(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+320,2+i] <- philentropy::hellinger(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+330,2+i] <- philentropy::divergence_sq(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+340,2+i] <- philentropy::jeffreys(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.AtOs[j+350,2+i] <- philentropy::jensen_difference(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.AtOs[j+360,2+i] <- philentropy::chebyshev(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+370,2+i] <- philentropy::kulczynski_d(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.AtOs[j+380,2+i] <- philentropy::harmonic_mean_dist(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+390,2+i] <- philentropy::matusita(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+400,2+i] <- philentropy::neyman_chi_sq(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.AtOs[j+410,2+i] <- philentropy::clark_sq(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+420,2+i] <- philentropy::k_divergence(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.AtOs[j+430,2+i] <- philentropy::taneja(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.AtOs[j+440,2+i] <- proxyC::dist(as.matrix(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = "maximum", margin = 2)[1,2]
    CS.Simmilarity.LECIF.AtOs[j+450,2+i] <- proxyC::dist(as.matrix(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = "hamming", margin = 2)[1,2]
    cat(paste0(j / 10 * 100, ' % LECIF bins completed \t'))
    }
  cat(paste0(i / 21 * 100, ' % CS groups completed \t'))
}

write.table(AtOs.df.F, "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/AtOs_CSsimmilarity.DF10.txt", col.names = T, row.names = F, sep = "\t", quote = F)
CS.Simmilarity.LECIF.AtOs.10 <- CS.Simmilarity.LECIF.AtOs
write.table(CS.Simmilarity.LECIF.AtOs, file = "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/AtOs_CSsimmilarity.DF.Distances10.txt", col.names = T, row.names = F, sep = "\t", quote = F)

## Example code to check if metrics are distance or simmilarity and max/min values for each
ade4::dist.binary(t(data.frame(CS1.At = AtOs.df.F3$CS2.At[which(AtOs.df.F$Percentile.Rank.Bin == "First")], CS1.Os = AtOs.df.F3$CS2.At[which(AtOs.df.F$Percentile.Rank.Bin == "First")])), method = 5)
ade4::dist.binary(t(data.frame(CS1.At = AtOs.df.F$CS2.At[which(AtOs.df.F$Percentile.Rank.Bin == "Second")], CS1.Os = AtOs.df.F$CS2.Os[which(AtOs.df.F$Percentile.Rank.Bin == "Second")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.At = AtOs.df.F$CS2.At[which(AtOs.df.F$Percentile.Rank.Bin == "Thirdth")], CS1.Os = AtOs.df.F$CS2.Os[which(AtOs.df.F$Percentile.Rank.Bin == "Thirdth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.At = AtOs.df.F$CS2.At[which(AtOs.df.F$Percentile.Rank.Bin == "Fourth")], CS1.Os = AtOs.df.F$CS2.Os[which(AtOs.df.F$Percentile.Rank.Bin == "Fourth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.At = AtOs.df.F$CS2.At[which(AtOs.df.F$Percentile.Rank.Bin == "Fifth")], CS1.Os = AtOs.df.F$CS2.Os[which(AtOs.df.F$Percentile.Rank.Bin == "Fifth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.At = AtOs.df.F$CS2.At[which(AtOs.df.F$Percentile.Rank.Bin == "Sixth")], CS1.Os = AtOs.df.F$CS2.Os[which(AtOs.df.F$Percentile.Rank.Bin == "Sixth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.At = AtOs.df.F$CS2.At[which(AtOs.df.F$Percentile.Rank.Bin == "Seventh")], CS1.Os = AtOs.df.F$CS2.Os[which(AtOs.df.F$Percentile.Rank.Bin == "Seventh")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.At = AtOs.df.F$CS2.At[which(AtOs.df.F$Percentile.Rank.Bin == "Eigth")], CS1.Os = AtOs.df.F$CS2.Os[which(AtOs.df.F$Percentile.Rank.Bin == "Eigth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.At = AtOs.df.F$CS2.At[which(AtOs.df.F$Percentile.Rank.Bin == "Nineth")], CS1.Os = AtOs.df.F$CS2.Os[which(AtOs.df.F$Percentile.Rank.Bin == "Nineth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.At = AtOs.df.F$CS2.At[which(AtOs.df.F$Percentile.Rank.Bin == "Tenth")], CS1.Os = AtOs.df.F$CS2.Os[which(AtOs.df.F$Percentile.Rank.Bin == "Tenth")])), method = 2)
philentropy::jaccard(P =  AtOs.df.F3$CS1.At[which(AtOs.df.F$Percentile.Rank.Bin == "First")], Q = AtOs.df.F3$CS1.At[which(AtOs.df.F$Percentile.Rank.Bin == "First")], testNA = F)
proxyC::dist(as.matrix(data.frame(Gene1 = AtOs.df.F$CS1.At[which(AtOs.df.F$Percentile.Rank.Bin == "First")], Gene2 = AtOs.df.F$CS1.At[which(AtOs.df.F$Percentile.Rank.Bin == "First")])), method = "hamming", margin = 2)[1,2]

## binned in 5 sets

AtOs.df.F$Percentile.Rank.Bin.Five <- NA
AtOs.df.F$Percentile.Rank.Bin.Five[which(AtOs.df.F$Percentile.Rank.Score < 0.2)] <- "First" # 197820 regions
AtOs.df.F$Percentile.Rank.Bin.Five[which(AtOs.df.F$Percentile.Rank.Score >= 0.2 & AtOs.df.F$Percentile.Rank.Score < 0.4)] <- "Second" # 200414
AtOs.df.F$Percentile.Rank.Bin.Five[which(AtOs.df.F$Percentile.Rank.Score >= 0.4 & AtOs.df.F$Percentile.Rank.Score < 0.6)] <- "Thirdth" # 195181
AtOs.df.F$Percentile.Rank.Bin.Five[which(AtOs.df.F$Percentile.Rank.Score >= 0.6 & AtOs.df.F$Percentile.Rank.Score < 0.8)] <- "Fourth" # 205375
AtOs.df.F$Percentile.Rank.Bin.Five[which(AtOs.df.F$Percentile.Rank.Score >= 0.8 & AtOs.df.F$Percentile.Rank.Score <= 1)] <- "Fifth" # 190222
(197820 + 200414 + 195181 + 205375 + 190222) == nrow(AtOs.df.F) # check


CS.Simmilarity.LECIF.AtOs <- data.frame(Bin = rep(c("First", "Second", "Thirdth", "Fourth", "Fifth"),46),
                                        Metric = c(rep(c("Jaccard"),5),rep(c("SMC"),5),rep(c("Sokal"),5), rep(c("Tanimoto"),5),
                                                   rep(c("Dice"),5), rep(c("Hanman"),5), rep(c("Ochiai"),5), rep(c("Sokal2"),5),
                                                   rep(c("Phi"),5), rep(c("S2"),5), 
                                                   rep(c("Eucledian"),5), rep(c("Sorensen"),5), rep(c("Canberra"),5), rep(c("Wavehedges"),5), 
                                                   rep(c("Tanimoto2"),5), rep(c("Fidelity"),5), rep(c("SquaredChord"),5), rep(c("SquaredChi"),5), 
                                                   rep(c("AdditiveSym"),5), rep(c("Topsoe"),5), rep(c("Kumarjhonson"),5), rep(c("Manhattan"),5),
                                                   rep(c("Lorentzian"),5), rep(c("Czekanowski"),5), rep(c("Hassebrook"),5),  rep(c("SquaredEucledian"),5),
                                                   rep(c("ProbSymm"),5), rep(c("Kullbackleiber"),5), rep(c("Jensenshannon"),5), rep(c("Minkowski"),5), 
                                                   rep(c("Soergel"),5), rep(c("Innerproduct"),5), rep(c("Hellinger"),5), rep(c("Divergence"),5), 
                                                   rep(c("Jeffreys"),5), rep(c("Jensendifference"),5), rep(c("Chebysev"),5), rep(c("Kulczynskid"),5), 
                                                   rep(c("Harmonicmean"),5), rep(c("Matusita"),5), rep(c("Neyman"),5), rep(c("Clark"),5), 
                                                   rep(c("Kdivergence"),5), rep(c("Taneja"),5), 
                                                   rep(c("Maximum"),5), rep(c("Hamming"),5)),
                                                   CS1 = NA, CS2 = NA, CS3 = NA, CS4 = NA,
                                                   CS5 = NA, CS6 = NA, CS7 = NA, CS8 = NA,
                                                   CS9 = NA, CS10 = NA, CS11 = NA, CS12 = NA,
                                                   CS13 = NA, CS14 = NA, CS15 = NA, CS16 = NA,
                                                   Bivalent = NA, Active = NA, Divergent = NA, Heterochromatin = NA, Quies = NA)

for(i in 1:21){# 16 CS and 5 CS functional groups
  cat(paste0(i / 21 * 100, ' % CS groups started \t'))
  for(j in 1:5){ # lecif score percentile bins
    cat(paste0(j / 5 * 100, ' % LECIF bins started \t'))
    CS.Simmilarity.LECIF.AtOs[j,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 1)
    CS.Simmilarity.LECIF.AtOs[j+5,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 2) 
    CS.Simmilarity.LECIF.AtOs[j+10,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 3)  
    CS.Simmilarity.LECIF.AtOs[j+15,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 4) 
    CS.Simmilarity.LECIF.AtOs[j+20,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 5) 
    CS.Simmilarity.LECIF.AtOs[j+25,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 6) 
    CS.Simmilarity.LECIF.AtOs[j+30,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 7) 
    CS.Simmilarity.LECIF.AtOs[j+35,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 8) 
    CS.Simmilarity.LECIF.AtOs[j+40,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 9) 
    CS.Simmilarity.LECIF.AtOs[j+45,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 10) 
    CS.Simmilarity.LECIF.AtOs[j+50,2+i] <- philentropy::euclidean(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F) 
    CS.Simmilarity.LECIF.AtOs[j+55,2+i] <- philentropy::sorensen(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+60,2+i] <- philentropy::canberra(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+65,2+i] <- philentropy::wave_hedges(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+70,2+i] <- philentropy::tanimoto(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+75,2+i] <- philentropy::fidelity(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+80,2+i] <- philentropy::squared_chord(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+85,2+i] <- philentropy::squared_chi_sq(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+90,2+i] <- philentropy::additive_symm_chi_sq(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+95,2+i] <- philentropy::topsoe(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F,unit = "log2")
    CS.Simmilarity.LECIF.AtOs[j+100,2+i] <- philentropy::kumar_johnson(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.AtOs[j+105,2+i] <- philentropy::manhattan(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+110,2+i] <- philentropy::lorentzian(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.AtOs[j+115,2+i] <- philentropy::czekanowski(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+120,2+i] <- philentropy::kumar_hassebrook(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+125,2+i] <- philentropy::squared_euclidean(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+130,2+i] <- philentropy::prob_symm_chi_sq(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+135,2+i] <- philentropy::kullback_leibler_distance(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.AtOs[j+140,2+i] <- philentropy::jensen_shannon(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.AtOs[j+145,2+i] <- philentropy::minkowski(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, n = 3)
    CS.Simmilarity.LECIF.AtOs[j+150,2+i] <- philentropy::soergel(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+155,2+i] <- philentropy::inner_product(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+160,2+i] <- philentropy::hellinger(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+165,2+i] <- philentropy::divergence_sq(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+170,2+i] <- philentropy::jeffreys(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.AtOs[j+175,2+i] <- philentropy::jensen_difference(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.AtOs[j+180,2+i] <- philentropy::chebyshev(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+185,2+i] <- philentropy::kulczynski_d(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.AtOs[j+190,2+i] <- philentropy::harmonic_mean_dist(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+195,2+i] <- philentropy::matusita(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+200,2+i] <- philentropy::neyman_chi_sq(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.AtOs[j+205,2+i] <- philentropy::clark_sq(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+210,2+i] <- philentropy::k_divergence(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.AtOs[j+215,2+i] <- philentropy::taneja(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.AtOs[j+220,2+i] <- proxyC::dist(as.matrix(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = "maximum", margin = 2)[1,2]
    CS.Simmilarity.LECIF.AtOs[j+225,2+i] <- proxyC::dist(as.matrix(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = "hamming", margin = 2)[1,2]
    cat(paste0(j / 5 * 100, ' % LECIF bins completed \t'))
  }
  cat(paste0(i / 21 * 100, ' % CS groups completed \t'))
}

write.table(AtOs.df.F, "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/AtOs_CSsimmilarity.DF5.txt", col.names = T, row.names = F, sep = "\t", quote = F)
CS.Simmilarity.LECIF.AtOs.5 <- CS.Simmilarity.LECIF.AtOs
write.table(CS.Simmilarity.LECIF.AtOs, file = "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/AtOs_CSsimmilarity.DF.Distances5.txt", col.names = T, row.names = F, sep = "\t", quote = F)
View(CS.Simmilarity.LECIF.AtOs.5)

AtOs.df.F$Percentile.Rank.Bin.Three <- NA
AtOs.df.F$Percentile.Rank.Bin.Three[which(AtOs.df.F$Percentile.Rank.Score < 0.333333)] <- "First" # 329870 regions
AtOs.df.F$Percentile.Rank.Bin.Three[which(AtOs.df.F$Percentile.Rank.Score >= 0.333333 & AtOs.df.F$Percentile.Rank.Score < 0.666666)] <- "Second" # 329497
AtOs.df.F$Percentile.Rank.Bin.Three[which(AtOs.df.F$Percentile.Rank.Score >= 0.666666 & AtOs.df.F$Percentile.Rank.Score <= 1)] <- "Thirdth" # 329645
( 329870 + 329497 + 329645) == nrow(AtOs.df.F) # check


CS.Simmilarity.LECIF.AtOs <- data.frame(Bin = rep(c("First", "Second", "Thirdth"),46),
                                        Metric = c(rep(c("Jaccard"),3),rep(c("SMC"),3),rep(c("Sokal"),3), rep(c("Tanimoto"),3),
                                                   rep(c("Dice"),3), rep(c("Hanman"),3), rep(c("Ochiai"),3), rep(c("Sokal2"),3),
                                                   rep(c("Phi"),3), rep(c("S2"),3), 
                                                   rep(c("Eucledian"),3), rep(c("Sorensen"),3), rep(c("Canberra"),3), rep(c("Wavehedges"),3), 
                                                   rep(c("Tanimoto2"),3), rep(c("Fidelity"),3), rep(c("SquaredChord"),3), rep(c("SquaredChi"),3), 
                                                   rep(c("AdditiveSym"),3), rep(c("Topsoe"),3), rep(c("Kumarjhonson"),3), rep(c("Manhattan"),3),
                                                   rep(c("Lorentzian"),3), rep(c("Czekanowski"),3), rep(c("Hassebrook"),3),  rep(c("SquaredEucledian"),3),
                                                   rep(c("ProbSymm"),3), rep(c("Kullbackleiber"),3), rep(c("Jensenshannon"),3), rep(c("Minkowski"),3), 
                                                   rep(c("Soergel"),3), rep(c("Innerproduct"),3), rep(c("Hellinger"),3), rep(c("Divergence"),3), 
                                                   rep(c("Jeffreys"),3), rep(c("Jensendifference"),3), rep(c("Chebysev"),3), rep(c("Kulczynskid"),3), 
                                                   rep(c("Harmonicmean"),3), rep(c("Matusita"),3), rep(c("Neyman"),3), rep(c("Clark"),3), 
                                                   rep(c("Kdivergence"),3), rep(c("Taneja"),3), 
                                                   rep(c("Maximum"),3), rep(c("Hamming"),3)),
                                        CS1 = NA, CS2 = NA, CS3 = NA, CS4 = NA,
                                        CS5 = NA, CS6 = NA, CS7 = NA, CS8 = NA,
                                        CS9 = NA, CS10 = NA, CS11 = NA, CS12 = NA,
                                        CS13 = NA, CS14 = NA, CS15 = NA, CS16 = NA,
                                        Bivalent = NA, Active = NA, Divergent = NA, Heterochromatin = NA, Quies = NA)

for(i in 1:21){# 16 CS and 5 CS functional groups
  cat(paste0(i / 21 * 100, ' % CS groups started \t'))
  for(j in 1:3){ # lecif score percentile bins
    cat(paste0(j / 3 * 100, ' % LECIF bins started \t'))
    CS.Simmilarity.LECIF.AtOs[j,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 1)
    CS.Simmilarity.LECIF.AtOs[j+3,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 2) 
    CS.Simmilarity.LECIF.AtOs[j+6,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 3)  
    CS.Simmilarity.LECIF.AtOs[j+9,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 4) 
    CS.Simmilarity.LECIF.AtOs[j+12,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 5) 
    CS.Simmilarity.LECIF.AtOs[j+15,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 6) 
    CS.Simmilarity.LECIF.AtOs[j+18,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 7) 
    CS.Simmilarity.LECIF.AtOs[j+21,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 8) 
    CS.Simmilarity.LECIF.AtOs[j+24,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 9) 
    CS.Simmilarity.LECIF.AtOs[j+27,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = 10) 
    CS.Simmilarity.LECIF.AtOs[j+30,2+i] <- philentropy::euclidean(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F) 
    CS.Simmilarity.LECIF.AtOs[j+33,2+i] <- philentropy::sorensen(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+36,2+i] <- philentropy::canberra(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+39,2+i] <- philentropy::wave_hedges(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+42,2+i] <- philentropy::tanimoto(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+45,2+i] <- philentropy::fidelity(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+48,2+i] <- philentropy::squared_chord(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+51,2+i] <- philentropy::squared_chi_sq(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+54,2+i] <- philentropy::additive_symm_chi_sq(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+57,2+i] <- philentropy::topsoe(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F,unit = "log2")
    CS.Simmilarity.LECIF.AtOs[j+60,2+i] <- philentropy::kumar_johnson(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.AtOs[j+63,2+i] <- philentropy::manhattan(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+66,2+i] <- philentropy::lorentzian(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.AtOs[j+69,2+i] <- philentropy::czekanowski(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+72,2+i] <- philentropy::kumar_hassebrook(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+75,2+i] <- philentropy::squared_euclidean(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+78,2+i] <- philentropy::prob_symm_chi_sq(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+81,2+i] <- philentropy::kullback_leibler_distance(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.AtOs[j+84,2+i] <- philentropy::jensen_shannon(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.AtOs[j+87,2+i] <- philentropy::minkowski(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, n = 3)
    CS.Simmilarity.LECIF.AtOs[j+90,2+i] <- philentropy::soergel(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+93,2+i] <- philentropy::inner_product(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+96,2+i] <- philentropy::hellinger(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+99,2+i] <- philentropy::divergence_sq(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+102,2+i] <- philentropy::jeffreys(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.AtOs[j+105,2+i] <- philentropy::jensen_difference(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.AtOs[j+108,2+i] <- philentropy::chebyshev(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+111,2+i] <- philentropy::kulczynski_d(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.AtOs[j+114,2+i] <- philentropy::harmonic_mean_dist(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+117,2+i] <- philentropy::matusita(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+120,2+i] <- philentropy::neyman_chi_sq(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.AtOs[j+123,2+i] <- philentropy::clark_sq(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtOs[j+126,2+i] <- philentropy::k_divergence(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.AtOs[j+129,2+i] <- philentropy::taneja(P = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], Q = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.AtOs[j+132,2+i] <- proxyC::dist(as.matrix(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = "maximum", margin = 2)[1,2]
    CS.Simmilarity.LECIF.AtOs[j+135,2+i] <- proxyC::dist(as.matrix(data.frame(CS.At = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),9+i], CS.Os = AtOs.df.F[which(AtOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtOs[j,1]),30+i])), method = "hamming", margin = 2)[1,2]
    cat(paste0(j / 3 * 100, ' % LECIF bins completed \t'))
  }
  cat(paste0(i / 21 * 100, ' % CS groups completed \t'))
}

write.table(AtOs.df.F, "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/AtOs_CSsimmilarity.DF3.txt", col.names = T, row.names = F, sep = "\t", quote = F)
CS.Simmilarity.LECIF.AtOs.3 <- CS.Simmilarity.LECIF.AtOs
write.table(CS.Simmilarity.LECIF.AtOs, file = "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/AtOs_CSsimmilarity.DF.Distances3.txt", col.names = T, row.names = F, sep = "\t", quote = F)
View(CS.Simmilarity.LECIF.AtOs.3) ## Dice-like would be the best option

###############################################
########## 1.2.2 Arabidopsis - Zea ########### 
##############################################

ZmAt_regions <- read.delim("H:/Taiotactical/Taiotactical_LECIF/example/check_at-zm//all.m", header = F, sep = "\t") 
AtZm_regions <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/At-Zm_LECIFscore.bedgraph", header = F, sep = "\t") 
AtZm_bed <- AtZm_regions[,c(1,2,3)]
ZmAt_bed <- ZmAt_regions[,c(1,2,3)]
ZmAt_bed <- ZmAt_bed[which(ZmAt_bed$V1 %in% c(1:10)),]
ZmAt_bed$V3 <- as.numeric(ZmAt_bed$V3) - 1
ZmAt_bed$V2 <- as.numeric(ZmAt_bed$V2)
nrow(AtZm_bed) == nrow(ZmAt_bed) # check the same number of at aligning regions to Zm
ZmAt_bed$V3 <- ZmAt_bed$V3 + (AtZm_bed$V3 - AtZm_bed$V2)
length(which(((AtZm_bed$V3 - AtZm_bed$V2) == (ZmAt_bed$V3) -  (ZmAt_bed$V2)) == TRUE)) # check same expansion of the first aligning base
write.table(ZmAt_bed,"H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/Zm_AtAligningRegions.bed", quote = F, sep = "\t", col.names = F, row.names = F)
# samtools sort -k1,1 -k2,2n Zm_AtAligningRegions.bed > Zm_AtAligningRegions.sorted.bed
# bedtools intersect -a Zm_AtAligningRegions.sorted.bed -b hihmm.model1.K15.Zm.ReMapped.ReNamed.Reduced.Reordered.bed -wa -wb | cut -f1,2,3,7 > Os_AtAligningRegions.hiHMM.sorted.bed
AtZm_hiHMMInfo <- LECIF_AtZm[,c(1,2,3,5)] 
ZmAt_hiHMMInfo <- read.table("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/Zm_AtAligningRegions.hiHMM.sorted.bed", header = F, sep = "\t")
AtZm.df.F <- data.frame(Chr.At = AtZm_bed$V1, Start.At = AtZm_bed$V2, End.At = AtZm_bed$V3,
                        Chr.Zm = ZmAt_bed$V1, Start.Zm = ZmAt_bed$V2, End.Zm = ZmAt_bed$V3,
                        LECIF.Score = AtZm_regions$V4, Percentile.Rank.Score = percent_rank(x = AtZm_regions$V4), Percentile.Rank.Bin = NA,
                        CS1.At = 0, CS2.At = 0, CS3.At = 0, CS4.At = 0, CS5.At = 0, CS6.At = 0,
                        CS7.At = 0, CS8.At = 0, CS9.At = 0, CS10.At = 0, CS11.At = 0, CS12.At = 0,
                        CS13.At = 0, CS14.At = 0, CS15.At = 0, CS16.At = 0,
                        CS1.Zm = 0, CS2.Zm = 0, CS3.Zm = 0, CS4.Zm = 0, CS5.Zm = 0, CS6.Zm = 0,
                        CS7.Zm = 0, CS8.Zm = 0, CS9.Zm = 0, CS10.Zm = 0, CS11.Zm = 0, CS12.Zm = 0,
                        CS13.Zm = 0, CS14.Zm = 0, CS15.Zm = 0, CS16.Zm = 0,
                        Bivalent.At = 0, Active.At = 0, Divergent.At = 0, Heterochromatin.At = 0, Quies.At = 0,
                        Bivalent.Zm = 0, Active.Zm = 0, Divergent.Zm = 0, Heterochromatin.Zm = 0, Quies.Zm = 0)
length(percent_rank(x = AtZm_regions$V4)) == nrow(AtZm.df.F) # check

AtZm.df.F$Percentile.Rank.Bin[which(AtZm.df.F$Percentile.Rank.Score < 0.1)] <- "First" # 82567 regions
AtZm.df.F$Percentile.Rank.Bin[which(AtZm.df.F$Percentile.Rank.Score >= 0.1 & AtZm.df.F$Percentile.Rank.Score < 0.2)] <- "Second" # 82538
AtZm.df.F$Percentile.Rank.Bin[which(AtZm.df.F$Percentile.Rank.Score >= 0.2 & AtZm.df.F$Percentile.Rank.Score < 0.3)] <- "Thirdth" # 82546
AtZm.df.F$Percentile.Rank.Bin[which(AtZm.df.F$Percentile.Rank.Score >= 0.3 & AtZm.df.F$Percentile.Rank.Score < 0.4)] <- "Fourth" # 82561
AtZm.df.F$Percentile.Rank.Bin[which(AtZm.df.F$Percentile.Rank.Score >= 0.4 & AtZm.df.F$Percentile.Rank.Score < 0.5)] <- "Fifth" # 82539
AtZm.df.F$Percentile.Rank.Bin[which(AtZm.df.F$Percentile.Rank.Score >= 0.5 & AtZm.df.F$Percentile.Rank.Score < 0.6)] <- "Sixth" # 82565
AtZm.df.F$Percentile.Rank.Bin[which(AtZm.df.F$Percentile.Rank.Score >= 0.6 & AtZm.df.F$Percentile.Rank.Score < 0.7)] <- "Seventh" # 82536
AtZm.df.F$Percentile.Rank.Bin[which(AtZm.df.F$Percentile.Rank.Score >= 0.7 & AtZm.df.F$Percentile.Rank.Score < 0.8)] <- "Eigth" # 82528
AtZm.df.F$Percentile.Rank.Bin[which(AtZm.df.F$Percentile.Rank.Score >= 0.8 & AtZm.df.F$Percentile.Rank.Score < 0.9)] <- "Nineth" # 82548
AtZm.df.F$Percentile.Rank.Bin[which(AtZm.df.F$Percentile.Rank.Score >= 0.9 & AtZm.df.F$Percentile.Rank.Score <= 1)] <- "Tenth" # 82546
(82567 + 82538 + 82546 + 82561 + 82539 + 82565 + 82536 + 82528 + 82548 + 82547) == nrow(AtZm.df.F) # check

for(i in 1:nrow(AtZm.df.F)){
  cat(paste0(i / nrow(AtZm.df.F) * 100, ' % completed \n'))
  AtCSinfo <- AtZm_hiHMMInfo$V5[which(AtZm_hiHMMInfo$V1 == AtZm.df.F$Chr.At[i] & AtZm_hiHMMInfo$V2 == AtZm.df.F$Start.At[i] & AtZm_hiHMMInfo$V3 == AtZm.df.F$End.At[i])]
  ZmCSinfo <- ZmAt_hiHMMInfo$V4[which(ZmAt_hiHMMInfo$V1 == AtZm.df.F$Chr.Zm[i] & ZmAt_hiHMMInfo$V2 == AtZm.df.F$Start.Zm[i] & ZmAt_hiHMMInfo$V3 == AtZm.df.F$End.Zm[i])]
  AtZm.df.F[i,which(c(paste0("CS",1:16)) %in% AtCSinfo)+9] <- 1
  AtZm.df.F[i,which(c(paste0("CS",1:16)) %in% ZmCSinfo)+25] <- 1
  if(length(which(c("CS1", "CS2", "CS3", "CS4", "CS5") %in% AtCSinfo))>0){
    AtZm.df.F$Bivalent.At[i] <- 1
  }
  if(length(which(c("CS1", "CS2", "CS3", "CS4", "CS5") %in% ZmCSinfo))>0){
    AtZm.df.F$Bivalent.Zm[i] <- 1
  }
  if(length(which(c("CS6", "CS7", "CS8", "CS9") %in% AtCSinfo))>0){
    AtZm.df.F$Active.At[i] <- 1
  }
  if(length(which(c("CS6", "CS7", "CS8", "CS9") %in% ZmCSinfo))>0){
    AtZm.df.F$Active.Zm[i] <- 1
  }
  if(length(which(c("CS10") %in% AtCSinfo))>0){
    AtZm.df.F$Divergent.At[i] <- 1
  }
  if(length(which(c("CS10") %in% ZmCSinfo))>0){
    AtZm.df.F$Divergent.Zm[i] <- 1
  }
  if(length(which(c("CS11", "CS12", "CS13") %in% AtCSinfo))>0){
    AtZm.df.F$Heterochromatin.At[i] <- 1
  }
  if(length(which(c("CS11", "CS12", "CS13") %in% ZmCSinfo))>0){
    AtZm.df.F$Heterochromatin.Zm[i] <- 1
  }
  if(length(which(c("CS14", "CS15", "CS16") %in% AtCSinfo))>0){
    AtZm.df.F$Quies.At[i] <- 1
  }
  if(length(which(c("CS14", "CS15", "CS16") %in% ZmCSinfo))>0){
    AtZm.df.F$Quies.Zm[i] <- 1
  }
}

#write.table(AtZm.df.F, "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/AtZm_CSsimmilarity.DF.txt", col.names = T, row.names = F, sep = "\t", quote = F)
#AtZm.df.F. <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/AtZm_CSsimmilarity.DF.txt",header = T, row.names = F, sep = "\t", quote = F)
# loop to compute simmilarity between this huge binary vectos || metric performing well in boolean data but not too influences by super-sparsed data (low coverage from some states depending on the socre bin)
CS.Simmilarity.LECIF.AtZm <- data.frame(Bin = rep(c("First", "Second", "Thirdth", "Fourth", "Fifth", "Sixth", "Seventh", "Eigth", "Nineth", "Tenth"),46),
                                        Metric = c(rep(c("Jaccard"),10),rep(c("SMC"),10),rep(c("Sokal"),10), rep(c("Tanimoto"),10),
                                                   rep(c("Dice"),10), rep(c("Hanman"),10), rep(c("Ochiai"),10), rep(c("Sokal2"),10),
                                                   rep(c("Phi"),10), rep(c("S2"),10), 
                                                   rep(c("Eucledian"),10), rep(c("Sorensen"),10), rep(c("Canberra"),10), rep(c("Wavehedges"),10), 
                                                   rep(c("Tanimoto2"),10), rep(c("Fidelity"),10), rep(c("SquaredChord"),10), rep(c("SquaredChi"),10), 
                                                   rep(c("AdditiveSym"),10), rep(c("Topsoe"),10), rep(c("Kumarjhonson"),10), rep(c("Manhattan"),10),
                                                   rep(c("Lorentzian"),10), rep(c("Czekanowski"),10), rep(c("Hassebrook"),10),  rep(c("SquaredEucledian"),10),
                                                   rep(c("ProbSymm"),10), rep(c("Kullbackleiber"),10), rep(c("Jensenshannon"),10), rep(c("Minkowski"),10), 
                                                   rep(c("Soergel"),10), rep(c("Innerproduct"),10), rep(c("Hellinger"),10), rep(c("Divergence"),10), 
                                                   rep(c("Jeffreys"),10), rep(c("Jensendifference"),10), rep(c("Chebysev"),10), rep(c("Kulczynskid"),10), 
                                                   rep(c("Harmonicmean"),10), rep(c("Matusita"),10), rep(c("Neyman"),10), rep(c("Clark"),10), 
                                                   rep(c("Kdivergence"),10), rep(c("Taneja"),10), 
                                                   rep(c("Maximum"),10), rep(c("Hamming"),10)),
                                        CS1 = NA, CS2 = NA, CS3 = NA, CS4 = NA,
                                        CS5 = NA, CS6 = NA, CS7 = NA, CS8 = NA,
                                        CS9 = NA, CS10 = NA, CS11 = NA, CS12 = NA,
                                        CS13 = NA, CS14 = NA, CS15 = NA, CS16 = NA,
                                        Bivalent = NA, Active = NA, Divergent = NA, Heterochromatin = NA, Quies = NA)
## Reorder prev dataframe to enter the loop
AtZm.df.F <- AtZm.df.F[,c(1:25, 42:46, 26:41, 47:51)]
## Loop to fill the distance/simmilarity metrics between CS (states and groups) splitted by LECIF percentile bins
for(i in 1:21){# 16 CS and 5 CS functional groups
  cat(paste0(i / 21 * 100, ' % CS groups started \n'))
  for(j in 1:10){ # lecif score percentile bins
    cat(paste0(j / 10 * 100, ' % LECIF bins started \n'))
    CS.Simmilarity.LECIF.AtZm[j,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 1)
    CS.Simmilarity.LECIF.AtZm[j+10,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 2) 
    CS.Simmilarity.LECIF.AtZm[j+20,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 3)  
    CS.Simmilarity.LECIF.AtZm[j+30,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 4) 
    CS.Simmilarity.LECIF.AtZm[j+40,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 5) 
    CS.Simmilarity.LECIF.AtZm[j+50,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 6) 
    CS.Simmilarity.LECIF.AtZm[j+60,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 7) 
    CS.Simmilarity.LECIF.AtZm[j+70,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 8) 
    CS.Simmilarity.LECIF.AtZm[j+80,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 9) 
    CS.Simmilarity.LECIF.AtZm[j+90,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 10) 
    CS.Simmilarity.LECIF.AtZm[j+100,2+i] <- philentropy::euclidean(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F) 
    CS.Simmilarity.LECIF.AtZm[j+110,2+i] <- philentropy::sorensen(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+120,2+i] <- philentropy::canberra(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+130,2+i] <- philentropy::wave_hedges(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+140,2+i] <- philentropy::tanimoto(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+150,2+i] <- philentropy::fidelity(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+160,2+i] <- philentropy::squared_chord(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+170,2+i] <- philentropy::squared_chi_sq(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+180,2+i] <- philentropy::additive_symm_chi_sq(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+190,2+i] <- philentropy::topsoe(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F,unit = "log2")
    CS.Simmilarity.LECIF.AtZm[j+200,2+i] <- philentropy::kumar_johnson(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.AtZm[j+210,2+i] <- philentropy::manhattan(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+220,2+i] <- philentropy::lorentzian(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.AtZm[j+230,2+i] <- philentropy::czekanowski(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+240,2+i] <- philentropy::kumar_hassebrook(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+250,2+i] <- philentropy::squared_euclidean(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+260,2+i] <- philentropy::prob_symm_chi_sq(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+270,2+i] <- philentropy::kullback_leibler_distance(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.AtZm[j+280,2+i] <- philentropy::jensen_shannon(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.AtZm[j+290,2+i] <- philentropy::minkowski(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, n = 3)
    CS.Simmilarity.LECIF.AtZm[j+300,2+i] <- philentropy::soergel(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+310,2+i] <- philentropy::inner_product(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+320,2+i] <- philentropy::hellinger(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+330,2+i] <- philentropy::divergence_sq(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+340,2+i] <- philentropy::jeffreys(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.AtZm[j+350,2+i] <- philentropy::jensen_difference(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.AtZm[j+360,2+i] <- philentropy::chebyshev(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+370,2+i] <- philentropy::kulczynski_d(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.AtZm[j+380,2+i] <- philentropy::harmonic_mean_dist(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+390,2+i] <- philentropy::matusita(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+400,2+i] <- philentropy::neyman_chi_sq(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.AtZm[j+410,2+i] <- philentropy::clark_sq(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+420,2+i] <- philentropy::k_divergence(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.AtZm[j+430,2+i] <- philentropy::taneja(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.AtZm[j+440,2+i] <- proxyC::dist(as.matrix(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = "maximum", margin = 2)[1,2]
    CS.Simmilarity.LECIF.AtZm[j+450,2+i] <- proxyC::dist(as.matrix(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = "hamming", margin = 2)[1,2]
    cat(paste0(j / 10 * 100, ' % LECIF bins completed \n'))
  }
  cat(paste0(i / 21 * 100, ' % CS groups completed \n'))
}

write.table(AtZm.df.F, "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/AtZm_CSsimmilarity.DF10.txt", col.names = T, row.names = F, sep = "\t", quote = F)
CS.Simmilarity.LECIF.AtZm.10 <- CS.Simmilarity.LECIF.AtZm
write.table(CS.Simmilarity.LECIF.AtZm, file = "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/AtZm_CSsimmilarity.DF.Distances10.txt", col.names = T, row.names = F, sep = "\t", quote = F)

## Example code to check if metrics are distance or simmilarity and max/min values for each
ade4::dist.binary(t(data.frame(CS1.At = AtZm.df.F3$CS2.At[which(AtZm.df.F$Percentile.Rank.Bin == "First")], CS1.Zm = AtZm.df.F3$CS2.At[which(AtZm.df.F$Percentile.Rank.Bin == "First")])), method = 5)
ade4::dist.binary(t(data.frame(CS1.At = AtZm.df.F$CS2.At[which(AtZm.df.F$Percentile.Rank.Bin == "Second")], CS1.Zm = AtZm.df.F$CS2.Zm[which(AtZm.df.F$Percentile.Rank.Bin == "Second")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.At = AtZm.df.F$CS2.At[which(AtZm.df.F$Percentile.Rank.Bin == "Thirdth")], CS1.Zm = AtZm.df.F$CS2.Zm[which(AtZm.df.F$Percentile.Rank.Bin == "Thirdth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.At = AtZm.df.F$CS2.At[which(AtZm.df.F$Percentile.Rank.Bin == "Fourth")], CS1.Zm = AtZm.df.F$CS2.Zm[which(AtZm.df.F$Percentile.Rank.Bin == "Fourth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.At = AtZm.df.F$CS2.At[which(AtZm.df.F$Percentile.Rank.Bin == "Fifth")], CS1.Zm = AtZm.df.F$CS2.Zm[which(AtZm.df.F$Percentile.Rank.Bin == "Fifth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.At = AtZm.df.F$CS2.At[which(AtZm.df.F$Percentile.Rank.Bin == "Sixth")], CS1.Zm = AtZm.df.F$CS2.Zm[which(AtZm.df.F$Percentile.Rank.Bin == "Sixth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.At = AtZm.df.F$CS2.At[which(AtZm.df.F$Percentile.Rank.Bin == "Seventh")], CS1.Zm = AtZm.df.F$CS2.Zm[which(AtZm.df.F$Percentile.Rank.Bin == "Seventh")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.At = AtZm.df.F$CS2.At[which(AtZm.df.F$Percentile.Rank.Bin == "Eigth")], CS1.Zm = AtZm.df.F$CS2.Zm[which(AtZm.df.F$Percentile.Rank.Bin == "Eigth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.At = AtZm.df.F$CS2.At[which(AtZm.df.F$Percentile.Rank.Bin == "Nineth")], CS1.Zm = AtZm.df.F$CS2.Zm[which(AtZm.df.F$Percentile.Rank.Bin == "Nineth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.At = AtZm.df.F$CS2.At[which(AtZm.df.F$Percentile.Rank.Bin == "Tenth")], CS1.Zm = AtZm.df.F$CS2.Zm[which(AtZm.df.F$Percentile.Rank.Bin == "Tenth")])), method = 2)
philentropy::jaccard(P =  AtZm.df.F3$CS1.At[which(AtZm.df.F$Percentile.Rank.Bin == "First")], Q = AtZm.df.F3$CS1.At[which(AtZm.df.F$Percentile.Rank.Bin == "First")], testNA = F)
proxyC::dist(as.matrix(data.frame(Gene1 = AtZm.df.F$CS1.At[which(AtZm.df.F$Percentile.Rank.Bin == "First")], Gene2 = AtZm.df.F$CS1.At[which(AtZm.df.F$Percentile.Rank.Bin == "First")])), method = "hamming", margin = 2)[1,2]

### with less bins (5 or 3)

AtZm.df.F$Percentile.Rank.Bin.Five <- NA
AtZm.df.F$Percentile.Rank.Bin.Five[which(AtZm.df.F$Percentile.Rank.Score < 0.2)] <- "First" # 165105 regions
AtZm.df.F$Percentile.Rank.Bin.Five[which(AtZm.df.F$Percentile.Rank.Score >= 0.2 & AtZm.df.F$Percentile.Rank.Score < 0.4)] <- "Second" # 165107
AtZm.df.F$Percentile.Rank.Bin.Five[which(AtZm.df.F$Percentile.Rank.Score >= 0.4 & AtZm.df.F$Percentile.Rank.Score < 0.6)] <- "Thirdth" # 165104
AtZm.df.F$Percentile.Rank.Bin.Five[which(AtZm.df.F$Percentile.Rank.Score >= 0.6 & AtZm.df.F$Percentile.Rank.Score < 0.8)] <- "Fourth" # 165064
AtZm.df.F$Percentile.Rank.Bin.Five[which(AtZm.df.F$Percentile.Rank.Score >= 0.8 & AtZm.df.F$Percentile.Rank.Score <= 1)] <- "Fifth" # 165095
(165105 + 165107 + 165104 + 165064 + 165095) == nrow(AtZm.df.F) # check


CS.Simmilarity.LECIF.AtZm <- data.frame(Bin = rep(c("First", "Second", "Thirdth", "Fourth", "Fifth"),46),
                                        Metric = c(rep(c("Jaccard"),5),rep(c("SMC"),5),rep(c("Sokal"),5), rep(c("Tanimoto"),5),
                                                   rep(c("Dice"),5), rep(c("Hanman"),5), rep(c("Ochiai"),5), rep(c("Sokal2"),5),
                                                   rep(c("Phi"),5), rep(c("S2"),5), 
                                                   rep(c("Eucledian"),5), rep(c("Sorensen"),5), rep(c("Canberra"),5), rep(c("Wavehedges"),5), 
                                                   rep(c("Tanimoto2"),5), rep(c("Fidelity"),5), rep(c("SquaredChord"),5), rep(c("SquaredChi"),5), 
                                                   rep(c("AdditiveSym"),5), rep(c("Topsoe"),5), rep(c("Kumarjhonson"),5), rep(c("Manhattan"),5),
                                                   rep(c("Lorentzian"),5), rep(c("Czekanowski"),5), rep(c("Hassebrook"),5),  rep(c("SquaredEucledian"),5),
                                                   rep(c("ProbSymm"),5), rep(c("Kullbackleiber"),5), rep(c("Jensenshannon"),5), rep(c("Minkowski"),5), 
                                                   rep(c("Soergel"),5), rep(c("Innerproduct"),5), rep(c("Hellinger"),5), rep(c("Divergence"),5), 
                                                   rep(c("Jeffreys"),5), rep(c("Jensendifference"),5), rep(c("Chebysev"),5), rep(c("Kulczynskid"),5), 
                                                   rep(c("Harmonicmean"),5), rep(c("Matusita"),5), rep(c("Neyman"),5), rep(c("Clark"),5), 
                                                   rep(c("Kdivergence"),5), rep(c("Taneja"),5), 
                                                   rep(c("Maximum"),5), rep(c("Hamming"),5)),
                                        CS1 = NA, CS2 = NA, CS3 = NA, CS4 = NA,
                                        CS5 = NA, CS6 = NA, CS7 = NA, CS8 = NA,
                                        CS9 = NA, CS10 = NA, CS11 = NA, CS12 = NA,
                                        CS13 = NA, CS14 = NA, CS15 = NA, CS16 = NA,
                                        Bivalent = NA, Active = NA, Divergent = NA, Heterochromatin = NA, Quies = NA)

for(i in 1:21){# 16 CS and 5 CS functional groups
  cat(paste0(i / 21 * 100, ' % CS groups started \n'))
  for(j in 1:5){ # lecif score percentile bins
    cat(paste0(j / 5 * 100, ' % LECIF bins started \n'))
    CS.Simmilarity.LECIF.AtZm[j,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 1)
    CS.Simmilarity.LECIF.AtZm[j+5,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 2) 
    CS.Simmilarity.LECIF.AtZm[j+10,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 3)  
    CS.Simmilarity.LECIF.AtZm[j+15,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 4) 
    CS.Simmilarity.LECIF.AtZm[j+20,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 5) 
    CS.Simmilarity.LECIF.AtZm[j+25,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 6) 
    CS.Simmilarity.LECIF.AtZm[j+30,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 7) 
    CS.Simmilarity.LECIF.AtZm[j+35,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 8) 
    CS.Simmilarity.LECIF.AtZm[j+40,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 9) 
    CS.Simmilarity.LECIF.AtZm[j+45,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 10) 
    CS.Simmilarity.LECIF.AtZm[j+50,2+i] <- philentropy::euclidean(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F) 
    CS.Simmilarity.LECIF.AtZm[j+55,2+i] <- philentropy::sorensen(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+60,2+i] <- philentropy::canberra(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+65,2+i] <- philentropy::wave_hedges(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+70,2+i] <- philentropy::tanimoto(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+75,2+i] <- philentropy::fidelity(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+80,2+i] <- philentropy::squared_chord(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+85,2+i] <- philentropy::squared_chi_sq(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+90,2+i] <- philentropy::additive_symm_chi_sq(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+95,2+i] <- philentropy::topsoe(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F,unit = "log2")
    CS.Simmilarity.LECIF.AtZm[j+100,2+i] <- philentropy::kumar_johnson(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.AtZm[j+105,2+i] <- philentropy::manhattan(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+110,2+i] <- philentropy::lorentzian(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.AtZm[j+115,2+i] <- philentropy::czekanowski(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+120,2+i] <- philentropy::kumar_hassebrook(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+125,2+i] <- philentropy::squared_euclidean(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+130,2+i] <- philentropy::prob_symm_chi_sq(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+135,2+i] <- philentropy::kullback_leibler_distance(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.AtZm[j+140,2+i] <- philentropy::jensen_shannon(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.AtZm[j+145,2+i] <- philentropy::minkowski(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, n = 3)
    CS.Simmilarity.LECIF.AtZm[j+150,2+i] <- philentropy::soergel(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+155,2+i] <- philentropy::inner_product(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+160,2+i] <- philentropy::hellinger(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+165,2+i] <- philentropy::divergence_sq(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+170,2+i] <- philentropy::jeffreys(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.AtZm[j+175,2+i] <- philentropy::jensen_difference(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.AtZm[j+180,2+i] <- philentropy::chebyshev(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+185,2+i] <- philentropy::kulczynski_d(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.AtZm[j+190,2+i] <- philentropy::harmonic_mean_dist(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+195,2+i] <- philentropy::matusita(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+200,2+i] <- philentropy::neyman_chi_sq(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.AtZm[j+205,2+i] <- philentropy::clark_sq(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+210,2+i] <- philentropy::k_divergence(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.AtZm[j+215,2+i] <- philentropy::taneja(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.AtZm[j+220,2+i] <- proxyC::dist(as.matrix(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = "maximum", margin = 2)[1,2]
    CS.Simmilarity.LECIF.AtZm[j+225,2+i] <- proxyC::dist(as.matrix(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = "hamming", margin = 2)[1,2]
    cat(paste0(j / 5 * 100, ' % LECIF bins completed \n'))
  }
  cat(paste0(i / 21 * 100, ' % CS groups completed \n'))
}

write.table(AtZm.df.F, "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/AtZm_CSsimmilarity.DF5.txt", col.names = T, row.names = F, sep = "\t", quote = F)
CS.Simmilarity.LECIF.AtZm.5 <- CS.Simmilarity.LECIF.AtZm
write.table(CS.Simmilarity.LECIF.AtZm, file = "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/AtZm_CSsimmilarity.DF.Distances5.txt", col.names = T, row.names = F, sep = "\t", quote = F)
View(CS.Simmilarity.LECIF.AtZm.5)

AtZm.df.F$Percentile.Rank.Bin.Three <- NA
AtZm.df.F$Percentile.Rank.Bin.Three[which(AtZm.df.F$Percentile.Rank.Score < 0.333333)] <- "First" # 275161 regions
AtZm.df.F$Percentile.Rank.Bin.Three[which(AtZm.df.F$Percentile.Rank.Score >= 0.333333 & AtZm.df.F$Percentile.Rank.Score < 0.666666)] <- "Second" # 275161
AtZm.df.F$Percentile.Rank.Bin.Three[which(AtZm.df.F$Percentile.Rank.Score >= 0.666666 & AtZm.df.F$Percentile.Rank.Score <= 1)] <- "Thirdth" # 275153
( 275161 + 275161 + 275153) == nrow(AtZm.df.F) # check


CS.Simmilarity.LECIF.AtZm <- data.frame(Bin = rep(c("First", "Second", "Thirdth"),46),
                                        Metric = c(rep(c("Jaccard"),3),rep(c("SMC"),3),rep(c("Sokal"),3), rep(c("Tanimoto"),3),
                                                   rep(c("Dice"),3), rep(c("Hanman"),3), rep(c("Ochiai"),3), rep(c("Sokal2"),3),
                                                   rep(c("Phi"),3), rep(c("S2"),3), 
                                                   rep(c("Eucledian"),3), rep(c("Sorensen"),3), rep(c("Canberra"),3), rep(c("Wavehedges"),3), 
                                                   rep(c("Tanimoto2"),3), rep(c("Fidelity"),3), rep(c("SquaredChord"),3), rep(c("SquaredChi"),3), 
                                                   rep(c("AdditiveSym"),3), rep(c("Topsoe"),3), rep(c("Kumarjhonson"),3), rep(c("Manhattan"),3),
                                                   rep(c("Lorentzian"),3), rep(c("Czekanowski"),3), rep(c("Hassebrook"),3),  rep(c("SquaredEucledian"),3),
                                                   rep(c("ProbSymm"),3), rep(c("Kullbackleiber"),3), rep(c("Jensenshannon"),3), rep(c("Minkowski"),3), 
                                                   rep(c("Soergel"),3), rep(c("Innerproduct"),3), rep(c("Hellinger"),3), rep(c("Divergence"),3), 
                                                   rep(c("Jeffreys"),3), rep(c("Jensendifference"),3), rep(c("Chebysev"),3), rep(c("Kulczynskid"),3), 
                                                   rep(c("Harmonicmean"),3), rep(c("Matusita"),3), rep(c("Neyman"),3), rep(c("Clark"),3), 
                                                   rep(c("Kdivergence"),3), rep(c("Taneja"),3), 
                                                   rep(c("Maximum"),3), rep(c("Hamming"),3)),
                                        CS1 = NA, CS2 = NA, CS3 = NA, CS4 = NA,
                                        CS5 = NA, CS6 = NA, CS7 = NA, CS8 = NA,
                                        CS9 = NA, CS10 = NA, CS11 = NA, CS12 = NA,
                                        CS13 = NA, CS14 = NA, CS15 = NA, CS16 = NA,
                                        Bivalent = NA, Active = NA, Divergent = NA, Heterochromatin = NA, Quies = NA)

for(i in 1:21){# 16 CS and 5 CS functional groups
  cat(paste0(i / 21 * 100, ' % CS groups started \n'))
  for(j in 1:3){ # lecif score percentile bins
    cat(paste0(j / 3 * 100, ' % LECIF bins started \n'))
    CS.Simmilarity.LECIF.AtZm[j,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 1)
    CS.Simmilarity.LECIF.AtZm[j+3,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 2) 
    CS.Simmilarity.LECIF.AtZm[j+6,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 3)  
    CS.Simmilarity.LECIF.AtZm[j+9,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 4) 
    CS.Simmilarity.LECIF.AtZm[j+12,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 5) 
    CS.Simmilarity.LECIF.AtZm[j+15,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 6) 
    CS.Simmilarity.LECIF.AtZm[j+18,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 7) 
    CS.Simmilarity.LECIF.AtZm[j+21,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 8) 
    CS.Simmilarity.LECIF.AtZm[j+24,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 9) 
    CS.Simmilarity.LECIF.AtZm[j+27,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = 10) 
    CS.Simmilarity.LECIF.AtZm[j+30,2+i] <- philentropy::euclidean(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F) 
    CS.Simmilarity.LECIF.AtZm[j+33,2+i] <- philentropy::sorensen(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+36,2+i] <- philentropy::canberra(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+39,2+i] <- philentropy::wave_hedges(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+42,2+i] <- philentropy::tanimoto(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+45,2+i] <- philentropy::fidelity(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+48,2+i] <- philentropy::squared_chord(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+51,2+i] <- philentropy::squared_chi_sq(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+54,2+i] <- philentropy::additive_symm_chi_sq(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+57,2+i] <- philentropy::topsoe(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F,unit = "log2")
    CS.Simmilarity.LECIF.AtZm[j+60,2+i] <- philentropy::kumar_johnson(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.AtZm[j+63,2+i] <- philentropy::manhattan(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+66,2+i] <- philentropy::lorentzian(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.AtZm[j+69,2+i] <- philentropy::czekanowski(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+72,2+i] <- philentropy::kumar_hassebrook(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+75,2+i] <- philentropy::squared_euclidean(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+78,2+i] <- philentropy::prob_symm_chi_sq(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+81,2+i] <- philentropy::kullback_leibler_distance(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.AtZm[j+84,2+i] <- philentropy::jensen_shannon(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.AtZm[j+87,2+i] <- philentropy::minkowski(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, n = 3)
    CS.Simmilarity.LECIF.AtZm[j+90,2+i] <- philentropy::soergel(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+93,2+i] <- philentropy::inner_product(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+96,2+i] <- philentropy::hellinger(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+99,2+i] <- philentropy::divergence_sq(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+102,2+i] <- philentropy::jeffreys(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.AtZm[j+105,2+i] <- philentropy::jensen_difference(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.AtZm[j+108,2+i] <- philentropy::chebyshev(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+111,2+i] <- philentropy::kulczynski_d(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.AtZm[j+114,2+i] <- philentropy::harmonic_mean_dist(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+117,2+i] <- philentropy::matusita(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+120,2+i] <- philentropy::neyman_chi_sq(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.AtZm[j+123,2+i] <- philentropy::clark_sq(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.AtZm[j+126,2+i] <- philentropy::k_divergence(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.AtZm[j+129,2+i] <- philentropy::taneja(P = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], Q = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.AtZm[j+132,2+i] <- proxyC::dist(as.matrix(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = "maximum", margin = 2)[1,2]
    CS.Simmilarity.LECIF.AtZm[j+135,2+i] <- proxyC::dist(as.matrix(data.frame(CS.At = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),9+i], CS.Zm = AtZm.df.F[which(AtZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.AtZm[j,1]),30+i])), method = "hamming", margin = 2)[1,2]
    cat(paste0(j / 3 * 100, ' % LECIF bins completed \n'))
  }
  cat(paste0(i / 21 * 100, ' % CS groups completed \n'))
}

write.table(AtZm.df.F, "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/AtZm_CSsimmilarity.DF3.txt", col.names = T, row.names = F, sep = "\t", quote = F)
CS.Simmilarity.LECIF.AtZm.3 <- CS.Simmilarity.LECIF.AtZm
write.table(CS.Simmilarity.LECIF.AtZm, file = "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/AtZm_CSsimmilarity.DF.Distances3.txt", col.names = T, row.names = F, sep = "\t", quote = F)
View(CS.Simmilarity.LECIF.AtZm.3) ## Dice-like would be the best option

###############################################
########## 1.2.3 Oryza - Arabidopsis ######### 
##############################################

AtOs_regions <- read.delim("H:/Taiotactical/Taiotactical_LECIF/example/check_os-at/all.m", header = F, sep = "\t") 
OsAt_regions <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/Os-At_LECIFscore.bedgraph", header = F, sep = "\t") 
OsAt_bed <- OsAt_regions[,c(1,2,3)]
AtOs_bed <- AtOs_regions[,c(1,2,3)]
AtOs_bed <- AtOs_bed[grep("Chr", AtOs_bed$V1),]
AtOs_bed$V3 <- as.numeric(AtOs_bed$V3) - 1
AtOs_bed$V2 <- as.numeric(AtOs_bed$V2)
nrow(OsAt_bed) == nrow(AtOs_bed) # check the same number of at aligning regions to os
AtOs_bed$V3 <- AtOs_bed$V3 + (OsAt_bed$V3 - OsAt_bed$V2)
length(which(((OsAt_bed$V3 - OsAt_bed$V2) == (AtOs_bed$V3) -  (AtOs_bed$V2)) == TRUE)) # check same expansion of the first aligning base
write.table(AtOs_bed,"H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/At_OsAligningRegions.bed", quote = F, sep = "\t", col.names = F, row.names = F)
# sort -k1,1 -k2,2n Os_AtAligningRegions.bed > Os_AtAligningRegions.sorted.bed
# bedtools intersect -a Os_AtAligningRegions.sorted.bed -b hihmm.model1.K15.Os.ReMapped.ReNamed.Reduced.Reordered.bed -wa -wb | cut -f1,2,3,7 > Os_AtAligningRegions.hiHMM.sorted.bed
OsAt_hiHMMInfo <- LECIF_OsAt[,c(1,2,3,5)] 
AtOs_hiHMMInfo <- read.table("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/At_OsAligningRegions.hiHMM.sorted.bed", header = F, sep = "\t")
OsAt.df.F <- data.frame(Chr.Os = OsAt_bed$V1, Start.Os = OsAt_bed$V2, End.Os = OsAt_bed$V3,
                        Chr.At = AtOs_bed$V1, Start.At = AtOs_bed$V2, End.At = AtOs_bed$V3,
                        LECIF.Score = OsAt_regions$V4, Percentile.Rank.Score = percent_rank(x = OsAt_regions$V4), Percentile.Rank.Bin = NA,
                        CS1.Os = 0, CS2.Os = 0, CS3.Os = 0, CS4.Os = 0, CS5.Os = 0, CS6.Os = 0,
                        CS7.Os = 0, CS8.Os = 0, CS9.Os = 0, CS10.Os = 0, CS11.Os = 0, CS12.Os = 0,
                        CS13.Os = 0, CS14.Os = 0, CS15.Os = 0, CS16.Os = 0,
                        CS1.At = 0, CS2.At = 0, CS3.At = 0, CS4.At = 0, CS5.At = 0, CS6.At = 0,
                        CS7.At = 0, CS8.At = 0, CS9.At = 0, CS10.At = 0, CS11.At = 0, CS12.At = 0,
                        CS13.At = 0, CS14.At = 0, CS15.At = 0, CS16.At = 0,
                        Bivalent.Os = 0, Active.Os = 0, Divergent.Os = 0, Heterochromatin.Os = 0, Quies.Os = 0,
                        Bivalent.At = 0, Active.At = 0, Divergent.At = 0, Heterochromatin.At = 0, Quies.At = 0)
length(percent_rank(x = OsAt_regions$V4)) == nrow(OsAt.df.F) # check


OsAt.df.F$Percentile.Rank.Bin[which(OsAt.df.F$Percentile.Rank.Score < 0.1)] <- "First" # 120212 regions
OsAt.df.F$Percentile.Rank.Bin[which(OsAt.df.F$Percentile.Rank.Score >= 0.1 & OsAt.df.F$Percentile.Rank.Score < 0.2)] <- "Second" # 120140
OsAt.df.F$Percentile.Rank.Bin[which(OsAt.df.F$Percentile.Rank.Score >= 0.2 & OsAt.df.F$Percentile.Rank.Score < 0.3)] <- "Thirdth" # 120190
OsAt.df.F$Percentile.Rank.Bin[which(OsAt.df.F$Percentile.Rank.Score >= 0.3 & OsAt.df.F$Percentile.Rank.Score < 0.4)] <- "Fourth" # 121832
OsAt.df.F$Percentile.Rank.Bin[which(OsAt.df.F$Percentile.Rank.Score >= 0.4 & OsAt.df.F$Percentile.Rank.Score < 0.5)] <- "Fifth" # 119512
OsAt.df.F$Percentile.Rank.Bin[which(OsAt.df.F$Percentile.Rank.Score >= 0.5 & OsAt.df.F$Percentile.Rank.Score < 0.6)] <- "Sixth" # 119528
OsAt.df.F$Percentile.Rank.Bin[which(OsAt.df.F$Percentile.Rank.Score >= 0.6 & OsAt.df.F$Percentile.Rank.Score < 0.7)] <- "Seventh" # 120058
OsAt.df.F$Percentile.Rank.Bin[which(OsAt.df.F$Percentile.Rank.Score >= 0.7 & OsAt.df.F$Percentile.Rank.Score < 0.8)] <- "Eigth" # 120116
OsAt.df.F$Percentile.Rank.Bin[which(OsAt.df.F$Percentile.Rank.Score >= 0.8 & OsAt.df.F$Percentile.Rank.Score < 0.9)] <- "Nineth" # 120004
OsAt.df.F$Percentile.Rank.Bin[which(OsAt.df.F$Percentile.Rank.Score >= 0.9 & OsAt.df.F$Percentile.Rank.Score <= 1)] <- "Tenth" # 120167
(120212 + 120140 + 120190 + 121832 + 119512 + 119528 + 120058 + 120116 + 120004 + 120167) == nrow(OsAt.df.F) # check

for(i in 1:nrow(OsAt.df.F)){
  cat(paste0(i / nrow(OsAt.df.F) * 100, ' % completed \n'))
  OsCSinfo <- OsAt_hiHMMInfo$V5[which(OsAt_hiHMMInfo$V1 == OsAt.df.F$Chr.Os[i] & OsAt_hiHMMInfo$V2 == OsAt.df.F$Start.Os[i] & OsAt_hiHMMInfo$V3 == OsAt.df.F$End.Os[i])]
  AtCSinfo <- AtOs_hiHMMInfo$V4[which(AtOs_hiHMMInfo$V1 == OsAt.df.F$Chr.At[i] & AtOs_hiHMMInfo$V2 == OsAt.df.F$Start.At[i] & AtOs_hiHMMInfo$V3 == OsAt.df.F$End.At[i])]
  OsAt.df.F[i,which(c(paste0("CS",1:16)) %in% OsCSinfo)+9] <- 1
  OsAt.df.F[i,which(c(paste0("CS",1:16)) %in% AtCSinfo)+25] <- 1
  if(length(which(c("CS1", "CS2", "CS3", "CS4", "CS5") %in% AtCSinfo))>0){
    OsAt.df.F$Bivalent.At[i] <- 1
  }
  if(length(which(c("CS1", "CS2", "CS3", "CS4", "CS5") %in% OsCSinfo))>0){
    OsAt.df.F$Bivalent.Os[i] <- 1
  }
  if(length(which(c("CS6", "CS7", "CS8", "CS9") %in% AtCSinfo))>0){
    OsAt.df.F$Active.At[i] <- 1
  }
  if(length(which(c("CS6", "CS7", "CS8", "CS9") %in% OsCSinfo))>0){
    OsAt.df.F$Active.Os[i] <- 1
  }
  if(length(which(c("CS10") %in% AtCSinfo))>0){
    OsAt.df.F$Divergent.At[i] <- 1
  }
  if(length(which(c("CS10") %in% OsCSinfo))>0){
    OsAt.df.F$Divergent.Os[i] <- 1
  }
  if(length(which(c("CS11", "CS12", "CS13") %in% AtCSinfo))>0){
    OsAt.df.F$Heterochromatin.At[i] <- 1
  }
  if(length(which(c("CS11", "CS12", "CS13") %in% OsCSinfo))>0){
    OsAt.df.F$Heterochromatin.Os[i] <- 1
  }
  if(length(which(c("CS14", "CS15", "CS16") %in% AtCSinfo))>0){
    OsAt.df.F$Quies.At[i] <- 1
  }
  if(length(which(c("CS14", "CS15", "CS16") %in% OsCSinfo))>0){
    OsAt.df.F$Quies.Os[i] <- 1
  }
}

#write.table(AtOs.df.F, "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/AtOs_CSsimmilarity.DF.txt", col.names = T, row.names = F, sep = "\t", quote = F)
#AtOs.df.F. <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/AtOs_CSsimmilarity.DF.txt",header = T, row.names = F, sep = "\t", quote = F)
# loop to compute simmilarity between this huge binary vectos || metric performing well in boolean data but not too influences by super-sparsed data (low coverage from some states)
CS.Simmilarity.LECIF.OsAt <- data.frame(Bin = rep(c("First", "Second", "Thirdth", "Fourth", "Fifth", "Sixth", "Seventh", "Eigth", "Nineth", "Tenth"),46),
                                        Metric = c(rep(c("Jaccard"),10),rep(c("SMC"),10),rep(c("Sokal"),10), rep(c("Tanimoto"),10),
                                                   rep(c("Dice"),10), rep(c("Hanman"),10), rep(c("Ochiai"),10), rep(c("Sokal2"),10),
                                                   rep(c("Phi"),10), rep(c("S2"),10), 
                                                   rep(c("Eucledian"),10), rep(c("Sorensen"),10), rep(c("Canberra"),10), rep(c("Wavehedges"),10), 
                                                   rep(c("Tanimoto2"),10), rep(c("Fidelity"),10), rep(c("SquaredChord"),10), rep(c("SquaredChi"),10), 
                                                   rep(c("AdditiveSym"),10), rep(c("Topsoe"),10), rep(c("Kumarjhonson"),10), rep(c("Manhattan"),10),
                                                   rep(c("Lorentzian"),10), rep(c("Czekanowski"),10), rep(c("Hassebrook"),10),  rep(c("SquaredEucledian"),10),
                                                   rep(c("ProbSymm"),10), rep(c("Kullbackleiber"),10), rep(c("Jensenshannon"),10), rep(c("Minkowski"),10), 
                                                   rep(c("Soergel"),10), rep(c("Innerproduct"),10), rep(c("Hellinger"),10), rep(c("Divergence"),10), 
                                                   rep(c("Jeffreys"),10), rep(c("Jensendifference"),10), rep(c("Chebysev"),10), rep(c("Kulczynskid"),10), 
                                                   rep(c("Harmonicmean"),10), rep(c("Matusita"),10), rep(c("Neyman"),10), rep(c("Clark"),10), 
                                                   rep(c("Kdivergence"),10), rep(c("Taneja"),10), 
                                                   rep(c("Maximum"),10), rep(c("Hamming"),10)),
                                        CS1 = NA, CS2 = NA, CS3 = NA, CS4 = NA,
                                        CS5 = NA, CS6 = NA, CS7 = NA, CS8 = NA,
                                        CS9 = NA, CS10 = NA, CS11 = NA, CS12 = NA,
                                        CS13 = NA, CS14 = NA, CS15 = NA, CS16 = NA,
                                        Bivalent = NA, Active = NA, Divergent = NA, Heterochromatin = NA, Quies = NA)
## Reorder prev dataframe to enter the loop
OsAt.df.F <- OsAt.df.F[,c(1:25, 42:46, 26:41, 47:51)]
## Loop to fill the distance/simmilarity metrics between CS (states and groups) splited by LECIF percentile bins
for(i in 1:21){# 16 CS and 5 CS functional groups
  cat(paste0(i / 21 * 100, ' % CS groups started \n'))
  for(j in 1:10){ # lecif score percentile bins
    cat(paste0(j / 10 * 100, ' % LECIF bins started \n'))
    CS.Simmilarity.LECIF.OsAt[j,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 1)
    CS.Simmilarity.LECIF.OsAt[j+10,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 2) 
    CS.Simmilarity.LECIF.OsAt[j+20,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 3)  
    CS.Simmilarity.LECIF.OsAt[j+30,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 4) 
    CS.Simmilarity.LECIF.OsAt[j+40,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 5) 
    CS.Simmilarity.LECIF.OsAt[j+50,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 6) 
    CS.Simmilarity.LECIF.OsAt[j+60,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 7) 
    CS.Simmilarity.LECIF.OsAt[j+70,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 8) 
    CS.Simmilarity.LECIF.OsAt[j+80,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 9) 
    CS.Simmilarity.LECIF.OsAt[j+90,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 10) 
    CS.Simmilarity.LECIF.OsAt[j+100,2+i] <- philentropy::euclidean(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F) 
    CS.Simmilarity.LECIF.OsAt[j+110,2+i] <- philentropy::sorensen(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+120,2+i] <- philentropy::canberra(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+130,2+i] <- philentropy::wave_hedges(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+140,2+i] <- philentropy::tanimoto(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+150,2+i] <- philentropy::fidelity(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+160,2+i] <- philentropy::squared_chord(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+170,2+i] <- philentropy::squared_chi_sq(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+180,2+i] <- philentropy::additive_symm_chi_sq(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+190,2+i] <- philentropy::topsoe(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F,unit = "log2")
    CS.Simmilarity.LECIF.OsAt[j+200,2+i] <- philentropy::kumar_johnson(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.OsAt[j+210,2+i] <- philentropy::manhattan(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+220,2+i] <- philentropy::lorentzian(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.OsAt[j+230,2+i] <- philentropy::czekanowski(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+240,2+i] <- philentropy::kumar_hassebrook(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+250,2+i] <- philentropy::squared_euclidean(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+260,2+i] <- philentropy::prob_symm_chi_sq(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+270,2+i] <- philentropy::kullback_leibler_distance(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.OsAt[j+280,2+i] <- philentropy::jensen_shannon(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.OsAt[j+290,2+i] <- philentropy::minkowski(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, n = 3)
    CS.Simmilarity.LECIF.OsAt[j+300,2+i] <- philentropy::soergel(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+310,2+i] <- philentropy::inner_product(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+320,2+i] <- philentropy::hellinger(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+330,2+i] <- philentropy::divergence_sq(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+340,2+i] <- philentropy::jeffreys(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.OsAt[j+350,2+i] <- philentropy::jensen_difference(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.OsAt[j+360,2+i] <- philentropy::chebyshev(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+370,2+i] <- philentropy::kulczynski_d(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.OsAt[j+380,2+i] <- philentropy::harmonic_mean_dist(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+390,2+i] <- philentropy::matusita(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+400,2+i] <- philentropy::neyman_chi_sq(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.OsAt[j+410,2+i] <- philentropy::clark_sq(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+420,2+i] <- philentropy::k_divergence(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.OsAt[j+430,2+i] <- philentropy::taneja(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.OsAt[j+440,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = "maximum", margin = 2)[1,2]
    CS.Simmilarity.LECIF.OsAt[j+450,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = "hamming", margin = 2)[1,2]
    cat(paste0(j / 10 * 100, ' % LECIF bins completed \n'))
  }
  cat(paste0(i / 21 * 100, ' % CS groups completed \n'))
}

write.table(OsAt.df.F, "H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/OsAt_CSsimmilarity.DF10.txt", col.names = T, row.names = F, sep = "\t", quote = F)
CS.Simmilarity.LECIF.OsAt.10 <- CS.Simmilarity.LECIF.OsAt
write.table(CS.Simmilarity.LECIF.OsAt, file = "H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/OsAt_CSsimmilarity.DF.Distances10.txt", col.names = T, row.names = F, sep = "\t", quote = F)

## Example code to check if metrics are distance or simmilarity and max/min values for each
ade4::dist.binary(t(data.frame(CS1.Os = OsAt.df.F3$CS2.At[which(OsAt.df.F$Percentile.Rank.Bin == "First")], CS1.At = OsAt.df.F3$CS2.At[which(OsAt.df.F$Percentile.Rank.Bin == "First")])), method = 5)
ade4::dist.binary(t(data.frame(CS1.Os = OsAt.df.F$CS2.At[which(OsAt.df.F$Percentile.Rank.Bin == "Second")], CS1.At = OsAt.df.F$CS2.Os[which(OsAt.df.F$Percentile.Rank.Bin == "Second")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Os = OsAt.df.F$CS2.At[which(OsAt.df.F$Percentile.Rank.Bin == "Thirdth")], CS1.At = OsAt.df.F$CS2.Os[which(OsAt.df.F$Percentile.Rank.Bin == "Thirdth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Os = OsAt.df.F$CS2.At[which(OsAt.df.F$Percentile.Rank.Bin == "Fourth")], CS1.At = OsAt.df.F$CS2.Os[which(OsAt.df.F$Percentile.Rank.Bin == "Fourth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Os = OsAt.df.F$CS2.At[which(OsAt.df.F$Percentile.Rank.Bin == "Fifth")], CS1.At = OsAt.df.F$CS2.Os[which(OsAt.df.F$Percentile.Rank.Bin == "Fifth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Os = OsAt.df.F$CS2.At[which(OsAt.df.F$Percentile.Rank.Bin == "Sixth")], CS1.At = OsAt.df.F$CS2.Os[which(OsAt.df.F$Percentile.Rank.Bin == "Sixth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Os = OsAt.df.F$CS2.At[which(OsAt.df.F$Percentile.Rank.Bin == "Seventh")], CS1.At = OsAt.df.F$CS2.Os[which(OsAt.df.F$Percentile.Rank.Bin == "Seventh")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Os = OsAt.df.F$CS2.At[which(OsAt.df.F$Percentile.Rank.Bin == "Eigth")], CS1.At = OsAt.df.F$CS2.Os[which(OsAt.df.F$Percentile.Rank.Bin == "Eigth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Os = OsAt.df.F$CS2.At[which(OsAt.df.F$Percentile.Rank.Bin == "Nineth")], CS1.At = OsAt.df.F$CS2.Os[which(OsAt.df.F$Percentile.Rank.Bin == "Nineth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Os = OsAt.df.F$CS2.At[which(OsAt.df.F$Percentile.Rank.Bin == "Tenth")], CS1.At = OsAt.df.F$CS2.Os[which(OsAt.df.F$Percentile.Rank.Bin == "Tenth")])), method = 2)
philentropy::jaccard(P =  OsAt.df.F3$CS1.At[which(OsAt.df.F$Percentile.Rank.Bin == "First")], Q = OsAt.df.F3$CS1.At[which(OsAt.df.F$Percentile.Rank.Bin == "First")], testNA = F)
proxyC::dist(as.matrix(data.frame(Gene1 = OsAt.df.F$CS1.At[which(OsAt.df.F$Percentile.Rank.Bin == "First")], Gene2 = OsAt.df.F$CS1.At[which(OsAt.df.F$Percentile.Rank.Bin == "First")])), method = "hamming", margin = 2)[1,2]

### with less bins (5 or 3)

OsAt.df.F$Percentile.Rank.Bin.Five <- NA
OsAt.df.F$Percentile.Rank.Bin.Five[which(OsAt.df.F$Percentile.Rank.Score < 0.2)] <- "First" # 240352 regions
OsAt.df.F$Percentile.Rank.Bin.Five[which(OsAt.df.F$Percentile.Rank.Score >= 0.2 & OsAt.df.F$Percentile.Rank.Score < 0.4)] <- "Second" # 242022
OsAt.df.F$Percentile.Rank.Bin.Five[which(OsAt.df.F$Percentile.Rank.Score >= 0.4 & OsAt.df.F$Percentile.Rank.Score < 0.6)] <- "Thirdth" # 239040
OsAt.df.F$Percentile.Rank.Bin.Five[which(OsAt.df.F$Percentile.Rank.Score >= 0.6 & OsAt.df.F$Percentile.Rank.Score < 0.8)] <- "Fourth" # 240174
OsAt.df.F$Percentile.Rank.Bin.Five[which(OsAt.df.F$Percentile.Rank.Score >= 0.8 & OsAt.df.F$Percentile.Rank.Score <= 1)] <- "Fifth" # 240171
(240352 + 242022 + 239040 + 240174 + 240171) == nrow(OsAt.df.F) # check


CS.Simmilarity.LECIF.OsAt <- data.frame(Bin = rep(c("First", "Second", "Thirdth", "Fourth", "Fifth"),46),
                                        Metric = c(rep(c("Jaccard"),5),rep(c("SMC"),5),rep(c("Sokal"),5), rep(c("Tanimoto"),5),
                                                   rep(c("Dice"),5), rep(c("Hanman"),5), rep(c("Ochiai"),5), rep(c("Sokal2"),5),
                                                   rep(c("Phi"),5), rep(c("S2"),5), 
                                                   rep(c("Eucledian"),5), rep(c("Sorensen"),5), rep(c("Canberra"),5), rep(c("Wavehedges"),5), 
                                                   rep(c("Tanimoto2"),5), rep(c("Fidelity"),5), rep(c("SquaredChord"),5), rep(c("SquaredChi"),5), 
                                                   rep(c("AdditiveSym"),5), rep(c("Topsoe"),5), rep(c("Kumarjhonson"),5), rep(c("Manhattan"),5),
                                                   rep(c("Lorentzian"),5), rep(c("Czekanowski"),5), rep(c("Hassebrook"),5),  rep(c("SquaredEucledian"),5),
                                                   rep(c("ProbSymm"),5), rep(c("Kullbackleiber"),5), rep(c("Jensenshannon"),5), rep(c("Minkowski"),5), 
                                                   rep(c("Soergel"),5), rep(c("Innerproduct"),5), rep(c("Hellinger"),5), rep(c("Divergence"),5), 
                                                   rep(c("Jeffreys"),5), rep(c("Jensendifference"),5), rep(c("Chebysev"),5), rep(c("Kulczynskid"),5), 
                                                   rep(c("Harmonicmean"),5), rep(c("Matusita"),5), rep(c("Neyman"),5), rep(c("Clark"),5), 
                                                   rep(c("Kdivergence"),5), rep(c("Taneja"),5), 
                                                   rep(c("Maximum"),5), rep(c("Hamming"),5)),
                                        CS1 = NA, CS2 = NA, CS3 = NA, CS4 = NA,
                                        CS5 = NA, CS6 = NA, CS7 = NA, CS8 = NA,
                                        CS9 = NA, CS10 = NA, CS11 = NA, CS12 = NA,
                                        CS13 = NA, CS14 = NA, CS15 = NA, CS16 = NA,
                                        Bivalent = NA, Active = NA, Divergent = NA, Heterochromatin = NA, Quies = NA)

for(i in 1:21){# 16 CS and 5 CS functional groups
  cat(paste0(i / 21 * 100, ' % CS groups started \n'))
  for(j in 1:5){ # lecif score percentile bins
    cat(paste0(j / 5 * 100, ' % LECIF bins started \n'))
    CS.Simmilarity.LECIF.OsAt[j,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 1)
    CS.Simmilarity.LECIF.OsAt[j+5,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 2) 
    CS.Simmilarity.LECIF.OsAt[j+10,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 3)  
    CS.Simmilarity.LECIF.OsAt[j+15,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 4) 
    CS.Simmilarity.LECIF.OsAt[j+20,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 5) 
    CS.Simmilarity.LECIF.OsAt[j+25,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 6) 
    CS.Simmilarity.LECIF.OsAt[j+30,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 7) 
    CS.Simmilarity.LECIF.OsAt[j+35,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 8) 
    CS.Simmilarity.LECIF.OsAt[j+40,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 9) 
    CS.Simmilarity.LECIF.OsAt[j+45,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 10) 
    CS.Simmilarity.LECIF.OsAt[j+50,2+i] <- philentropy::euclidean(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F) 
    CS.Simmilarity.LECIF.OsAt[j+55,2+i] <- philentropy::sorensen(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+60,2+i] <- philentropy::canberra(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+65,2+i] <- philentropy::wave_hedges(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+70,2+i] <- philentropy::tanimoto(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+75,2+i] <- philentropy::fidelity(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+80,2+i] <- philentropy::squared_chord(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+85,2+i] <- philentropy::squared_chi_sq(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+90,2+i] <- philentropy::additive_symm_chi_sq(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+95,2+i] <- philentropy::topsoe(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F,unit = "log2")
    CS.Simmilarity.LECIF.OsAt[j+100,2+i] <- philentropy::kumar_johnson(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.OsAt[j+105,2+i] <- philentropy::manhattan(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+110,2+i] <- philentropy::lorentzian(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.OsAt[j+115,2+i] <- philentropy::czekanowski(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+120,2+i] <- philentropy::kumar_hassebrook(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+125,2+i] <- philentropy::squared_euclidean(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+130,2+i] <- philentropy::prob_symm_chi_sq(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+135,2+i] <- philentropy::kullback_leibler_distance(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.OsAt[j+140,2+i] <- philentropy::jensen_shannon(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.OsAt[j+145,2+i] <- philentropy::minkowski(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, n = 3)
    CS.Simmilarity.LECIF.OsAt[j+150,2+i] <- philentropy::soergel(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+155,2+i] <- philentropy::inner_product(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+160,2+i] <- philentropy::hellinger(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+165,2+i] <- philentropy::divergence_sq(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+170,2+i] <- philentropy::jeffreys(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.OsAt[j+175,2+i] <- philentropy::jensen_difference(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.OsAt[j+180,2+i] <- philentropy::chebyshev(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+185,2+i] <- philentropy::kulczynski_d(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.OsAt[j+190,2+i] <- philentropy::harmonic_mean_dist(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+195,2+i] <- philentropy::matusita(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+200,2+i] <- philentropy::neyman_chi_sq(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.OsAt[j+205,2+i] <- philentropy::clark_sq(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+210,2+i] <- philentropy::k_divergence(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.OsAt[j+215,2+i] <- philentropy::taneja(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.OsAt[j+220,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = "maximum", margin = 2)[1,2]
    CS.Simmilarity.LECIF.OsAt[j+225,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = "hamming", margin = 2)[1,2]
    cat(paste0(j / 5 * 100, ' % LECIF bins completed \n'))
  }
  cat(paste0(i / 21 * 100, ' % CS groups completed \n'))
}

write.table(OsAt.df.F, "H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/OsAt_CSsimmilarity.DF5.txt", col.names = T, row.names = F, sep = "\t", quote = F)
CS.Simmilarity.LECIF.OsAt.5 <- CS.Simmilarity.LECIF.OsAt
write.table(CS.Simmilarity.LECIF.OsAt, file = "H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/OsAt_CSsimmilarity.DF.Distances5.txt", col.names = T, row.names = F, sep = "\t", quote = F)
View(CS.Simmilarity.LECIF.OsAt.5)

OsAt.df.F$Percentile.Rank.Bin.Three <- NA
OsAt.df.F$Percentile.Rank.Bin.Three[which(OsAt.df.F$Percentile.Rank.Score < 0.333333)] <- "First" # 400590 regions
OsAt.df.F$Percentile.Rank.Bin.Three[which(OsAt.df.F$Percentile.Rank.Score >= 0.333333 & OsAt.df.F$Percentile.Rank.Score < 0.666666)] <- "Second" # 400606
OsAt.df.F$Percentile.Rank.Bin.Three[which(OsAt.df.F$Percentile.Rank.Score >= 0.666666 & OsAt.df.F$Percentile.Rank.Score <= 1)] <- "Thirdth" # 400563
( 400590 + 400606 + 400563) == nrow(OsAt.df.F) # check


CS.Simmilarity.LECIF.OsAt <- data.frame(Bin = rep(c("First", "Second", "Thirdth"),46),
                                        Metric = c(rep(c("Jaccard"),3),rep(c("SMC"),3),rep(c("Sokal"),3), rep(c("Tanimoto"),3),
                                                   rep(c("Dice"),3), rep(c("Hanman"),3), rep(c("Ochiai"),3), rep(c("Sokal2"),3),
                                                   rep(c("Phi"),3), rep(c("S2"),3), 
                                                   rep(c("Eucledian"),3), rep(c("Sorensen"),3), rep(c("Canberra"),3), rep(c("Wavehedges"),3), 
                                                   rep(c("Tanimoto2"),3), rep(c("Fidelity"),3), rep(c("SquaredChord"),3), rep(c("SquaredChi"),3), 
                                                   rep(c("AdditiveSym"),3), rep(c("Topsoe"),3), rep(c("Kumarjhonson"),3), rep(c("Manhattan"),3),
                                                   rep(c("Lorentzian"),3), rep(c("Czekanowski"),3), rep(c("Hassebrook"),3),  rep(c("SquaredEucledian"),3),
                                                   rep(c("ProbSymm"),3), rep(c("Kullbackleiber"),3), rep(c("Jensenshannon"),3), rep(c("Minkowski"),3), 
                                                   rep(c("Soergel"),3), rep(c("Innerproduct"),3), rep(c("Hellinger"),3), rep(c("Divergence"),3), 
                                                   rep(c("Jeffreys"),3), rep(c("Jensendifference"),3), rep(c("Chebysev"),3), rep(c("Kulczynskid"),3), 
                                                   rep(c("Harmonicmean"),3), rep(c("Matusita"),3), rep(c("Neyman"),3), rep(c("Clark"),3), 
                                                   rep(c("Kdivergence"),3), rep(c("Taneja"),3), 
                                                   rep(c("Maximum"),3), rep(c("Hamming"),3)),
                                        CS1 = NA, CS2 = NA, CS3 = NA, CS4 = NA,
                                        CS5 = NA, CS6 = NA, CS7 = NA, CS8 = NA,
                                        CS9 = NA, CS10 = NA, CS11 = NA, CS12 = NA,
                                        CS13 = NA, CS14 = NA, CS15 = NA, CS16 = NA,
                                        Bivalent = NA, Active = NA, Divergent = NA, Heterochromatin = NA, Quies = NA)

for(i in 1:21){# 16 CS and 5 CS functional groups
  cat(paste0(i / 21 * 100, ' % CS groups started \n'))
  for(j in 1:3){ # lecif score percentile bins
    cat(paste0(j / 3 * 100, ' % LECIF bins started \n'))
    CS.Simmilarity.LECIF.OsAt[j,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 1)
    CS.Simmilarity.LECIF.OsAt[j+3,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 2) 
    CS.Simmilarity.LECIF.OsAt[j+6,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 3)  
    CS.Simmilarity.LECIF.OsAt[j+9,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 4) 
    CS.Simmilarity.LECIF.OsAt[j+12,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 5) 
    CS.Simmilarity.LECIF.OsAt[j+15,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 6) 
    CS.Simmilarity.LECIF.OsAt[j+18,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 7) 
    CS.Simmilarity.LECIF.OsAt[j+21,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 8) 
    CS.Simmilarity.LECIF.OsAt[j+24,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 9) 
    CS.Simmilarity.LECIF.OsAt[j+27,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = 10) 
    CS.Simmilarity.LECIF.OsAt[j+30,2+i] <- philentropy::euclidean(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F) 
    CS.Simmilarity.LECIF.OsAt[j+33,2+i] <- philentropy::sorensen(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+36,2+i] <- philentropy::canberra(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+39,2+i] <- philentropy::wave_hedges(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+42,2+i] <- philentropy::tanimoto(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+45,2+i] <- philentropy::fidelity(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+48,2+i] <- philentropy::squared_chord(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+51,2+i] <- philentropy::squared_chi_sq(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+54,2+i] <- philentropy::additive_symm_chi_sq(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+57,2+i] <- philentropy::topsoe(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F,unit = "log2")
    CS.Simmilarity.LECIF.OsAt[j+60,2+i] <- philentropy::kumar_johnson(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.OsAt[j+63,2+i] <- philentropy::manhattan(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+66,2+i] <- philentropy::lorentzian(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.OsAt[j+69,2+i] <- philentropy::czekanowski(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+72,2+i] <- philentropy::kumar_hassebrook(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+75,2+i] <- philentropy::squared_euclidean(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+78,2+i] <- philentropy::prob_symm_chi_sq(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+81,2+i] <- philentropy::kullback_leibler_distance(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.OsAt[j+84,2+i] <- philentropy::jensen_shannon(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.OsAt[j+87,2+i] <- philentropy::minkowski(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, n = 3)
    CS.Simmilarity.LECIF.OsAt[j+90,2+i] <- philentropy::soergel(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+93,2+i] <- philentropy::inner_product(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+96,2+i] <- philentropy::hellinger(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+99,2+i] <- philentropy::divergence_sq(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+102,2+i] <- philentropy::jeffreys(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.OsAt[j+105,2+i] <- philentropy::jensen_difference(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.OsAt[j+108,2+i] <- philentropy::chebyshev(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+111,2+i] <- philentropy::kulczynski_d(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.OsAt[j+114,2+i] <- philentropy::harmonic_mean_dist(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+117,2+i] <- philentropy::matusita(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+120,2+i] <- philentropy::neyman_chi_sq(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.OsAt[j+123,2+i] <- philentropy::clark_sq(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsAt[j+126,2+i] <- philentropy::k_divergence(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.OsAt[j+129,2+i] <- philentropy::taneja(P = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], Q = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.OsAt[j+132,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = "maximum", margin = 2)[1,2]
    CS.Simmilarity.LECIF.OsAt[j+135,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Os = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),9+i], CS.At = OsAt.df.F[which(OsAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsAt[j,1]),30+i])), method = "hamming", margin = 2)[1,2]
    cat(paste0(j / 3 * 100, ' % LECIF bins completed \n'))
  }
  cat(paste0(i / 21 * 100, ' % CS groups completed \n'))
}

write.table(OsAt.df.F, "H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/OsAt_CSsimmilarity.DF3.txt", col.names = T, row.names = F, sep = "\t", quote = F)
CS.Simmilarity.LECIF.OsAt.3 <- CS.Simmilarity.LECIF.OsAt
write.table(CS.Simmilarity.LECIF.OsAt, file = "H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/OsAt_CSsimmilarity.DF.Distances3.txt", col.names = T, row.names = F, sep = "\t", quote = F)
View(CS.Simmilarity.LECIF.OsAt.3) ## Dice-like would be the best option

###############################################
############# 1.2.4 Oryza - Zea ############## 
##############################################

ZmOs_regions <- read.delim("H:/Taiotactical/Taiotactical_LECIF/example/check_os-zm/all.m", header = F, sep = "\t") 
OsZm_regions <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/Os-Zm_LECIFscore.bedgraph", header = F, sep = "\t")
OsZm_bed <- OsZm_regions[,c(1,2,3)]
ZmOs_bed <- ZmOs_regions[,c(1,2,3)]
ZmOs_bed <- ZmOs_bed[which(ZmOs_bed$V1 %in% c(1:10)),]
ZmOs_bed$V3 <- as.numeric(ZmOs_bed$V3) - 1
ZmOs_bed$V2 <- as.numeric(ZmOs_bed$V2)
nrow(OsZm_bed) == nrow(ZmOs_bed) # check the same number of Os aligning regions to Zm
ZmOs_bed$V3 <- ZmOs_bed$V3 + (OsZm_bed$V3 - OsZm_bed$V2)
length(which(((OsZm_bed$V3 - OsZm_bed$V2) == (ZmOs_bed$V3) -  (ZmOs_bed$V2)) == TRUE)) # check same expansion of the first aligning base
write.table(ZmOs_bed,"H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/Zm_OsAligningRegions.bed", quote = F, sep = "\t", col.names = F, row.names = F)
# samtools sort -k1,1 -k2,2n Zm_OsAligningRegions.bed > Zm_OsAligningRegions.sorted.bed
# bedtools intersect -a Zm_OsAligningRegions.sorted.bed -b hihmm.model1.K15.Zm.ReMapped.ReNamed.Reduced.Reordered.bed -wa -wb | cut -f1,2,3,7 > Os_OsAligningRegions.hiHMM.sorted.bed
OsZm_hiHMMInfo <- LECIF_OsZm[,c(1,2,3,5)] 
ZmOs_hiHMMInfo <- read.table("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/Zm_OsAligningRegions.hiHMM.sorted.bed", header = F, sep = "\t")
OsZm.df.F <- data.frame(Chr.Os = OsZm_bed$V1, Start.Os = OsZm_bed$V2, End.Os = OsZm_bed$V3,
                        Chr.Zm = ZmOs_bed$V1, Start.Zm = ZmOs_bed$V2, End.Zm = ZmOs_bed$V3,
                        LECIF.Score = OsZm_regions$V4, Percentile.Rank.Score = percent_rank(x = OsZm_regions$V4), Percentile.Rank.Bin = NA,
                        CS1.Os = 0, CS2.Os = 0, CS3.Os = 0, CS4.Os = 0, CS5.Os = 0, CS6.Os = 0,
                        CS7.Os = 0, CS8.Os = 0, CS9.Os = 0, CS10.Os = 0, CS11.Os = 0, CS12.Os = 0,
                        CS13.Os = 0, CS14.Os = 0, CS15.Os = 0, CS16.Os = 0,
                        CS1.Zm = 0, CS2.Zm = 0, CS3.Zm = 0, CS4.Zm = 0, CS5.Zm = 0, CS6.Zm = 0,
                        CS7.Zm = 0, CS8.Zm = 0, CS9.Zm = 0, CS10.Zm = 0, CS11.Zm = 0, CS12.Zm = 0,
                        CS13.Zm = 0, CS14.Zm = 0, CS15.Zm = 0, CS16.Zm = 0,
                        Bivalent.Os = 0, Active.Os = 0, Divergent.Os = 0, Heterochromatin.Os = 0, Quies.Os = 0,
                        Bivalent.Zm = 0, Active.Zm = 0, Divergent.Zm = 0, Heterochromatin.Zm = 0, Quies.Zm = 0)
length(percent_rank(x = OsZm_regions$V4)) == nrow(OsZm.df.F) # check

OsZm.df.F$Percentile.Rank.Bin[which(OsZm.df.F$Percentile.Rank.Score < 0.1)] <- "First" # 196364 regions
OsZm.df.F$Percentile.Rank.Bin[which(OsZm.df.F$Percentile.Rank.Score >= 0.1 & OsZm.df.F$Percentile.Rank.Score < 0.2)] <- "Second" # 194782
OsZm.df.F$Percentile.Rank.Bin[which(OsZm.df.F$Percentile.Rank.Score >= 0.2 & OsZm.df.F$Percentile.Rank.Score < 0.3)] <- "Thirdth" # 195581
OsZm.df.F$Percentile.Rank.Bin[which(OsZm.df.F$Percentile.Rank.Score >= 0.3 & OsZm.df.F$Percentile.Rank.Score < 0.4)] <- "Fourth" # 195836
OsZm.df.F$Percentile.Rank.Bin[which(OsZm.df.F$Percentile.Rank.Score >= 0.4 & OsZm.df.F$Percentile.Rank.Score < 0.5)] <- "Fifth" # 195452
OsZm.df.F$Percentile.Rank.Bin[which(OsZm.df.F$Percentile.Rank.Score >= 0.5 & OsZm.df.F$Percentile.Rank.Score < 0.6)] <- "Sixth" # 195356
OsZm.df.F$Percentile.Rank.Bin[which(OsZm.df.F$Percentile.Rank.Score >= 0.6 & OsZm.df.F$Percentile.Rank.Score < 0.7)] <- "Seventh" # 195694
OsZm.df.F$Percentile.Rank.Bin[which(OsZm.df.F$Percentile.Rank.Score >= 0.7 & OsZm.df.F$Percentile.Rank.Score < 0.8)] <- "Eigth" # 195670
OsZm.df.F$Percentile.Rank.Bin[which(OsZm.df.F$Percentile.Rank.Score >= 0.8 & OsZm.df.F$Percentile.Rank.Score < 0.9)] <- "Nineth" # 195322
OsZm.df.F$Percentile.Rank.Bin[which(OsZm.df.F$Percentile.Rank.Score >= 0.9 & OsZm.df.F$Percentile.Rank.Score <= 1)] <- "Tenth" # 195559
(196364 + 194782 + 195581 + 195836 + 195452 + 195356 + 195694 + 195670 + 195322 + 195559) == nrow(OsZm.df.F) # check

for(i in 1:nrow(OsZm.df.F)){
  cat(paste0(i / nrow(OsZm.df.F) * 100, ' % completed \n'))
  OsCSinfo <- OsZm_hiHMMInfo$V5[which(OsZm_hiHMMInfo$V1 == OsZm.df.F$Chr.Os[i] & OsZm_hiHMMInfo$V2 == OsZm.df.F$Start.Os[i] & OsZm_hiHMMInfo$V3 == OsZm.df.F$End.Os[i])]
  ZmCSinfo <- ZmOs_hiHMMInfo$V4[which(ZmOs_hiHMMInfo$V1 == OsZm.df.F$Chr.Zm[i] & ZmOs_hiHMMInfo$V2 == OsZm.df.F$Start.Zm[i] & ZmOs_hiHMMInfo$V3 == OsZm.df.F$End.Zm[i])]
  OsZm.df.F[i,which(c(paste0("CS",1:16)) %in% OsCSinfo)+9] <- 1
  OsZm.df.F[i,which(c(paste0("CS",1:16)) %in% ZmCSinfo)+25] <- 1
  if(length(which(c("CS1", "CS2", "CS3", "CS4", "CS5") %in% OsCSinfo))>0){
    OsZm.df.F$Bivalent.Os[i] <- 1
  }
  if(length(which(c("CS1", "CS2", "CS3", "CS4", "CS5") %in% ZmCSinfo))>0){
    OsZm.df.F$Bivalent.Zm[i] <- 1
  }
  if(length(which(c("CS6", "CS7", "CS8", "CS9") %in% OsCSinfo))>0){
    OsZm.df.F$Active.Os[i] <- 1
  }
  if(length(which(c("CS6", "CS7", "CS8", "CS9") %in% ZmCSinfo))>0){
    OsZm.df.F$Active.Zm[i] <- 1
  }
  if(length(which(c("CS10") %in% OsCSinfo))>0){
    OsZm.df.F$Divergent.Os[i] <- 1
  }
  if(length(which(c("CS10") %in% ZmCSinfo))>0){
    OsZm.df.F$Divergent.Zm[i] <- 1
  }
  if(length(which(c("CS11", "CS12", "CS13") %in% OsCSinfo))>0){
    OsZm.df.F$Heterochromatin.Os[i] <- 1
  }
  if(length(which(c("CS11", "CS12", "CS13") %in% ZmCSinfo))>0){
    OsZm.df.F$Heterochromatin.Zm[i] <- 1
  }
  if(length(which(c("CS14", "CS15", "CS16") %in% OsCSinfo))>0){
    OsZm.df.F$Quies.Os[i] <- 1
  }
  if(length(which(c("CS14", "CS15", "CS16") %in% ZmCSinfo))>0){
    OsZm.df.F$Quies.Zm[i] <- 1
  }
}

#write.table(AtZm.df.F, "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/OsZm_CSsimmilarity.DF.txt", col.names = T, row.names = F, sep = "\t", quote = F)
#OsZm.df.F. <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/OsZm_CSsimmilarity.DF.txt",header = T, row.names = F, sep = "\t", quote = F)
# loop to compute simmilarity between this huge binary vectos || metric performing well in boolean data but not too influences by super-sparsed data (low coverage from some states depending on the socre bin)
CS.Simmilarity.LECIF.OsZm <- data.frame(Bin = rep(c("First", "Second", "Thirdth", "Fourth", "Fifth", "Sixth", "Seventh", "Eigth", "Nineth", "Tenth"),46),
                                        Metric = c(rep(c("Jaccard"),10),rep(c("SMC"),10),rep(c("Sokal"),10), rep(c("Tanimoto"),10),
                                                   rep(c("Dice"),10), rep(c("Hanman"),10), rep(c("Ochiai"),10), rep(c("Sokal2"),10),
                                                   rep(c("Phi"),10), rep(c("S2"),10), 
                                                   rep(c("Eucledian"),10), rep(c("Sorensen"),10), rep(c("Canberra"),10), rep(c("Wavehedges"),10), 
                                                   rep(c("Tanimoto2"),10), rep(c("Fidelity"),10), rep(c("SquaredChord"),10), rep(c("SquaredChi"),10), 
                                                   rep(c("AdditiveSym"),10), rep(c("Topsoe"),10), rep(c("Kumarjhonson"),10), rep(c("Manhattan"),10),
                                                   rep(c("Lorentzian"),10), rep(c("Czekanowski"),10), rep(c("Hassebrook"),10),  rep(c("SquaredEucledian"),10),
                                                   rep(c("ProbSymm"),10), rep(c("Kullbackleiber"),10), rep(c("Jensenshannon"),10), rep(c("Minkowski"),10), 
                                                   rep(c("Soergel"),10), rep(c("Innerproduct"),10), rep(c("Hellinger"),10), rep(c("Divergence"),10), 
                                                   rep(c("Jeffreys"),10), rep(c("Jensendifference"),10), rep(c("Chebysev"),10), rep(c("Kulczynskid"),10), 
                                                   rep(c("Harmonicmean"),10), rep(c("Matusita"),10), rep(c("Neyman"),10), rep(c("Clark"),10), 
                                                   rep(c("Kdivergence"),10), rep(c("Taneja"),10), 
                                                   rep(c("Maximum"),10), rep(c("Hamming"),10)),
                                        CS1 = NA, CS2 = NA, CS3 = NA, CS4 = NA,
                                        CS5 = NA, CS6 = NA, CS7 = NA, CS8 = NA,
                                        CS9 = NA, CS10 = NA, CS11 = NA, CS12 = NA,
                                        CS13 = NA, CS14 = NA, CS15 = NA, CS16 = NA,
                                        Bivalent = NA, Active = NA, Divergent = NA, Heterochromatin = NA, Quies = NA)
## Reorder prev dataframe to enter the loop
OsZm.df.F <- OsZm.df.F[,c(1:25, 42:46, 26:41, 47:51)]
## Loop to fill the distance/simmilarity metrics between CS (states and groups) splitted by LECIF percentile bins
for(i in 1:21){# 16 CS and 5 CS functional groups
  cat(paste0(i / 21 * 100, ' % CS groups started \n'))
  for(j in 1:10){ # lecif score percentile bins
    cat(paste0(j / 10 * 100, ' % LECIF bins started \n'))
    CS.Simmilarity.LECIF.OsZm[j,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 1)
    CS.Simmilarity.LECIF.OsZm[j+10,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 2) 
    CS.Simmilarity.LECIF.OsZm[j+20,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 3)  
    CS.Simmilarity.LECIF.OsZm[j+30,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 4) 
    CS.Simmilarity.LECIF.OsZm[j+40,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 5) 
    CS.Simmilarity.LECIF.OsZm[j+50,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 6) 
    CS.Simmilarity.LECIF.OsZm[j+60,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 7) 
    CS.Simmilarity.LECIF.OsZm[j+70,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 8) 
    CS.Simmilarity.LECIF.OsZm[j+80,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 9) 
    CS.Simmilarity.LECIF.OsZm[j+90,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 10) 
    CS.Simmilarity.LECIF.OsZm[j+100,2+i] <- philentropy::euclidean(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F) 
    CS.Simmilarity.LECIF.OsZm[j+110,2+i] <- philentropy::sorensen(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+120,2+i] <- philentropy::canberra(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+130,2+i] <- philentropy::wave_hedges(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+140,2+i] <- philentropy::tanimoto(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+150,2+i] <- philentropy::fidelity(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+160,2+i] <- philentropy::squared_chord(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+170,2+i] <- philentropy::squared_chi_sq(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+180,2+i] <- philentropy::additive_symm_chi_sq(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+190,2+i] <- philentropy::topsoe(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F,unit = "log2")
    CS.Simmilarity.LECIF.OsZm[j+200,2+i] <- philentropy::kumar_johnson(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.OsZm[j+210,2+i] <- philentropy::manhattan(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+220,2+i] <- philentropy::lorentzian(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.OsZm[j+230,2+i] <- philentropy::czekanowski(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+240,2+i] <- philentropy::kumar_hassebrook(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+250,2+i] <- philentropy::squared_euclidean(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+260,2+i] <- philentropy::prob_symm_chi_sq(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+270,2+i] <- philentropy::kullback_leibler_distance(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.OsZm[j+280,2+i] <- philentropy::jensen_shannon(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.OsZm[j+290,2+i] <- philentropy::minkowski(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, n = 3)
    CS.Simmilarity.LECIF.OsZm[j+300,2+i] <- philentropy::soergel(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+310,2+i] <- philentropy::inner_product(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+320,2+i] <- philentropy::hellinger(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+330,2+i] <- philentropy::divergence_sq(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+340,2+i] <- philentropy::jeffreys(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.OsZm[j+350,2+i] <- philentropy::jensen_difference(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.OsZm[j+360,2+i] <- philentropy::chebyshev(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+370,2+i] <- philentropy::kulczynski_d(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.OsZm[j+380,2+i] <- philentropy::harmonic_mean_dist(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+390,2+i] <- philentropy::matusita(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+400,2+i] <- philentropy::neyman_chi_sq(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.OsZm[j+410,2+i] <- philentropy::clark_sq(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+420,2+i] <- philentropy::k_divergence(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.OsZm[j+430,2+i] <- philentropy::taneja(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.OsZm[j+440,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = "maximum", margin = 2)[1,2]
    CS.Simmilarity.LECIF.OsZm[j+450,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = "hamming", margin = 2)[1,2]
    cat(paste0(j / 10 * 100, ' % LECIF bins completed \n'))
  }
  cat(paste0(i / 21 * 100, ' % CS groups completed \n'))
}

write.table(OsZm.df.F, "H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/OsZm_CSsimmilarity.DF10.txt", col.names = T, row.names = F, sep = "\t", quote = F)
CS.Simmilarity.LECIF.OsZm.10 <- CS.Simmilarity.LECIF.OsZm
write.table(CS.Simmilarity.LECIF.OsZm, file = "H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/OsZm_CSsimmilarity.DF.Distances10.txt", col.names = T, row.names = F, sep = "\t", quote = F)

## Example code to check if metrics are distance or simmilarity and max/min values for each
ade4::dist.binary(t(data.frame(CS1.Os = OsZm.df.F3$CS2.At[which(OsZm.df.F$Percentile.Rank.Bin == "First")], CS1.Zm = OsZm.df.F3$CS2.At[which(OsZm.df.F$Percentile.Rank.Bin == "First")])), method = 5)
ade4::dist.binary(t(data.frame(CS1.Os = OsZm.df.F$CS2.At[which(OsZm.df.F$Percentile.Rank.Bin == "Second")], CS1.Zm = OsZm.df.F$CS2.Zm[which(OsZm.df.F$Percentile.Rank.Bin == "Second")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Os = OsZm.df.F$CS2.At[which(OsZm.df.F$Percentile.Rank.Bin == "Thirdth")], CS1.Zm = OsZm.df.F$CS2.Zm[which(OsZm.df.F$Percentile.Rank.Bin == "Thirdth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Os = OsZm.df.F$CS2.At[which(OsZm.df.F$Percentile.Rank.Bin == "Fourth")], CS1.Zm = OsZm.df.F$CS2.Zm[which(OsZm.df.F$Percentile.Rank.Bin == "Fourth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Os = OsZm.df.F$CS2.At[which(OsZm.df.F$Percentile.Rank.Bin == "Fifth")], CS1.Zm = OsZm.df.F$CS2.Zm[which(OsZm.df.F$Percentile.Rank.Bin == "Fifth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Os = OsZm.df.F$CS2.At[which(OsZm.df.F$Percentile.Rank.Bin == "Sixth")], CS1.Zm = OsZm.df.F$CS2.Zm[which(OsZm.df.F$Percentile.Rank.Bin == "Sixth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Os = OsZm.df.F$CS2.At[which(OsZm.df.F$Percentile.Rank.Bin == "Seventh")], CS1.Zm = OsZm.df.F$CS2.Zm[which(OsZm.df.F$Percentile.Rank.Bin == "Seventh")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Os = OsZm.df.F$CS2.At[which(OsZm.df.F$Percentile.Rank.Bin == "Eigth")], CS1.Zm = OsZm.df.F$CS2.Zm[which(OsZm.df.F$Percentile.Rank.Bin == "Eigth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Os = OsZm.df.F$CS2.At[which(OsZm.df.F$Percentile.Rank.Bin == "Nineth")], CS1.Zm = OsZm.df.F$CS2.Zm[which(OsZm.df.F$Percentile.Rank.Bin == "Nineth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Os = OsZm.df.F$CS2.At[which(OsZm.df.F$Percentile.Rank.Bin == "Tenth")], CS1.Zm = OsZm.df.F$CS2.Zm[which(OsZm.df.F$Percentile.Rank.Bin == "Tenth")])), method = 2)
philentropy::jaccard(P =  OsZm.df.F3$CS1.At[which(OsZm.df.F$Percentile.Rank.Bin == "First")], Q = OsZm.df.F3$CS1.At[which(OsZm.df.F$Percentile.Rank.Bin == "First")], testNA = F)
proxyC::dist(as.matrix(data.frame(Gene1 = OsZm.df.F$CS1.At[which(OsZm.df.F$Percentile.Rank.Bin == "First")], Gene2 = OsZm.df.F$CS1.At[which(OsZm.df.F$Percentile.Rank.Bin == "First")])), method = "hamming", margin = 2)[1,2]

### with less bins (5 or 3)

OsZm.df.F$Percentile.Rank.Bin.Five <- NA
OsZm.df.F$Percentile.Rank.Bin.Five[which(OsZm.df.F$Percentile.Rank.Score < 0.2)] <- "First" # 391146 regions
OsZm.df.F$Percentile.Rank.Bin.Five[which(OsZm.df.F$Percentile.Rank.Score >= 0.2 & OsZm.df.F$Percentile.Rank.Score < 0.4)] <- "Second" # 391417
OsZm.df.F$Percentile.Rank.Bin.Five[which(OsZm.df.F$Percentile.Rank.Score >= 0.4 & OsZm.df.F$Percentile.Rank.Score < 0.6)] <- "Thirdth" # 390808
OsZm.df.F$Percentile.Rank.Bin.Five[which(OsZm.df.F$Percentile.Rank.Score >= 0.6 & OsZm.df.F$Percentile.Rank.Score < 0.8)] <- "Fourth" # 391364
OsZm.df.F$Percentile.Rank.Bin.Five[which(OsZm.df.F$Percentile.Rank.Score >= 0.8 & OsZm.df.F$Percentile.Rank.Score <= 1)] <- "Fifth" # 390881
(391146 + 391417 + 390808 + 391364 + 390881) == nrow(OsZm.df.F) # check


CS.Simmilarity.LECIF.OsZm <- data.frame(Bin = rep(c("First", "Second", "Thirdth", "Fourth", "Fifth"),46),
                                        Metric = c(rep(c("Jaccard"),5),rep(c("SMC"),5),rep(c("Sokal"),5), rep(c("Tanimoto"),5),
                                                   rep(c("Dice"),5), rep(c("Hanman"),5), rep(c("Ochiai"),5), rep(c("Sokal2"),5),
                                                   rep(c("Phi"),5), rep(c("S2"),5), 
                                                   rep(c("Eucledian"),5), rep(c("Sorensen"),5), rep(c("Canberra"),5), rep(c("Wavehedges"),5), 
                                                   rep(c("Tanimoto2"),5), rep(c("Fidelity"),5), rep(c("SquaredChord"),5), rep(c("SquaredChi"),5), 
                                                   rep(c("AdditiveSym"),5), rep(c("Topsoe"),5), rep(c("Kumarjhonson"),5), rep(c("Manhattan"),5),
                                                   rep(c("Lorentzian"),5), rep(c("Czekanowski"),5), rep(c("Hassebrook"),5),  rep(c("SquaredEucledian"),5),
                                                   rep(c("ProbSymm"),5), rep(c("Kullbackleiber"),5), rep(c("Jensenshannon"),5), rep(c("Minkowski"),5), 
                                                   rep(c("Soergel"),5), rep(c("Innerproduct"),5), rep(c("Hellinger"),5), rep(c("Divergence"),5), 
                                                   rep(c("Jeffreys"),5), rep(c("Jensendifference"),5), rep(c("Chebysev"),5), rep(c("Kulczynskid"),5), 
                                                   rep(c("Harmonicmean"),5), rep(c("Matusita"),5), rep(c("Neyman"),5), rep(c("Clark"),5), 
                                                   rep(c("Kdivergence"),5), rep(c("Taneja"),5), 
                                                   rep(c("Maximum"),5), rep(c("Hamming"),5)),
                                        CS1 = NA, CS2 = NA, CS3 = NA, CS4 = NA,
                                        CS5 = NA, CS6 = NA, CS7 = NA, CS8 = NA,
                                        CS9 = NA, CS10 = NA, CS11 = NA, CS12 = NA,
                                        CS13 = NA, CS14 = NA, CS15 = NA, CS16 = NA,
                                        Bivalent = NA, Active = NA, Divergent = NA, Heterochromatin = NA, Quies = NA)

for(i in 1:21){# 16 CS and 5 CS functional groups
  cat(paste0(i / 21 * 100, ' % CS groups started \n'))
  for(j in 1:5){ # lecif score percentile bins
    cat(paste0(j / 5 * 100, ' % LECIF bins started \n'))
    CS.Simmilarity.LECIF.OsZm[j,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 1)
    CS.Simmilarity.LECIF.OsZm[j+5,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 2) 
    CS.Simmilarity.LECIF.OsZm[j+10,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 3)  
    CS.Simmilarity.LECIF.OsZm[j+15,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 4) 
    CS.Simmilarity.LECIF.OsZm[j+20,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 5) 
    CS.Simmilarity.LECIF.OsZm[j+25,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 6) 
    CS.Simmilarity.LECIF.OsZm[j+30,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 7) 
    CS.Simmilarity.LECIF.OsZm[j+35,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 8) 
    CS.Simmilarity.LECIF.OsZm[j+40,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 9) 
    CS.Simmilarity.LECIF.OsZm[j+45,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 10) 
    CS.Simmilarity.LECIF.OsZm[j+50,2+i] <- philentropy::euclidean(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F) 
    CS.Simmilarity.LECIF.OsZm[j+55,2+i] <- philentropy::sorensen(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+60,2+i] <- philentropy::canberra(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+65,2+i] <- philentropy::wave_hedges(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+70,2+i] <- philentropy::tanimoto(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+75,2+i] <- philentropy::fidelity(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+80,2+i] <- philentropy::squared_chord(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+85,2+i] <- philentropy::squared_chi_sq(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+90,2+i] <- philentropy::additive_symm_chi_sq(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+95,2+i] <- philentropy::topsoe(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F,unit = "log2")
    CS.Simmilarity.LECIF.OsZm[j+100,2+i] <- philentropy::kumar_johnson(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.OsZm[j+105,2+i] <- philentropy::manhattan(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+110,2+i] <- philentropy::lorentzian(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.OsZm[j+115,2+i] <- philentropy::czekanowski(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+120,2+i] <- philentropy::kumar_hassebrook(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+125,2+i] <- philentropy::squared_euclidean(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+130,2+i] <- philentropy::prob_symm_chi_sq(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+135,2+i] <- philentropy::kullback_leibler_distance(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.OsZm[j+140,2+i] <- philentropy::jensen_shannon(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.OsZm[j+145,2+i] <- philentropy::minkowski(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, n = 3)
    CS.Simmilarity.LECIF.OsZm[j+150,2+i] <- philentropy::soergel(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+155,2+i] <- philentropy::inner_product(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+160,2+i] <- philentropy::hellinger(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+165,2+i] <- philentropy::divergence_sq(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+170,2+i] <- philentropy::jeffreys(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.OsZm[j+175,2+i] <- philentropy::jensen_difference(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.OsZm[j+180,2+i] <- philentropy::chebyshev(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+185,2+i] <- philentropy::kulczynski_d(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.OsZm[j+190,2+i] <- philentropy::harmonic_mean_dist(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+195,2+i] <- philentropy::matusita(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+200,2+i] <- philentropy::neyman_chi_sq(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.OsZm[j+205,2+i] <- philentropy::clark_sq(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+210,2+i] <- philentropy::k_divergence(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.OsZm[j+215,2+i] <- philentropy::taneja(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.OsZm[j+220,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = "maximum", margin = 2)[1,2]
    CS.Simmilarity.LECIF.OsZm[j+225,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = "hamming", margin = 2)[1,2]
    cat(paste0(j / 5 * 100, ' % LECIF bins completed \n'))
  }
  cat(paste0(i / 21 * 100, ' % CS groups completed \n'))
}

write.table(OsZm.df.F, "H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/OsZm_CSsimmilarity.DF5.txt", col.names = T, row.names = F, sep = "\t", quote = F)
CS.Simmilarity.LECIF.OsZm.5 <- CS.Simmilarity.LECIF.OsZm
write.table(CS.Simmilarity.LECIF.OsZm, file = "H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/OsZm_CSsimmilarity.DF.Distances5.txt", col.names = T, row.names = F, sep = "\t", quote = F)
View(CS.Simmilarity.LECIF.OsZm.5)

OsZm.df.F$Percentile.Rank.Bin.Three <- NA
OsZm.df.F$Percentile.Rank.Bin.Three[which(OsZm.df.F$Percentile.Rank.Score < 0.333333)] <- "First" # 651876 regions
OsZm.df.F$Percentile.Rank.Bin.Three[which(OsZm.df.F$Percentile.Rank.Score >= 0.333333 & OsZm.df.F$Percentile.Rank.Score < 0.666666)] <- "Second" # 651873
OsZm.df.F$Percentile.Rank.Bin.Three[which(OsZm.df.F$Percentile.Rank.Score >= 0.666666 & OsZm.df.F$Percentile.Rank.Score <= 1)] <- "Thirdth" # 651867
( 651876 + 651873 + 651867) == nrow(OsZm.df.F) # check


CS.Simmilarity.LECIF.OsZm <- data.frame(Bin = rep(c("First", "Second", "Thirdth"),46),
                                        Metric = c(rep(c("Jaccard"),3),rep(c("SMC"),3),rep(c("Sokal"),3), rep(c("Tanimoto"),3),
                                                   rep(c("Dice"),3), rep(c("Hanman"),3), rep(c("Ochiai"),3), rep(c("Sokal2"),3),
                                                   rep(c("Phi"),3), rep(c("S2"),3), 
                                                   rep(c("Eucledian"),3), rep(c("Sorensen"),3), rep(c("Canberra"),3), rep(c("Wavehedges"),3), 
                                                   rep(c("Tanimoto2"),3), rep(c("Fidelity"),3), rep(c("SquaredChord"),3), rep(c("SquaredChi"),3), 
                                                   rep(c("AdditiveSym"),3), rep(c("Topsoe"),3), rep(c("Kumarjhonson"),3), rep(c("Manhattan"),3),
                                                   rep(c("Lorentzian"),3), rep(c("Czekanowski"),3), rep(c("Hassebrook"),3),  rep(c("SquaredEucledian"),3),
                                                   rep(c("ProbSymm"),3), rep(c("Kullbackleiber"),3), rep(c("Jensenshannon"),3), rep(c("Minkowski"),3), 
                                                   rep(c("Soergel"),3), rep(c("Innerproduct"),3), rep(c("Hellinger"),3), rep(c("Divergence"),3), 
                                                   rep(c("Jeffreys"),3), rep(c("Jensendifference"),3), rep(c("Chebysev"),3), rep(c("Kulczynskid"),3), 
                                                   rep(c("Harmonicmean"),3), rep(c("Matusita"),3), rep(c("Neyman"),3), rep(c("Clark"),3), 
                                                   rep(c("Kdivergence"),3), rep(c("Taneja"),3), 
                                                   rep(c("Maximum"),3), rep(c("Hamming"),3)),
                                        CS1 = NA, CS2 = NA, CS3 = NA, CS4 = NA,
                                        CS5 = NA, CS6 = NA, CS7 = NA, CS8 = NA,
                                        CS9 = NA, CS10 = NA, CS11 = NA, CS12 = NA,
                                        CS13 = NA, CS14 = NA, CS15 = NA, CS16 = NA,
                                        Bivalent = NA, Active = NA, Divergent = NA, Heterochromatin = NA, Quies = NA)

for(i in 1:21){# 16 CS and 5 CS functional groups
  cat(paste0(i / 21 * 100, ' % CS groups started \n'))
  for(j in 1:3){ # lecif score percentile bins
    cat(paste0(j / 3 * 100, ' % LECIF bins started \n'))
    CS.Simmilarity.LECIF.OsZm[j,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 1)
    CS.Simmilarity.LECIF.OsZm[j+3,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 2) 
    CS.Simmilarity.LECIF.OsZm[j+6,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 3)  
    CS.Simmilarity.LECIF.OsZm[j+9,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 4) 
    CS.Simmilarity.LECIF.OsZm[j+12,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 5) 
    CS.Simmilarity.LECIF.OsZm[j+15,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 6) 
    CS.Simmilarity.LECIF.OsZm[j+18,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 7) 
    CS.Simmilarity.LECIF.OsZm[j+21,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 8) 
    CS.Simmilarity.LECIF.OsZm[j+24,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 9) 
    CS.Simmilarity.LECIF.OsZm[j+27,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = 10) 
    CS.Simmilarity.LECIF.OsZm[j+30,2+i] <- philentropy::euclidean(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F) 
    CS.Simmilarity.LECIF.OsZm[j+33,2+i] <- philentropy::sorensen(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+36,2+i] <- philentropy::canberra(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+39,2+i] <- philentropy::wave_hedges(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+42,2+i] <- philentropy::tanimoto(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+45,2+i] <- philentropy::fidelity(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+48,2+i] <- philentropy::squared_chord(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+51,2+i] <- philentropy::squared_chi_sq(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+54,2+i] <- philentropy::additive_symm_chi_sq(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+57,2+i] <- philentropy::topsoe(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F,unit = "log2")
    CS.Simmilarity.LECIF.OsZm[j+60,2+i] <- philentropy::kumar_johnson(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.OsZm[j+63,2+i] <- philentropy::manhattan(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+66,2+i] <- philentropy::lorentzian(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.OsZm[j+69,2+i] <- philentropy::czekanowski(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+72,2+i] <- philentropy::kumar_hassebrook(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+75,2+i] <- philentropy::squared_euclidean(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+78,2+i] <- philentropy::prob_symm_chi_sq(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+81,2+i] <- philentropy::kullback_leibler_distance(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.OsZm[j+84,2+i] <- philentropy::jensen_shannon(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.OsZm[j+87,2+i] <- philentropy::minkowski(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, n = 3)
    CS.Simmilarity.LECIF.OsZm[j+90,2+i] <- philentropy::soergel(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+93,2+i] <- philentropy::inner_product(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+96,2+i] <- philentropy::hellinger(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+99,2+i] <- philentropy::divergence_sq(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+102,2+i] <- philentropy::jeffreys(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.OsZm[j+105,2+i] <- philentropy::jensen_difference(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.OsZm[j+108,2+i] <- philentropy::chebyshev(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+111,2+i] <- philentropy::kulczynski_d(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.OsZm[j+114,2+i] <- philentropy::harmonic_mean_dist(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+117,2+i] <- philentropy::matusita(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+120,2+i] <- philentropy::neyman_chi_sq(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.OsZm[j+123,2+i] <- philentropy::clark_sq(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.OsZm[j+126,2+i] <- philentropy::k_divergence(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.OsZm[j+129,2+i] <- philentropy::taneja(P = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], Q = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.OsZm[j+132,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = "maximum", margin = 2)[1,2]
    CS.Simmilarity.LECIF.OsZm[j+135,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Os = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),9+i], CS.Zm = OsZm.df.F[which(OsZm.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.OsZm[j,1]),30+i])), method = "hamming", margin = 2)[1,2]
    cat(paste0(j / 3 * 100, ' % LECIF bins completed \n'))
  }
  cat(paste0(i / 21 * 100, ' % CS groups completed \n'))
}

write.table(OsZm.df.F, "H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/OsZm_CSsimmilarity.DF3.txt", col.names = T, row.names = F, sep = "\t", quote = F)
CS.Simmilarity.LECIF.OsZm.3 <- CS.Simmilarity.LECIF.OsZm
write.table(CS.Simmilarity.LECIF.OsZm, file = "H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/OsZm_CSsimmilarity.DF.Distances3.txt", col.names = T, row.names = F, sep = "\t", quote = F)
View(CS.Simmilarity.LECIF.OsZm.3) ## Dice-like would be the best option

###############################################
########### 1.2.5 Zea - Arabidopsis ########### 
##############################################

AtZm_regions <- read.delim("H:/Taiotactical/Taiotactical_LECIF/example/check_zm-at/all.m", header = F, sep = "\t") 
ZmAt_regions <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/Zm-At_LECIFscore.bedgraph", header = F, sep = "\t") 
ZmAt_bed <- ZmAt_regions[,c(1,2,3)]
AtZm_bed <- AtZm_regions[,c(1,2,3)]
AtZm_bed <- AtZm_bed[grep("Chr", AtZm_bed$V1),]
AtZm_bed$V3 <- as.numeric(AtZm_bed$V3) - 1
AtZm_bed$V2 <- as.numeric(AtZm_bed$V2)
nrow(ZmAt_bed) == nrow(AtZm_bed) # check the same number of at aligning regions to os
AtZm_bed$V3 <- AtZm_bed$V3 + (ZmAt_bed$V3 - ZmAt_bed$V2)
length(which(((ZmAt_bed$V3 - ZmAt_bed$V2) == (AtZm_bed$V3) -  (AtZm_bed$V2)) == TRUE)) # check same expansion of the first aligning base
write.table(AtZm_bed,"H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/At_ZmAligningRegions.bed", quote = F, sep = "\t", col.names = F, row.names = F)
# sort -k1,1 -k2,2n At_ZmAligningRegions.bed > At_ZmAligningRegions.sorted.bed
# bedtools intersect -a At_ZmAligningRegions.sorted.bed -b hihmm.model1.K15.At.ReMapped.ReNamed.Reduced.Reordered.bed -wa -wb | cut -f1,2,3,7 > At_ZmAligningRegions.hiHMM.sorted.bed
ZmAt_hiHMMInfo <- LECIF_ZmAt[,c(1,2,3,5)] 
AtZm_hiHMMInfo <- read.table("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/At_ZmAligningRegions.hiHMM.sorted.bed", header = F, sep = "\t")
ZmAt.df.F <- data.frame(Chr.Zm = ZmAt_bed$V1, Start.Zm = ZmAt_bed$V2, End.Zm = ZmAt_bed$V3,
                        Chr.At = AtZm_bed$V1, Start.At = AtZm_bed$V2, End.At = AtZm_bed$V3,
                        LECIF.Score = ZmAt_regions$V4, Percentile.Rank.Score = percent_rank(x = ZmAt_regions$V4), Percentile.Rank.Bin = NA,
                        CS1.Zm = 0, CS2.Zm = 0, CS3.Zm = 0, CS4.Zm = 0, CS5.Zm = 0, CS6.Zm = 0,
                        CS7.Zm = 0, CS8.Zm = 0, CS9.Zm = 0, CS10.Zm = 0, CS11.Zm = 0, CS12.Zm = 0,
                        CS13.Zm = 0, CS14.Zm = 0, CS15.Zm = 0, CS16.Zm = 0,
                        CS1.At = 0, CS2.At = 0, CS3.At = 0, CS4.At = 0, CS5.At = 0, CS6.At = 0,
                        CS7.At = 0, CS8.At = 0, CS9.At = 0, CS10.At = 0, CS11.At = 0, CS12.At = 0,
                        CS13.At = 0, CS14.At = 0, CS15.At = 0, CS16.At = 0,
                        Bivalent.Zm = 0, Active.Zm = 0, Divergent.Zm = 0, Heterochromatin.Zm = 0, Quies.Zm = 0,
                        Bivalent.At = 0, Active.At = 0, Divergent.At = 0, Heterochromatin.At = 0, Quies.At = 0)
length(percent_rank(x = ZmAt_regions$V4)) == nrow(ZmAt.df.F) # check

ZmAt.df.F$Percentile.Rank.Bin[which(ZmAt.df.F$Percentile.Rank.Score < 0.1)] <- "First" # 134359 regions
ZmAt.df.F$Percentile.Rank.Bin[which(ZmAt.df.F$Percentile.Rank.Score >= 0.1 & ZmAt.df.F$Percentile.Rank.Score < 0.2)] <- "Second" # 134410
ZmAt.df.F$Percentile.Rank.Bin[which(ZmAt.df.F$Percentile.Rank.Score >= 0.2 & ZmAt.df.F$Percentile.Rank.Score < 0.3)] <- "Thirdth" # 134315
ZmAt.df.F$Percentile.Rank.Bin[which(ZmAt.df.F$Percentile.Rank.Score >= 0.3 & ZmAt.df.F$Percentile.Rank.Score < 0.4)] <- "Fourth" # 134401
ZmAt.df.F$Percentile.Rank.Bin[which(ZmAt.df.F$Percentile.Rank.Score >= 0.4 & ZmAt.df.F$Percentile.Rank.Score < 0.5)] <- "Fifth" # 134348
ZmAt.df.F$Percentile.Rank.Bin[which(ZmAt.df.F$Percentile.Rank.Score >= 0.5 & ZmAt.df.F$Percentile.Rank.Score < 0.6)] <- "Sixth" # 134336
ZmAt.df.F$Percentile.Rank.Bin[which(ZmAt.df.F$Percentile.Rank.Score >= 0.6 & ZmAt.df.F$Percentile.Rank.Score < 0.7)] <- "Seventh" # 134347
ZmAt.df.F$Percentile.Rank.Bin[which(ZmAt.df.F$Percentile.Rank.Score >= 0.7 & ZmAt.df.F$Percentile.Rank.Score < 0.8)] <- "Eigth" # 134411
ZmAt.df.F$Percentile.Rank.Bin[which(ZmAt.df.F$Percentile.Rank.Score >= 0.8 & ZmAt.df.F$Percentile.Rank.Score < 0.9)] <- "Nineth" # 134309
ZmAt.df.F$Percentile.Rank.Bin[which(ZmAt.df.F$Percentile.Rank.Score >= 0.9 & ZmAt.df.F$Percentile.Rank.Score <= 1)] <- "Tenth" # 134346
(134359 + 134410 + 134315 + 134401 + 134348 + 134336 + 134347 + 134411 + 134309 + 134346) == nrow(ZmAt.df.F) # check

for(i in 1:nrow(ZmAt.df.F)){
  cat(paste0(i / nrow(ZmAt.df.F) * 100, ' % completed \n'))
  ZmCSinfo <- ZmAt_hiHMMInfo$V5[which(ZmAt_hiHMMInfo$V1 == ZmAt.df.F$Chr.Zm[i] & ZmAt_hiHMMInfo$V2 == ZmAt.df.F$Start.Zm[i] & ZmAt_hiHMMInfo$V3 == ZmAt.df.F$End.Zm[i])]
  AtCSinfo <- AtZm_hiHMMInfo$V4[which(AtZm_hiHMMInfo$V1 == ZmAt.df.F$Chr.At[i] & AtZm_hiHMMInfo$V2 == ZmAt.df.F$Start.At[i] & AtZm_hiHMMInfo$V3 == ZmAt.df.F$End.At[i])]
  ZmAt.df.F[i,which(c(paste0("CS",1:16)) %in% ZmCSinfo)+9] <- 1
  ZmAt.df.F[i,which(c(paste0("CS",1:16)) %in% AtCSinfo)+25] <- 1
  if(length(which(c("CS1", "CS2", "CS3", "CS4", "CS5") %in% ZmCSinfo))>0){
    ZmAt.df.F$Bivalent.Zm[i] <- 1
  }
  if(length(which(c("CS1", "CS2", "CS3", "CS4", "CS5") %in% AtCSinfo))>0){
    ZmAt.df.F$Bivalent.At[i] <- 1
  }
  if(length(which(c("CS6", "CS7", "CS8", "CS9") %in% ZmCSinfo))>0){
    ZmAt.df.F$Active.Zm[i] <- 1
  }
  if(length(which(c("CS6", "CS7", "CS8", "CS9") %in% AtCSinfo))>0){
    ZmAt.df.F$Active.At[i] <- 1
  }
  if(length(which(c("CS10") %in% ZmCSinfo))>0){
    ZmAt.df.F$Divergent.Zm[i] <- 1
  }
  if(length(which(c("CS10") %in% AtCSinfo))>0){
    ZmAt.df.F$Divergent.At[i] <- 1
  }
  if(length(which(c("CS11", "CS12", "CS13") %in% ZmCSinfo))>0){
    ZmAt.df.F$Heterochromatin.Zm[i] <- 1
  }
  if(length(which(c("CS11", "CS12", "CS13") %in% AtCSinfo))>0){
    ZmAt.df.F$Heterochromatin.At[i] <- 1
  }
  if(length(which(c("CS14", "CS15", "CS16") %in% ZmCSinfo))>0){
    ZmAt.df.F$Quies.Zm[i] <- 1
  }
  if(length(which(c("CS14", "CS15", "CS16") %in% AtCSinfo))>0){
    ZmAt.df.F$Quies.At[i] <- 1
  }
}

#write.table(AtOs.df.F, "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/AtOs_CSsimmilarity.DF.txt", col.names = T, row.names = F, sep = "\t", quote = F)
#AtOs.df.F. <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/AtOs_CSsimmilarity.DF.txt",header = T, row.names = F, sep = "\t", quote = F)
# loop to compute simmilarity between this huge binary vectos || metric performing well in boolean data but not too influences by super-sparsed data (low coverage from some states)
CS.Simmilarity.LECIF.ZmAt <- data.frame(Bin = rep(c("First", "Second", "Thirdth", "Fourth", "Fifth", "Sixth", "Seventh", "Eigth", "Nineth", "Tenth"),46),
                                        Metric = c(rep(c("Jaccard"),10),rep(c("SMC"),10),rep(c("Sokal"),10), rep(c("Tanimoto"),10),
                                                   rep(c("Dice"),10), rep(c("Hanman"),10), rep(c("Ochiai"),10), rep(c("Sokal2"),10),
                                                   rep(c("Phi"),10), rep(c("S2"),10), 
                                                   rep(c("Eucledian"),10), rep(c("Sorensen"),10), rep(c("Canberra"),10), rep(c("Wavehedges"),10), 
                                                   rep(c("Tanimoto2"),10), rep(c("Fidelity"),10), rep(c("SquaredChord"),10), rep(c("SquaredChi"),10), 
                                                   rep(c("AdditiveSym"),10), rep(c("Topsoe"),10), rep(c("Kumarjhonson"),10), rep(c("Manhattan"),10),
                                                   rep(c("Lorentzian"),10), rep(c("Czekanowski"),10), rep(c("Hassebrook"),10),  rep(c("SquaredEucledian"),10),
                                                   rep(c("ProbSymm"),10), rep(c("Kullbackleiber"),10), rep(c("Jensenshannon"),10), rep(c("Minkowski"),10), 
                                                   rep(c("Soergel"),10), rep(c("Innerproduct"),10), rep(c("Hellinger"),10), rep(c("Divergence"),10), 
                                                   rep(c("Jeffreys"),10), rep(c("Jensendifference"),10), rep(c("Chebysev"),10), rep(c("Kulczynskid"),10), 
                                                   rep(c("Harmonicmean"),10), rep(c("Matusita"),10), rep(c("Neyman"),10), rep(c("Clark"),10), 
                                                   rep(c("Kdivergence"),10), rep(c("Taneja"),10), 
                                                   rep(c("Maximum"),10), rep(c("Hamming"),10)),
                                        CS1 = NA, CS2 = NA, CS3 = NA, CS4 = NA,
                                        CS5 = NA, CS6 = NA, CS7 = NA, CS8 = NA,
                                        CS9 = NA, CS10 = NA, CS11 = NA, CS12 = NA,
                                        CS13 = NA, CS14 = NA, CS15 = NA, CS16 = NA,
                                        Bivalent = NA, Active = NA, Divergent = NA, Heterochromatin = NA, Quies = NA)
## Reorder prev dataframe to enter the loop
ZmAt.df.F <- ZmAt.df.F[,c(1:25, 42:46, 26:41, 47:51)]
## Loop to fill the distance/simmilarity metrics between CS (states and groups) splited by LECIF percentile bins
for(i in 1:21){# 16 CS and 5 CS functional groups
  cat(paste0(i / 21 * 100, ' % CS groups started \n'))
  for(j in 1:10){ # lecif score percentile bins
    cat(paste0(j / 10 * 100, ' % LECIF bins started \n'))
    CS.Simmilarity.LECIF.ZmAt[j,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 1)
    CS.Simmilarity.LECIF.ZmAt[j+10,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 2) 
    CS.Simmilarity.LECIF.ZmAt[j+20,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 3)  
    CS.Simmilarity.LECIF.ZmAt[j+30,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 4) 
    CS.Simmilarity.LECIF.ZmAt[j+40,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 5) 
    CS.Simmilarity.LECIF.ZmAt[j+50,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 6) 
    CS.Simmilarity.LECIF.ZmAt[j+60,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 7) 
    CS.Simmilarity.LECIF.ZmAt[j+70,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 8) 
    CS.Simmilarity.LECIF.ZmAt[j+80,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 9) 
    CS.Simmilarity.LECIF.ZmAt[j+90,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 10) 
    CS.Simmilarity.LECIF.ZmAt[j+100,2+i] <- philentropy::euclidean(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F) 
    CS.Simmilarity.LECIF.ZmAt[j+110,2+i] <- philentropy::sorensen(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+120,2+i] <- philentropy::canberra(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+130,2+i] <- philentropy::wave_hedges(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+140,2+i] <- philentropy::tanimoto(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+150,2+i] <- philentropy::fidelity(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+160,2+i] <- philentropy::squared_chord(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+170,2+i] <- philentropy::squared_chi_sq(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+180,2+i] <- philentropy::additive_symm_chi_sq(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+190,2+i] <- philentropy::topsoe(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F,unit = "log2")
    CS.Simmilarity.LECIF.ZmAt[j+200,2+i] <- philentropy::kumar_johnson(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmAt[j+210,2+i] <- philentropy::manhattan(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+220,2+i] <- philentropy::lorentzian(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.ZmAt[j+230,2+i] <- philentropy::czekanowski(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+240,2+i] <- philentropy::kumar_hassebrook(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+250,2+i] <- philentropy::squared_euclidean(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+260,2+i] <- philentropy::prob_symm_chi_sq(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+270,2+i] <- philentropy::kullback_leibler_distance(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmAt[j+280,2+i] <- philentropy::jensen_shannon(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.ZmAt[j+290,2+i] <- philentropy::minkowski(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, n = 3)
    CS.Simmilarity.LECIF.ZmAt[j+300,2+i] <- philentropy::soergel(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+310,2+i] <- philentropy::inner_product(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+320,2+i] <- philentropy::hellinger(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+330,2+i] <- philentropy::divergence_sq(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+340,2+i] <- philentropy::jeffreys(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmAt[j+350,2+i] <- philentropy::jensen_difference(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.ZmAt[j+360,2+i] <- philentropy::chebyshev(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+370,2+i] <- philentropy::kulczynski_d(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmAt[j+380,2+i] <- philentropy::harmonic_mean_dist(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+390,2+i] <- philentropy::matusita(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+400,2+i] <- philentropy::neyman_chi_sq(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmAt[j+410,2+i] <- philentropy::clark_sq(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+420,2+i] <- philentropy::k_divergence(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.ZmAt[j+430,2+i] <- philentropy::taneja(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmAt[j+440,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = "maximum", margin = 2)[1,2]
    CS.Simmilarity.LECIF.ZmAt[j+450,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = "hamming", margin = 2)[1,2]
    cat(paste0(j / 10 * 100, ' % LECIF bins completed \n'))
  }
  cat(paste0(i / 21 * 100, ' % CS groups completed \n'))
}

write.table(ZmAt.df.F, "H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/ZmAt_CSsimmilarity.DF10.txt", col.names = T, row.names = F, sep = "\t", quote = F)
CS.Simmilarity.LECIF.ZmAt.10 <- CS.Simmilarity.LECIF.ZmAt
write.table(CS.Simmilarity.LECIF.ZmAt, file = "H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/ZmAt_CSsimmilarity.DF.Distances10.txt", col.names = T, row.names = F, sep = "\t", quote = F)

## Example code to check if metrics are distance or simmilarity and max/min values for each
ade4::dist.binary(t(data.frame(CS1.Zm = ZmAt.df.F3$CS2.At[which(ZmAt.df.F$Percentile.Rank.Bin == "First")], CS1.At = ZmAt.df.F3$CS2.At[which(ZmAt.df.F$Percentile.Rank.Bin == "First")])), method = 5)
ade4::dist.binary(t(data.frame(CS1.Zm = ZmAt.df.F$CS2.At[which(ZmAt.df.F$Percentile.Rank.Bin == "Second")], CS1.At = ZmAt.df.F$CS2.Os[which(ZmAt.df.F$Percentile.Rank.Bin == "Second")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Zm = ZmAt.df.F$CS2.At[which(ZmAt.df.F$Percentile.Rank.Bin == "Thirdth")], CS1.At = ZmAt.df.F$CS2.Os[which(ZmAt.df.F$Percentile.Rank.Bin == "Thirdth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Zm = ZmAt.df.F$CS2.At[which(ZmAt.df.F$Percentile.Rank.Bin == "Fourth")], CS1.At = ZmAt.df.F$CS2.Os[which(ZmAt.df.F$Percentile.Rank.Bin == "Fourth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Zm = ZmAt.df.F$CS2.At[which(ZmAt.df.F$Percentile.Rank.Bin == "Fifth")], CS1.At = ZmAt.df.F$CS2.Os[which(ZmAt.df.F$Percentile.Rank.Bin == "Fifth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Zm = ZmAt.df.F$CS2.At[which(ZmAt.df.F$Percentile.Rank.Bin == "Sixth")], CS1.At = ZmAt.df.F$CS2.Os[which(ZmAt.df.F$Percentile.Rank.Bin == "Sixth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Zm = ZmAt.df.F$CS2.At[which(ZmAt.df.F$Percentile.Rank.Bin == "Seventh")], CS1.At = ZmAt.df.F$CS2.Os[which(ZmAt.df.F$Percentile.Rank.Bin == "Seventh")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Zm = ZmAt.df.F$CS2.At[which(ZmAt.df.F$Percentile.Rank.Bin == "Eigth")], CS1.At = ZmAt.df.F$CS2.Os[which(ZmAt.df.F$Percentile.Rank.Bin == "Eigth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Zm = ZmAt.df.F$CS2.At[which(ZmAt.df.F$Percentile.Rank.Bin == "Nineth")], CS1.At = ZmAt.df.F$CS2.Os[which(ZmAt.df.F$Percentile.Rank.Bin == "Nineth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Zm = ZmAt.df.F$CS2.At[which(ZmAt.df.F$Percentile.Rank.Bin == "Tenth")], CS1.At = ZmAt.df.F$CS2.Os[which(ZmAt.df.F$Percentile.Rank.Bin == "Tenth")])), method = 2)
philentropy::jaccard(P =  ZmAt.df.F3$CS1.At[which(ZmAt.df.F$Percentile.Rank.Bin == "First")], Q = ZmAt.df.F3$CS1.At[which(ZmAt.df.F$Percentile.Rank.Bin == "First")], testNA = F)
proxyC::dist(as.matrix(data.frame(Gene1 = ZmAt.df.F$CS1.At[which(ZmAt.df.F$Percentile.Rank.Bin == "First")], Gene2 = ZmAt.df.F$CS1.At[which(ZmAt.df.F$Percentile.Rank.Bin == "First")])), method = "hamming", margin = 2)[1,2]

### with less bins (5 or 3)

ZmAt.df.F$Percentile.Rank.Bin.Five <- NA
ZmAt.df.F$Percentile.Rank.Bin.Five[which(ZmAt.df.F$Percentile.Rank.Score < 0.2)] <- "First" # 268769 regions
ZmAt.df.F$Percentile.Rank.Bin.Five[which(ZmAt.df.F$Percentile.Rank.Score >= 0.2 & ZmAt.df.F$Percentile.Rank.Score < 0.4)] <- "Second" # 268716
ZmAt.df.F$Percentile.Rank.Bin.Five[which(ZmAt.df.F$Percentile.Rank.Score >= 0.4 & ZmAt.df.F$Percentile.Rank.Score < 0.6)] <- "Thirdth" # 268684
ZmAt.df.F$Percentile.Rank.Bin.Five[which(ZmAt.df.F$Percentile.Rank.Score >= 0.6 & ZmAt.df.F$Percentile.Rank.Score < 0.8)] <- "Fourth" # 268758
ZmAt.df.F$Percentile.Rank.Bin.Five[which(ZmAt.df.F$Percentile.Rank.Score >= 0.8 & ZmAt.df.F$Percentile.Rank.Score <= 1)] <- "Fifth" # 268655
(268769 + 268716 + 268684 + 268758 + 268655) == nrow(ZmAt.df.F) # check


CS.Simmilarity.LECIF.ZmAt <- data.frame(Bin = rep(c("First", "Second", "Thirdth", "Fourth", "Fifth"),46),
                                        Metric = c(rep(c("Jaccard"),5),rep(c("SMC"),5),rep(c("Sokal"),5), rep(c("Tanimoto"),5),
                                                   rep(c("Dice"),5), rep(c("Hanman"),5), rep(c("Ochiai"),5), rep(c("Sokal2"),5),
                                                   rep(c("Phi"),5), rep(c("S2"),5), 
                                                   rep(c("Eucledian"),5), rep(c("Sorensen"),5), rep(c("Canberra"),5), rep(c("Wavehedges"),5), 
                                                   rep(c("Tanimoto2"),5), rep(c("Fidelity"),5), rep(c("SquaredChord"),5), rep(c("SquaredChi"),5), 
                                                   rep(c("AdditiveSym"),5), rep(c("Topsoe"),5), rep(c("Kumarjhonson"),5), rep(c("Manhattan"),5),
                                                   rep(c("Lorentzian"),5), rep(c("Czekanowski"),5), rep(c("Hassebrook"),5),  rep(c("SquaredEucledian"),5),
                                                   rep(c("ProbSymm"),5), rep(c("Kullbackleiber"),5), rep(c("Jensenshannon"),5), rep(c("Minkowski"),5), 
                                                   rep(c("Soergel"),5), rep(c("Innerproduct"),5), rep(c("Hellinger"),5), rep(c("Divergence"),5), 
                                                   rep(c("Jeffreys"),5), rep(c("Jensendifference"),5), rep(c("Chebysev"),5), rep(c("Kulczynskid"),5), 
                                                   rep(c("Harmonicmean"),5), rep(c("Matusita"),5), rep(c("Neyman"),5), rep(c("Clark"),5), 
                                                   rep(c("Kdivergence"),5), rep(c("Taneja"),5), 
                                                   rep(c("Maximum"),5), rep(c("Hamming"),5)),
                                        CS1 = NA, CS2 = NA, CS3 = NA, CS4 = NA,
                                        CS5 = NA, CS6 = NA, CS7 = NA, CS8 = NA,
                                        CS9 = NA, CS10 = NA, CS11 = NA, CS12 = NA,
                                        CS13 = NA, CS14 = NA, CS15 = NA, CS16 = NA,
                                        Bivalent = NA, Active = NA, Divergent = NA, Heterochromatin = NA, Quies = NA)

for(i in 1:21){# 16 CS and 5 CS functional groups
  cat(paste0(i / 21 * 100, ' % CS groups started \n'))
  for(j in 1:5){ # lecif score percentile bins
    cat(paste0(j / 5 * 100, ' % LECIF bins started \n'))
    CS.Simmilarity.LECIF.ZmAt[j,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 1)
    CS.Simmilarity.LECIF.ZmAt[j+5,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 2) 
    CS.Simmilarity.LECIF.ZmAt[j+10,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 3)  
    CS.Simmilarity.LECIF.ZmAt[j+15,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 4) 
    CS.Simmilarity.LECIF.ZmAt[j+20,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 5) 
    CS.Simmilarity.LECIF.ZmAt[j+25,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 6) 
    CS.Simmilarity.LECIF.ZmAt[j+30,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 7) 
    CS.Simmilarity.LECIF.ZmAt[j+35,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 8) 
    CS.Simmilarity.LECIF.ZmAt[j+40,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 9) 
    CS.Simmilarity.LECIF.ZmAt[j+45,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 10) 
    CS.Simmilarity.LECIF.ZmAt[j+50,2+i] <- philentropy::euclidean(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F) 
    CS.Simmilarity.LECIF.ZmAt[j+55,2+i] <- philentropy::sorensen(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+60,2+i] <- philentropy::canberra(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+65,2+i] <- philentropy::wave_hedges(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+70,2+i] <- philentropy::tanimoto(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+75,2+i] <- philentropy::fidelity(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+80,2+i] <- philentropy::squared_chord(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+85,2+i] <- philentropy::squared_chi_sq(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+90,2+i] <- philentropy::additive_symm_chi_sq(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+95,2+i] <- philentropy::topsoe(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F,unit = "log2")
    CS.Simmilarity.LECIF.ZmAt[j+100,2+i] <- philentropy::kumar_johnson(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmAt[j+105,2+i] <- philentropy::manhattan(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+110,2+i] <- philentropy::lorentzian(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.ZmAt[j+115,2+i] <- philentropy::czekanowski(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+120,2+i] <- philentropy::kumar_hassebrook(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+125,2+i] <- philentropy::squared_euclidean(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+130,2+i] <- philentropy::prob_symm_chi_sq(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+135,2+i] <- philentropy::kullback_leibler_distance(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmAt[j+140,2+i] <- philentropy::jensen_shannon(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.ZmAt[j+145,2+i] <- philentropy::minkowski(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, n = 3)
    CS.Simmilarity.LECIF.ZmAt[j+150,2+i] <- philentropy::soergel(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+155,2+i] <- philentropy::inner_product(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+160,2+i] <- philentropy::hellinger(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+165,2+i] <- philentropy::divergence_sq(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+170,2+i] <- philentropy::jeffreys(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmAt[j+175,2+i] <- philentropy::jensen_difference(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.ZmAt[j+180,2+i] <- philentropy::chebyshev(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+185,2+i] <- philentropy::kulczynski_d(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmAt[j+190,2+i] <- philentropy::harmonic_mean_dist(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+195,2+i] <- philentropy::matusita(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+200,2+i] <- philentropy::neyman_chi_sq(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmAt[j+205,2+i] <- philentropy::clark_sq(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+210,2+i] <- philentropy::k_divergence(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.ZmAt[j+215,2+i] <- philentropy::taneja(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmAt[j+220,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = "maximum", margin = 2)[1,2]
    CS.Simmilarity.LECIF.ZmAt[j+225,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = "hamming", margin = 2)[1,2]
    cat(paste0(j / 5 * 100, ' % LECIF bins completed \n'))
  }
  cat(paste0(i / 21 * 100, ' % CS groups completed \n'))
}

write.table(ZmAt.df.F, "H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/ZmAt_CSsimmilarity.DF5.txt", col.names = T, row.names = F, sep = "\t", quote = F)
CS.Simmilarity.LECIF.ZmAt.5 <- CS.Simmilarity.LECIF.ZmAt
write.table(CS.Simmilarity.LECIF.ZmAt, file = "H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/ZmAt_CSsimmilarity.DF.Distances5.txt", col.names = T, row.names = F, sep = "\t", quote = F)
View(CS.Simmilarity.LECIF.ZmAt.5)

ZmAt.df.F$Percentile.Rank.Bin.Three <- NA
ZmAt.df.F$Percentile.Rank.Bin.Three[which(ZmAt.df.F$Percentile.Rank.Score < 0.333333)] <- "First" # 447976 regions
ZmAt.df.F$Percentile.Rank.Bin.Three[which(ZmAt.df.F$Percentile.Rank.Score >= 0.333333 & ZmAt.df.F$Percentile.Rank.Score < 0.666666)] <- "Second" # 447790
ZmAt.df.F$Percentile.Rank.Bin.Three[which(ZmAt.df.F$Percentile.Rank.Score >= 0.666666 & ZmAt.df.F$Percentile.Rank.Score <= 1)] <- "Thirdth" # 447816
( 447976 + 447790 + 447816) == nrow(ZmAt.df.F) # check


CS.Simmilarity.LECIF.ZmAt <- data.frame(Bin = rep(c("First", "Second", "Thirdth"),46),
                                        Metric = c(rep(c("Jaccard"),3),rep(c("SMC"),3),rep(c("Sokal"),3), rep(c("Tanimoto"),3),
                                                   rep(c("Dice"),3), rep(c("Hanman"),3), rep(c("Ochiai"),3), rep(c("Sokal2"),3),
                                                   rep(c("Phi"),3), rep(c("S2"),3), 
                                                   rep(c("Eucledian"),3), rep(c("Sorensen"),3), rep(c("Canberra"),3), rep(c("Wavehedges"),3), 
                                                   rep(c("Tanimoto2"),3), rep(c("Fidelity"),3), rep(c("SquaredChord"),3), rep(c("SquaredChi"),3), 
                                                   rep(c("AdditiveSym"),3), rep(c("Topsoe"),3), rep(c("Kumarjhonson"),3), rep(c("Manhattan"),3),
                                                   rep(c("Lorentzian"),3), rep(c("Czekanowski"),3), rep(c("Hassebrook"),3),  rep(c("SquaredEucledian"),3),
                                                   rep(c("ProbSymm"),3), rep(c("Kullbackleiber"),3), rep(c("Jensenshannon"),3), rep(c("Minkowski"),3), 
                                                   rep(c("Soergel"),3), rep(c("Innerproduct"),3), rep(c("Hellinger"),3), rep(c("Divergence"),3), 
                                                   rep(c("Jeffreys"),3), rep(c("Jensendifference"),3), rep(c("Chebysev"),3), rep(c("Kulczynskid"),3), 
                                                   rep(c("Harmonicmean"),3), rep(c("Matusita"),3), rep(c("Neyman"),3), rep(c("Clark"),3), 
                                                   rep(c("Kdivergence"),3), rep(c("Taneja"),3), 
                                                   rep(c("Maximum"),3), rep(c("Hamming"),3)),
                                        CS1 = NA, CS2 = NA, CS3 = NA, CS4 = NA,
                                        CS5 = NA, CS6 = NA, CS7 = NA, CS8 = NA,
                                        CS9 = NA, CS10 = NA, CS11 = NA, CS12 = NA,
                                        CS13 = NA, CS14 = NA, CS15 = NA, CS16 = NA,
                                        Bivalent = NA, Active = NA, Divergent = NA, Heterochromatin = NA, Quies = NA)

for(i in 1:21){# 16 CS and 5 CS functional groups
  cat(paste0(i / 21 * 100, ' % CS groups started \n'))
  for(j in 1:3){ # lecif score percentile bins
    cat(paste0(j / 3 * 100, ' % LECIF bins started \n'))
    CS.Simmilarity.LECIF.ZmAt[j,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 1)
    CS.Simmilarity.LECIF.ZmAt[j+3,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 2) 
    CS.Simmilarity.LECIF.ZmAt[j+6,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 3)  
    CS.Simmilarity.LECIF.ZmAt[j+9,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 4) 
    CS.Simmilarity.LECIF.ZmAt[j+12,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 5) 
    CS.Simmilarity.LECIF.ZmAt[j+15,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 6) 
    CS.Simmilarity.LECIF.ZmAt[j+18,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 7) 
    CS.Simmilarity.LECIF.ZmAt[j+21,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 8) 
    CS.Simmilarity.LECIF.ZmAt[j+24,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 9) 
    CS.Simmilarity.LECIF.ZmAt[j+27,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = 10) 
    CS.Simmilarity.LECIF.ZmAt[j+30,2+i] <- philentropy::euclidean(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F) 
    CS.Simmilarity.LECIF.ZmAt[j+33,2+i] <- philentropy::sorensen(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+36,2+i] <- philentropy::canberra(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+39,2+i] <- philentropy::wave_hedges(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+42,2+i] <- philentropy::tanimoto(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+45,2+i] <- philentropy::fidelity(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+48,2+i] <- philentropy::squared_chord(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+51,2+i] <- philentropy::squared_chi_sq(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+54,2+i] <- philentropy::additive_symm_chi_sq(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+57,2+i] <- philentropy::topsoe(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F,unit = "log2")
    CS.Simmilarity.LECIF.ZmAt[j+60,2+i] <- philentropy::kumar_johnson(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmAt[j+63,2+i] <- philentropy::manhattan(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+66,2+i] <- philentropy::lorentzian(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.ZmAt[j+69,2+i] <- philentropy::czekanowski(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+72,2+i] <- philentropy::kumar_hassebrook(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+75,2+i] <- philentropy::squared_euclidean(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+78,2+i] <- philentropy::prob_symm_chi_sq(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+81,2+i] <- philentropy::kullback_leibler_distance(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmAt[j+84,2+i] <- philentropy::jensen_shannon(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.ZmAt[j+87,2+i] <- philentropy::minkowski(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, n = 3)
    CS.Simmilarity.LECIF.ZmAt[j+90,2+i] <- philentropy::soergel(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+93,2+i] <- philentropy::inner_product(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+96,2+i] <- philentropy::hellinger(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+99,2+i] <- philentropy::divergence_sq(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+102,2+i] <- philentropy::jeffreys(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmAt[j+105,2+i] <- philentropy::jensen_difference(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.ZmAt[j+108,2+i] <- philentropy::chebyshev(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+111,2+i] <- philentropy::kulczynski_d(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmAt[j+114,2+i] <- philentropy::harmonic_mean_dist(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+117,2+i] <- philentropy::matusita(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+120,2+i] <- philentropy::neyman_chi_sq(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmAt[j+123,2+i] <- philentropy::clark_sq(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmAt[j+126,2+i] <- philentropy::k_divergence(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.ZmAt[j+129,2+i] <- philentropy::taneja(P = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], Q = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmAt[j+132,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = "maximum", margin = 2)[1,2]
    CS.Simmilarity.LECIF.ZmAt[j+135,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Zm = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),9+i], CS.At = ZmAt.df.F[which(ZmAt.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmAt[j,1]),30+i])), method = "hamming", margin = 2)[1,2]
    cat(paste0(j / 3 * 100, ' % LECIF bins completed \n'))
  }
  cat(paste0(i / 21 * 100, ' % CS groups completed \n'))
}

write.table(ZmAt.df.F, "H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/ZmAt_CSsimmilarity.DF3.txt", col.names = T, row.names = F, sep = "\t", quote = F)
CS.Simmilarity.LECIF.ZmAt.3 <- CS.Simmilarity.LECIF.ZmAt
write.table(CS.Simmilarity.LECIF.ZmAt, file = "H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/ZmAt_CSsimmilarity.DF.Distances3.txt", col.names = T, row.names = F, sep = "\t", quote = F)
View(CS.Simmilarity.LECIF.ZmAt.3) ## Dice-like would be the best option

###############################################
############# 1.2.6 Zea - Oryza ############## 
##############################################

OsZm_regions <- read.delim("H:/Taiotactical/Taiotactical_LECIF/example/check_zm-os/all.m", header = F, sep = "\t") 
ZmOs_regions <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/Zm-Os_LECIFscore.bedgraph", header = F, sep = "\t") 
ZmOs_bed <- ZmOs_regions[,c(1,2,3)]
OsZm_bed <- OsZm_regions[,c(1,2,3)]
OsZm_bed <- OsZm_bed[grep("Chr", OsZm_bed$V1),]
OsZm_bed$V3 <- as.numeric(OsZm_bed$V3) - 1
OsZm_bed$V2 <- as.numeric(OsZm_bed$V2)
nrow(ZmOs_bed) == nrow(OsZm_bed) # check the same number of at aligning regions to Zm
OsZm_bed$V3 <- OsZm_bed$V3 + (ZmOs_bed$V3 - ZmOs_bed$V2)
length(which(((ZmOs_bed$V3 - ZmOs_bed$V2) == (OsZm_bed$V3) -  (OsZm_bed$V2)) == TRUE)) # check same expansion of the first aligning base
write.table(OsZm_bed,"H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/Os_ZmAligningRegions.bed", quote = F, sep = "\t", col.names = F, row.names = F)
# sort -k1,1 -k2,2n Os_ZmAligningRegions.bed > Os_ZmAligningRegions.sorted.bed
# bedtools intersect -a Os_ZmAligningRegions.sorted.bed -b hihmm.model1.K15.Os.ReMapped.ReNamed.Reduced.Reordered.bed -wa -wb | cut -f1,2,3,7 > Os_ZmAligningRegions.hiHMM.sorted.bed
ZmOs_hiHMMInfo <- LECIF_ZmOs[,c(1,2,3,5)] 
OsZm_hiHMMInfo <- read.table("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/Os_ZmAligningRegions.hiHMM.sorted.bed", header = F, sep = "\t")
ZmOs.df.F <- data.frame(Chr.Zm = ZmOs_bed$V1, Start.Zm = ZmOs_bed$V2, End.Zm = ZmOs_bed$V3,
                        Chr.Os = OsZm_bed$V1, Start.Os = OsZm_bed$V2, End.Os = OsZm_bed$V3,
                        LECIF.Score = ZmOs_regions$V4, Percentile.Rank.Score = percent_rank(x = ZmOs_regions$V4), Percentile.Rank.Bin = NA,
                        CS1.Zm = 0, CS2.Zm = 0, CS3.Zm = 0, CS4.Zm = 0, CS5.Zm = 0, CS6.Zm = 0,
                        CS7.Zm = 0, CS8.Zm = 0, CS9.Zm = 0, CS10.Zm = 0, CS11.Zm = 0, CS12.Zm = 0,
                        CS13.Zm = 0, CS14.Zm = 0, CS15.Zm = 0, CS16.Zm = 0,
                        CS1.Os = 0, CS2.Os = 0, CS3.Os = 0, CS4.Os = 0, CS5.Os = 0, CS6.Os = 0,
                        CS7.Os = 0, CS8.Os = 0, CS9.Os = 0, CS10.Os = 0, CS11.Os = 0, CS12.Os = 0,
                        CS13.Os = 0, CS14.Os = 0, CS15.Os = 0, CS16.Os = 0,
                        Bivalent.Zm = 0, Active.Zm = 0, Divergent.Zm = 0, Heterochromatin.Zm = 0, Quies.Zm = 0,
                        Bivalent.Os = 0, Active.Os = 0, Divergent.Os = 0, Heterochromatin.Os = 0, Quies.Os = 0)
length(percent_rank(x = ZmOs_regions$V4)) == nrow(ZmOs.df.F) # check

ZmOs.df.F$Percentile.Rank.Bin[which(ZmOs.df.F$Percentile.Rank.Score < 0.1)] <- "First" # 330426 regions
ZmOs.df.F$Percentile.Rank.Bin[which(ZmOs.df.F$Percentile.Rank.Score >= 0.1 & ZmOs.df.F$Percentile.Rank.Score < 0.2)] <- "Second" # 330090
ZmOs.df.F$Percentile.Rank.Bin[which(ZmOs.df.F$Percentile.Rank.Score >= 0.2 & ZmOs.df.F$Percentile.Rank.Score < 0.3)] <- "Thirdth" # 330692
ZmOs.df.F$Percentile.Rank.Bin[which(ZmOs.df.F$Percentile.Rank.Score >= 0.3 & ZmOs.df.F$Percentile.Rank.Score < 0.4)] <- "Fourth" # 329765
ZmOs.df.F$Percentile.Rank.Bin[which(ZmOs.df.F$Percentile.Rank.Score >= 0.4 & ZmOs.df.F$Percentile.Rank.Score < 0.5)] <- "Fifth" # 330674
ZmOs.df.F$Percentile.Rank.Bin[which(ZmOs.df.F$Percentile.Rank.Score >= 0.5 & ZmOs.df.F$Percentile.Rank.Score < 0.6)] <- "Sixth" # 329892
ZmOs.df.F$Percentile.Rank.Bin[which(ZmOs.df.F$Percentile.Rank.Score >= 0.6 & ZmOs.df.F$Percentile.Rank.Score < 0.7)] <- "Seventh" # 330140
ZmOs.df.F$Percentile.Rank.Bin[which(ZmOs.df.F$Percentile.Rank.Score >= 0.7 & ZmOs.df.F$Percentile.Rank.Score < 0.8)] <- "Eigth" # 330727
ZmOs.df.F$Percentile.Rank.Bin[which(ZmOs.df.F$Percentile.Rank.Score >= 0.8 & ZmOs.df.F$Percentile.Rank.Score < 0.9)] <- "Nineth" # 329752
ZmOs.df.F$Percentile.Rank.Bin[which(ZmOs.df.F$Percentile.Rank.Score >= 0.9 & ZmOs.df.F$Percentile.Rank.Score <= 1)] <- "Tenth" # 330231
(330426 + 330090 + 330692 + 329765 + 330674 + 329892 + 330140 + 330727 + 329752 + 330231) == nrow(ZmOs.df.F) # check

for(i in 1:nrow(ZmOs.df.F)){
  cat(paste0(i / nrow(ZmOs.df.F) * 100, ' % completed \n'))
  ZmCSinfo <- ZmOs_hiHMMInfo$V5[which(ZmOs_hiHMMInfo$V1 == ZmOs.df.F$Chr.Zm[i] & ZmOs_hiHMMInfo$V2 == ZmOs.df.F$Start.Zm[i] & ZmOs_hiHMMInfo$V3 == ZmOs.df.F$End.Zm[i])]
  OsCSinfo <- OsZm_hiHMMInfo$V4[which(OsZm_hiHMMInfo$V1 == ZmOs.df.F$Chr.Os[i] & OsZm_hiHMMInfo$V2 == ZmOs.df.F$Start.Os[i] & OsZm_hiHMMInfo$V3 == ZmOs.df.F$End.Os[i])]
  ZmOs.df.F[i,which(c(paste0("CS",1:16)) %in% ZmCSinfo)+9] <- 1
  ZmOs.df.F[i,which(c(paste0("CS",1:16)) %in% OsCSinfo)+25] <- 1
  if(length(which(c("CS1", "CS2", "CS3", "CS4", "CS5") %in% ZmCSinfo))>0){
    ZmOs.df.F$Bivalent.Zm[i] <- 1
  }
  if(length(which(c("CS1", "CS2", "CS3", "CS4", "CS5") %in% OsCSinfo))>0){
    ZmOs.df.F$Bivalent.Os[i] <- 1
  }
  if(length(which(c("CS6", "CS7", "CS8", "CS9") %in% ZmCSinfo))>0){
    ZmOs.df.F$Active.Zm[i] <- 1
  }
  if(length(which(c("CS6", "CS7", "CS8", "CS9") %in% OsCSinfo))>0){
    ZmOs.df.F$Active.Os[i] <- 1
  }
  if(length(which(c("CS10") %in% ZmCSinfo))>0){
    ZmOs.df.F$Divergent.Zm[i] <- 1
  }
  if(length(which(c("CS10") %in% OsCSinfo))>0){
    ZmOs.df.F$Divergent.Os[i] <- 1
  }
  if(length(which(c("CS11", "CS12", "CS13") %in% ZmCSinfo))>0){
    ZmOs.df.F$Heterochromatin.Zm[i] <- 1
  }
  if(length(which(c("CS11", "CS12", "CS13") %in% OsCSinfo))>0){
    ZmOs.df.F$Heterochromatin.Os[i] <- 1
  }
  if(length(which(c("CS14", "CS15", "CS16") %in% ZmCSinfo))>0){
    ZmOs.df.F$Quies.Zm[i] <- 1
  }
  if(length(which(c("CS14", "CS15", "CS16") %in% OsCSinfo))>0){
    ZmOs.df.F$Quies.Os[i] <- 1
  }
}

#write.table(ZmOs.df.F, "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/ZmOs_CSsimmilarity.DF.txt", col.names = T, row.names = F, sep = "\t", quote = F)
#ZmOs.df.F. <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/ZmOs_CSsimmilarity.DF.txt",header = T, row.names = F, sep = "\t", quote = F)
# loop to compute simmilarity between this huge binary vectos || metric performing well in boolean data but not too influences by super-sparsed data (low coverage from some states depending on the socre bin)
CS.Simmilarity.LECIF.ZmOs <- data.frame(Bin = rep(c("First", "Second", "Thirdth", "Fourth", "Fifth", "Sixth", "Seventh", "Eigth", "Nineth", "Tenth"),46),
                                        Metric = c(rep(c("Jaccard"),10),rep(c("SMC"),10),rep(c("Sokal"),10), rep(c("Tanimoto"),10),
                                                   rep(c("Dice"),10), rep(c("Hanman"),10), rep(c("Ochiai"),10), rep(c("Sokal2"),10),
                                                   rep(c("Phi"),10), rep(c("S2"),10), 
                                                   rep(c("Eucledian"),10), rep(c("Sorensen"),10), rep(c("Canberra"),10), rep(c("Wavehedges"),10), 
                                                   rep(c("Tanimoto2"),10), rep(c("Fidelity"),10), rep(c("SquaredChord"),10), rep(c("SquaredChi"),10), 
                                                   rep(c("AdditiveSym"),10), rep(c("Topsoe"),10), rep(c("Kumarjhonson"),10), rep(c("Manhattan"),10),
                                                   rep(c("Lorentzian"),10), rep(c("Czekanowski"),10), rep(c("Hassebrook"),10),  rep(c("SquaredEucledian"),10),
                                                   rep(c("ProbSymm"),10), rep(c("Kullbackleiber"),10), rep(c("Jensenshannon"),10), rep(c("Minkowski"),10), 
                                                   rep(c("Soergel"),10), rep(c("Innerproduct"),10), rep(c("Hellinger"),10), rep(c("Divergence"),10), 
                                                   rep(c("Jeffreys"),10), rep(c("Jensendifference"),10), rep(c("Chebysev"),10), rep(c("Kulczynskid"),10), 
                                                   rep(c("Harmonicmean"),10), rep(c("Matusita"),10), rep(c("Neyman"),10), rep(c("Clark"),10), 
                                                   rep(c("Kdivergence"),10), rep(c("Taneja"),10), 
                                                   rep(c("Maximum"),10), rep(c("Hamming"),10)),
                                        CS1 = NA, CS2 = NA, CS3 = NA, CS4 = NA,
                                        CS5 = NA, CS6 = NA, CS7 = NA, CS8 = NA,
                                        CS9 = NA, CS10 = NA, CS11 = NA, CS12 = NA,
                                        CS13 = NA, CS14 = NA, CS15 = NA, CS16 = NA,
                                        Bivalent = NA, Active = NA, Divergent = NA, Heterochromatin = NA, Quies = NA)
## Reorder prev dataframe to enter the loop
ZmOs.df.F <- ZmOs.df.F[,c(1:25, 42:46, 26:41, 47:51)]
## Loop to fill the distance/simmilarity metrics between CS (states and groups) splitted by LECIF percentile bins
for(i in 1:21){# 16 CS and 5 CS functional groups
  cat(paste0(i / 21 * 100, ' % CS groups started \n'))
  for(j in 1:10){ # lecif score percentile bins
    cat(paste0(j / 10 * 100, ' % LECIF bins started \n'))
    CS.Simmilarity.LECIF.ZmOs[j,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 1)
    CS.Simmilarity.LECIF.ZmOs[j+10,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 2) 
    CS.Simmilarity.LECIF.ZmOs[j+20,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 3)  
    CS.Simmilarity.LECIF.ZmOs[j+30,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 4) 
    CS.Simmilarity.LECIF.ZmOs[j+40,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 5) 
    CS.Simmilarity.LECIF.ZmOs[j+50,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 6) 
    CS.Simmilarity.LECIF.ZmOs[j+60,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 7) 
    CS.Simmilarity.LECIF.ZmOs[j+70,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 8) 
    CS.Simmilarity.LECIF.ZmOs[j+80,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 9) 
    CS.Simmilarity.LECIF.ZmOs[j+90,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 10) 
    CS.Simmilarity.LECIF.ZmOs[j+100,2+i] <- philentropy::euclidean(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F) 
    CS.Simmilarity.LECIF.ZmOs[j+110,2+i] <- philentropy::sorensen(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+120,2+i] <- philentropy::canberra(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+130,2+i] <- philentropy::wave_hedges(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+140,2+i] <- philentropy::tanimoto(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+150,2+i] <- philentropy::fidelity(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+160,2+i] <- philentropy::squared_chord(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+170,2+i] <- philentropy::squared_chi_sq(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+180,2+i] <- philentropy::additive_symm_chi_sq(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+190,2+i] <- philentropy::topsoe(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F,unit = "log2")
    CS.Simmilarity.LECIF.ZmOs[j+200,2+i] <- philentropy::kumar_johnson(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmOs[j+210,2+i] <- philentropy::manhattan(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+220,2+i] <- philentropy::lorentzian(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.ZmOs[j+230,2+i] <- philentropy::czekanowski(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+240,2+i] <- philentropy::kumar_hassebrook(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+250,2+i] <- philentropy::squared_euclidean(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+260,2+i] <- philentropy::prob_symm_chi_sq(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+270,2+i] <- philentropy::kullback_leibler_distance(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmOs[j+280,2+i] <- philentropy::jensen_shannon(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.ZmOs[j+290,2+i] <- philentropy::minkowski(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, n = 3)
    CS.Simmilarity.LECIF.ZmOs[j+300,2+i] <- philentropy::soergel(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+310,2+i] <- philentropy::inner_product(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+320,2+i] <- philentropy::hellinger(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+330,2+i] <- philentropy::divergence_sq(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+340,2+i] <- philentropy::jeffreys(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmOs[j+350,2+i] <- philentropy::jensen_difference(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.ZmOs[j+360,2+i] <- philentropy::chebyshev(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+370,2+i] <- philentropy::kulczynski_d(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmOs[j+380,2+i] <- philentropy::harmonic_mean_dist(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+390,2+i] <- philentropy::matusita(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+400,2+i] <- philentropy::neyman_chi_sq(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmOs[j+410,2+i] <- philentropy::clark_sq(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+420,2+i] <- philentropy::k_divergence(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.ZmOs[j+430,2+i] <- philentropy::taneja(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmOs[j+440,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = "maximum", margin = 2)[1,2]
    CS.Simmilarity.LECIF.ZmOs[j+450,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = "hamming", margin = 2)[1,2]
    cat(paste0(j / 10 * 100, ' % LECIF bins completed \n'))
  }
  cat(paste0(i / 21 * 100, ' % CS groups completed \n'))
}

write.table(ZmOs.df.F, "H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/ZmOs_CSsimmilarity.DF10.txt", col.names = T, row.names = F, sep = "\t", quote = F)
CS.Simmilarity.LECIF.ZmOs.10 <- CS.Simmilarity.LECIF.ZmOs
write.table(CS.Simmilarity.LECIF.ZmOs, file = "H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/ZmOs_CSsimmilarity.DF.Distances10.txt", col.names = T, row.names = F, sep = "\t", quote = F)

## Example code to check if metrics are distance or simmilarity and max/min values for each
ade4::dist.binary(t(data.frame(CS1.Zm = ZmOs.df.F3$CS2.Zm[which(ZmOs.df.F$Percentile.Rank.Bin == "First")], CS1.Os = ZmOs.df.F3$CS2.Zm[which(ZmOs.df.F$Percentile.Rank.Bin == "First")])), method = 5)
ade4::dist.binary(t(data.frame(CS1.Zm = ZmOs.df.F$CS2.At[which(ZmOs.df.F$Percentile.Rank.Bin == "Second")], CS1.Os = ZmOs.df.F$CS2.Zm[which(ZmOs.df.F$Percentile.Rank.Bin == "Second")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Zm = ZmOs.df.F$CS2.At[which(ZmOs.df.F$Percentile.Rank.Bin == "Thirdth")], CS1.Os = ZmOs.df.F$CS2.Zm[which(ZmOs.df.F$Percentile.Rank.Bin == "Thirdth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Zm = ZmOs.df.F$CS2.At[which(ZmOs.df.F$Percentile.Rank.Bin == "Fourth")], CS1.Os = ZmOs.df.F$CS2.Zm[which(ZmOs.df.F$Percentile.Rank.Bin == "Fourth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Zm = ZmOs.df.F$CS2.At[which(ZmOs.df.F$Percentile.Rank.Bin == "Fifth")], CS1.Os = ZmOs.df.F$CS2.Zm[which(ZmOs.df.F$Percentile.Rank.Bin == "Fifth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Zm = ZmOs.df.F$CS2.At[which(ZmOs.df.F$Percentile.Rank.Bin == "Sixth")], CS1.Os = ZmOs.df.F$CS2.Zm[which(ZmOs.df.F$Percentile.Rank.Bin == "Sixth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Zm = ZmOs.df.F$CS2.At[which(ZmOs.df.F$Percentile.Rank.Bin == "Seventh")], CS1.Os = ZmOs.df.F$CS2.Zm[which(ZmOs.df.F$Percentile.Rank.Bin == "Seventh")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Zm = ZmOs.df.F$CS2.At[which(ZmOs.df.F$Percentile.Rank.Bin == "Eigth")], CS1.Os = ZmOs.df.F$CS2.Zm[which(ZmOs.df.F$Percentile.Rank.Bin == "Eigth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Zm = ZmOs.df.F$CS2.At[which(ZmOs.df.F$Percentile.Rank.Bin == "Nineth")], CS1.Os = ZmOs.df.F$CS2.Zm[which(ZmOs.df.F$Percentile.Rank.Bin == "Nineth")])), method = 2)
ade4::dist.binary(t(data.frame(CS1.Zm = ZmOs.df.F$CS2.At[which(ZmOs.df.F$Percentile.Rank.Bin == "Tenth")], CS1.Os = ZmOs.df.F$CS2.Zm[which(ZmOs.df.F$Percentile.Rank.Bin == "Tenth")])), method = 2)
philentropy::jaccard(P =  ZmOs.df.F3$CS1.At[which(ZmOs.df.F$Percentile.Rank.Bin == "First")], Q = ZmOs.df.F3$CS1.At[which(ZmOs.df.F$Percentile.Rank.Bin == "First")], testNA = F)
proxyC::dist(as.matrix(data.frame(Gene1 = ZmOs.df.F$CS1.At[which(ZmOs.df.F$Percentile.Rank.Bin == "First")], Gene2 = ZmOs.df.F$CS1.At[which(ZmOs.df.F$Percentile.Rank.Bin == "First")])), method = "hamming", margin = 2)[1,2]

### with less bins (5 and 3)

ZmOs.df.F$Percentile.Rank.Bin.Five <- NA
ZmOs.df.F$Percentile.Rank.Bin.Five[which(ZmOs.df.F$Percentile.Rank.Score < 0.2)] <- "First" # 660516 regions
ZmOs.df.F$Percentile.Rank.Bin.Five[which(ZmOs.df.F$Percentile.Rank.Score >= 0.2 & ZmOs.df.F$Percentile.Rank.Score < 0.4)] <- "Second" # 660457
ZmOs.df.F$Percentile.Rank.Bin.Five[which(ZmOs.df.F$Percentile.Rank.Score >= 0.4 & ZmOs.df.F$Percentile.Rank.Score < 0.6)] <- "Thirdth" # 660566
ZmOs.df.F$Percentile.Rank.Bin.Five[which(ZmOs.df.F$Percentile.Rank.Score >= 0.6 & ZmOs.df.F$Percentile.Rank.Score < 0.8)] <- "Fourth" # 660867
ZmOs.df.F$Percentile.Rank.Bin.Five[which(ZmOs.df.F$Percentile.Rank.Score >= 0.8 & ZmOs.df.F$Percentile.Rank.Score <= 1)] <- "Fifth" # 659983
(660516 + 660457 + 660566 + 660867 + 659983) == nrow(ZmOs.df.F) # check


CS.Simmilarity.LECIF.ZmOs <- data.frame(Bin = rep(c("First", "Second", "Thirdth", "Fourth", "Fifth"),46),
                                        Metric = c(rep(c("Jaccard"),5),rep(c("SMC"),5),rep(c("Sokal"),5), rep(c("Tanimoto"),5),
                                                   rep(c("Dice"),5), rep(c("Hanman"),5), rep(c("Ochiai"),5), rep(c("Sokal2"),5),
                                                   rep(c("Phi"),5), rep(c("S2"),5), 
                                                   rep(c("Eucledian"),5), rep(c("Sorensen"),5), rep(c("Canberra"),5), rep(c("Wavehedges"),5), 
                                                   rep(c("Tanimoto2"),5), rep(c("Fidelity"),5), rep(c("SquaredChord"),5), rep(c("SquaredChi"),5), 
                                                   rep(c("AdditiveSym"),5), rep(c("Topsoe"),5), rep(c("Kumarjhonson"),5), rep(c("Manhattan"),5),
                                                   rep(c("Lorentzian"),5), rep(c("Czekanowski"),5), rep(c("Hassebrook"),5),  rep(c("SquaredEucledian"),5),
                                                   rep(c("ProbSymm"),5), rep(c("Kullbackleiber"),5), rep(c("Jensenshannon"),5), rep(c("Minkowski"),5), 
                                                   rep(c("Soergel"),5), rep(c("Innerproduct"),5), rep(c("Hellinger"),5), rep(c("Divergence"),5), 
                                                   rep(c("Jeffreys"),5), rep(c("Jensendifference"),5), rep(c("Chebysev"),5), rep(c("Kulczynskid"),5), 
                                                   rep(c("Harmonicmean"),5), rep(c("Matusita"),5), rep(c("Neyman"),5), rep(c("Clark"),5), 
                                                   rep(c("Kdivergence"),5), rep(c("Taneja"),5), 
                                                   rep(c("Maximum"),5), rep(c("Hamming"),5)),
                                        CS1 = NA, CS2 = NA, CS3 = NA, CS4 = NA,
                                        CS5 = NA, CS6 = NA, CS7 = NA, CS8 = NA,
                                        CS9 = NA, CS10 = NA, CS11 = NA, CS12 = NA,
                                        CS13 = NA, CS14 = NA, CS15 = NA, CS16 = NA,
                                        Bivalent = NA, Active = NA, Divergent = NA, Heterochromatin = NA, Quies = NA)

for(i in 1:21){# 16 CS and 5 CS functional groups
  cat(paste0(i / 21 * 100, ' % CS groups started \n'))
  for(j in 1:5){ # lecif score percentile bins
    cat(paste0(j / 5 * 100, ' % LECIF bins started \n'))
    CS.Simmilarity.LECIF.ZmOs[j,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 1)
    CS.Simmilarity.LECIF.ZmOs[j+5,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 2) 
    CS.Simmilarity.LECIF.ZmOs[j+10,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 3)  
    CS.Simmilarity.LECIF.ZmOs[j+15,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 4) 
    CS.Simmilarity.LECIF.ZmOs[j+20,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 5) 
    CS.Simmilarity.LECIF.ZmOs[j+25,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 6) 
    CS.Simmilarity.LECIF.ZmOs[j+30,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 7) 
    CS.Simmilarity.LECIF.ZmOs[j+35,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 8) 
    CS.Simmilarity.LECIF.ZmOs[j+40,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 9) 
    CS.Simmilarity.LECIF.ZmOs[j+45,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 10) 
    CS.Simmilarity.LECIF.ZmOs[j+50,2+i] <- philentropy::euclidean(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F) 
    CS.Simmilarity.LECIF.ZmOs[j+55,2+i] <- philentropy::sorensen(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+60,2+i] <- philentropy::canberra(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+65,2+i] <- philentropy::wave_hedges(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+70,2+i] <- philentropy::tanimoto(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+75,2+i] <- philentropy::fidelity(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+80,2+i] <- philentropy::squared_chord(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+85,2+i] <- philentropy::squared_chi_sq(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+90,2+i] <- philentropy::additive_symm_chi_sq(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+95,2+i] <- philentropy::topsoe(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F,unit = "log2")
    CS.Simmilarity.LECIF.ZmOs[j+100,2+i] <- philentropy::kumar_johnson(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmOs[j+105,2+i] <- philentropy::manhattan(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+110,2+i] <- philentropy::lorentzian(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.ZmOs[j+115,2+i] <- philentropy::czekanowski(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+120,2+i] <- philentropy::kumar_hassebrook(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+125,2+i] <- philentropy::squared_euclidean(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+130,2+i] <- philentropy::prob_symm_chi_sq(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+135,2+i] <- philentropy::kullback_leibler_distance(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmOs[j+140,2+i] <- philentropy::jensen_shannon(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.ZmOs[j+145,2+i] <- philentropy::minkowski(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, n = 3)
    CS.Simmilarity.LECIF.ZmOs[j+150,2+i] <- philentropy::soergel(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+155,2+i] <- philentropy::inner_product(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+160,2+i] <- philentropy::hellinger(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+165,2+i] <- philentropy::divergence_sq(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+170,2+i] <- philentropy::jeffreys(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmOs[j+175,2+i] <- philentropy::jensen_difference(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.ZmOs[j+180,2+i] <- philentropy::chebyshev(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+185,2+i] <- philentropy::kulczynski_d(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmOs[j+190,2+i] <- philentropy::harmonic_mean_dist(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+195,2+i] <- philentropy::matusita(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+200,2+i] <- philentropy::neyman_chi_sq(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmOs[j+205,2+i] <- philentropy::clark_sq(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+210,2+i] <- philentropy::k_divergence(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.ZmOs[j+215,2+i] <- philentropy::taneja(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmOs[j+220,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = "maximum", margin = 2)[1,2]
    CS.Simmilarity.LECIF.ZmOs[j+225,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Five == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = "hamming", margin = 2)[1,2]
    cat(paste0(j / 5 * 100, ' % LECIF bins completed \n'))
  }
  cat(paste0(i / 21 * 100, ' % CS groups completed \n'))
}

write.table(ZmOs.df.F, "H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/ZmOs_CSsimmilarity.DF5.txt", col.names = T, row.names = F, sep = "\t", quote = F)
CS.Simmilarity.LECIF.ZmOs.5 <- CS.Simmilarity.LECIF.ZmOs
write.table(CS.Simmilarity.LECIF.ZmOs, file = "H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/ZmOs_CSsimmilarity.DF.Distances5.txt", col.names = T, row.names = F, sep = "\t", quote = F)
View(CS.Simmilarity.LECIF.ZmOs.5)

ZmOs.df.F$Percentile.Rank.Bin.Three <- NA
ZmOs.df.F$Percentile.Rank.Bin.Three[which(ZmOs.df.F$Percentile.Rank.Score < 0.333333)] <- "First" # 1100826 regions
ZmOs.df.F$Percentile.Rank.Bin.Three[which(ZmOs.df.F$Percentile.Rank.Score >= 0.333333 & ZmOs.df.F$Percentile.Rank.Score < 0.666666)] <- "Second" # 1100799
ZmOs.df.F$Percentile.Rank.Bin.Three[which(ZmOs.df.F$Percentile.Rank.Score >= 0.666666 & ZmOs.df.F$Percentile.Rank.Score <= 1)] <- "Thirdth" # 1100764
( 1100826 + 1100799 + 1100764) == nrow(ZmOs.df.F) # check


CS.Simmilarity.LECIF.ZmOs <- data.frame(Bin = rep(c("First", "Second", "Thirdth"),46),
                                        Metric = c(rep(c("Jaccard"),3),rep(c("SMC"),3),rep(c("Sokal"),3), rep(c("Tanimoto"),3),
                                                   rep(c("Dice"),3), rep(c("Hanman"),3), rep(c("Ochiai"),3), rep(c("Sokal2"),3),
                                                   rep(c("Phi"),3), rep(c("S2"),3), 
                                                   rep(c("Eucledian"),3), rep(c("Sorensen"),3), rep(c("Canberra"),3), rep(c("Wavehedges"),3), 
                                                   rep(c("Tanimoto2"),3), rep(c("Fidelity"),3), rep(c("SquaredChord"),3), rep(c("SquaredChi"),3), 
                                                   rep(c("AdditiveSym"),3), rep(c("Topsoe"),3), rep(c("Kumarjhonson"),3), rep(c("Manhattan"),3),
                                                   rep(c("Lorentzian"),3), rep(c("Czekanowski"),3), rep(c("Hassebrook"),3),  rep(c("SquaredEucledian"),3),
                                                   rep(c("ProbSymm"),3), rep(c("Kullbackleiber"),3), rep(c("Jensenshannon"),3), rep(c("Minkowski"),3), 
                                                   rep(c("Soergel"),3), rep(c("Innerproduct"),3), rep(c("Hellinger"),3), rep(c("Divergence"),3), 
                                                   rep(c("Jeffreys"),3), rep(c("Jensendifference"),3), rep(c("Chebysev"),3), rep(c("Kulczynskid"),3), 
                                                   rep(c("Harmonicmean"),3), rep(c("Matusita"),3), rep(c("Neyman"),3), rep(c("Clark"),3), 
                                                   rep(c("Kdivergence"),3), rep(c("Taneja"),3), 
                                                   rep(c("Maximum"),3), rep(c("Hamming"),3)),
                                        CS1 = NA, CS2 = NA, CS3 = NA, CS4 = NA,
                                        CS5 = NA, CS6 = NA, CS7 = NA, CS8 = NA,
                                        CS9 = NA, CS10 = NA, CS11 = NA, CS12 = NA,
                                        CS13 = NA, CS14 = NA, CS15 = NA, CS16 = NA,
                                        Bivalent = NA, Active = NA, Divergent = NA, Heterochromatin = NA, Quies = NA)

for(i in 1:21){# 16 CS and 5 CS functional groups
  cat(paste0(i / 21 * 100, ' % CS groups started \n'))
  for(j in 1:3){ # lecif score percentile bins
    cat(paste0(j / 3 * 100, ' % LECIF bins started \n'))
    CS.Simmilarity.LECIF.ZmOs[j,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 1)
    CS.Simmilarity.LECIF.ZmOs[j+3,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 2) 
    CS.Simmilarity.LECIF.ZmOs[j+6,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 3)  
    CS.Simmilarity.LECIF.ZmOs[j+9,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 4) 
    CS.Simmilarity.LECIF.ZmOs[j+12,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 5) 
    CS.Simmilarity.LECIF.ZmOs[j+15,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 6) 
    CS.Simmilarity.LECIF.ZmOs[j+18,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 7) 
    CS.Simmilarity.LECIF.ZmOs[j+21,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 8) 
    CS.Simmilarity.LECIF.ZmOs[j+24,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 9) 
    CS.Simmilarity.LECIF.ZmOs[j+27,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = 10) 
    CS.Simmilarity.LECIF.ZmOs[j+30,2+i] <- philentropy::euclidean(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F) 
    CS.Simmilarity.LECIF.ZmOs[j+33,2+i] <- philentropy::sorensen(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+36,2+i] <- philentropy::canberra(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+39,2+i] <- philentropy::wave_hedges(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+42,2+i] <- philentropy::tanimoto(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+45,2+i] <- philentropy::fidelity(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+48,2+i] <- philentropy::squared_chord(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+51,2+i] <- philentropy::squared_chi_sq(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+54,2+i] <- philentropy::additive_symm_chi_sq(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+57,2+i] <- philentropy::topsoe(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F,unit = "log2")
    CS.Simmilarity.LECIF.ZmOs[j+60,2+i] <- philentropy::kumar_johnson(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmOs[j+63,2+i] <- philentropy::manhattan(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+66,2+i] <- philentropy::lorentzian(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.ZmOs[j+69,2+i] <- philentropy::czekanowski(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+72,2+i] <- philentropy::kumar_hassebrook(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+75,2+i] <- philentropy::squared_euclidean(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+78,2+i] <- philentropy::prob_symm_chi_sq(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+81,2+i] <- philentropy::kullback_leibler_distance(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmOs[j+84,2+i] <- philentropy::jensen_shannon(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.ZmOs[j+87,2+i] <- philentropy::minkowski(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, n = 3)
    CS.Simmilarity.LECIF.ZmOs[j+90,2+i] <- philentropy::soergel(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+93,2+i] <- philentropy::inner_product(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+96,2+i] <- philentropy::hellinger(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+99,2+i] <- philentropy::divergence_sq(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+102,2+i] <- philentropy::jeffreys(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmOs[j+105,2+i] <- philentropy::jensen_difference(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.ZmOs[j+108,2+i] <- philentropy::chebyshev(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+111,2+i] <- philentropy::kulczynski_d(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmOs[j+114,2+i] <- philentropy::harmonic_mean_dist(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+117,2+i] <- philentropy::matusita(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+120,2+i] <- philentropy::neyman_chi_sq(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmOs[j+123,2+i] <- philentropy::clark_sq(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F)
    CS.Simmilarity.LECIF.ZmOs[j+126,2+i] <- philentropy::k_divergence(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, unit = "log2")
    CS.Simmilarity.LECIF.ZmOs[j+129,2+i] <- philentropy::taneja(P = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], Q = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i], testNA = F, unit = "log2", epsilon = 0.01)
    CS.Simmilarity.LECIF.ZmOs[j+132,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = "maximum", margin = 2)[1,2]
    CS.Simmilarity.LECIF.ZmOs[j+135,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Zm = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),9+i], CS.Os = ZmOs.df.F[which(ZmOs.df.F$Percentile.Rank.Bin.Three == CS.Simmilarity.LECIF.ZmOs[j,1]),30+i])), method = "hamming", margin = 2)[1,2]
    cat(paste0(j / 3 * 100, ' % LECIF bins completed \n'))
  }
  cat(paste0(i / 21 * 100, ' % CS groups completed \n'))
}

write.table(ZmOs.df.F, "H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/ZmOs_CSsimmilarity.DF3.txt", col.names = T, row.names = F, sep = "\t", quote = F)
CS.Simmilarity.LECIF.ZmOs.3 <- CS.Simmilarity.LECIF.ZmOs
write.table(CS.Simmilarity.LECIF.ZmOs, file = "H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/ZmOs_CSsimmilarity.DF.Distances3.txt", col.names = T, row.names = F, sep = "\t", quote = F)
View(CS.Simmilarity.LECIF.ZmOs.3) ## Dice-like would be the best option


### 1.3 Simmilarity of Chromatin States (hiHMM) over high/low LECIF and high/low PhyloP sets  
## LECIF score high/low and PhyloP score high/low in target (i.e Arabidopsis) and check CS simmilarity on that regions with source (i.e Oryza, Zea) 

###############################################
########## 1.3.1 Arabidopsis - Oryza ######### 
##############################################

## Condition PhyloP scores for all 100% bases in LECIF region

AtOs.LECIFscore.PhyloPsinglent.mean <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/At-Os_LECIFscore_PhyloP_singlent_100overlap_mean.bedgraph", header = F, sep = "\t")
AtOs.df.F3 <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/AtOs_CSsimmilarity.DF3.txt", header = T, sep = "\t")

colnames(AtOs.LECIFscore.PhyloPsinglent.mean) <- c("Chr", "Start", "End", "LECIF", "PhyloP")
AtOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score <- percent_rank(AtOs.LECIFscore.PhyloPsinglent.mean$LECIF)
AtOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score <- percent_rank(AtOs.LECIFscore.PhyloPsinglent.mean$PhyloP)

# After_50;50
AtOs.LECIFscore.PhyloPsinglent.mean$After_5050 <- NA
# low,low 109669 regions
AtOs.LECIFscore.PhyloPsinglent.mean$After_5050[which(AtOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.5 & AtOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.5)] <- "Low,Low"
# high,low 114629 regions
AtOs.LECIFscore.PhyloPsinglent.mean$After_5050[which(AtOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.5 & AtOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.5)] <- "High,Low"
# high,high 109507 regions
AtOs.LECIFscore.PhyloPsinglent.mean$After_5050[which(AtOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.5 & AtOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.5)] <- "High,High"
# low,high 114791 regions
AtOs.LECIFscore.PhyloPsinglent.mean$After_5050[which(AtOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.5 & AtOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.5)] <- "Low,High"
(109669 + 114629 + 109507 + 114791) == nrow(AtOs.LECIFscore.PhyloPsinglent.mean) # check
length(which(is.na(AtOs.LECIFscore.PhyloPsinglent.mean$After_5050) == TRUE)) == 0 # check

# After_60;40
AtOs.LECIFscore.PhyloPsinglent.mean$After_6040 <- NA
# low,low 69742 regions
AtOs.LECIFscore.PhyloPsinglent.mean$After_6040[which(AtOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.4 & AtOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.4)] <- "Low,Low"
# high,low 76138 regions
AtOs.LECIFscore.PhyloPsinglent.mean$After_6040[which(AtOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.6 & AtOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.4)] <- "High,Low"
# high,high 68109 regions
AtOs.LECIFscore.PhyloPsinglent.mean$After_6040[which(AtOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.6 & AtOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.6)] <- "High,High"
# low,high 73305 regions
AtOs.LECIFscore.PhyloPsinglent.mean$After_6040[which(AtOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.4 & AtOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.6)] <- "Low,High"
length(which(is.na(AtOs.LECIFscore.PhyloPsinglent.mean$After_6040) == TRUE)) # check
(69742+76138+68109+73305+161302) == nrow(AtOs.LECIFscore.PhyloPsinglent.mean)

# After_70;30
AtOs.LECIFscore.PhyloPsinglent.mean$After_7030 <- NA
# low,low 39354 regions
AtOs.LECIFscore.PhyloPsinglent.mean$After_7030[which(AtOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.3 & AtOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.3)] <- "Low,Low"
# high,low 43509 regions
AtOs.LECIFscore.PhyloPsinglent.mean$After_7030[which(AtOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.7 & AtOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.3)] <- "High,Low"
# high,high 37379 regions
AtOs.LECIFscore.PhyloPsinglent.mean$After_7030[which(AtOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.7 & AtOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.7)] <- "High,High"
# low,high 39747 regions
AtOs.LECIFscore.PhyloPsinglent.mean$After_7030[which(AtOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.3 & AtOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.7)] <- "Low,High"
length(which(is.na(AtOs.LECIFscore.PhyloPsinglent.mean$After_7030) == TRUE)) # check
(39354+43509+37379+39747+288607) == nrow(AtOs.LECIFscore.PhyloPsinglent.mean)

# After_80;20
AtOs.LECIFscore.PhyloPsinglent.mean$After_8020 <- NA
# low,low 17456 regions
AtOs.LECIFscore.PhyloPsinglent.mean$After_8020[which(AtOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.2 & AtOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.2)] <- "Low,Low"
# high,low 20342 regions
AtOs.LECIFscore.PhyloPsinglent.mean$After_8020[which(AtOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.8 & AtOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.2)] <- "High,Low"
# high,high 16835 regions
AtOs.LECIFscore.PhyloPsinglent.mean$After_8020[which(AtOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.8 & AtOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.8)] <- "High,High"
# low,high 16463 regions
AtOs.LECIFscore.PhyloPsinglent.mean$After_8020[which(AtOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.2 & AtOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.8)] <- "Low,High"
length(which(is.na(AtOs.LECIFscore.PhyloPsinglent.mean$After_8020) == TRUE)) # check
(17456+20342+16835+16463+377500) == nrow(AtOs.LECIFscore.PhyloPsinglent.mean)

# After_90;10
AtOs.LECIFscore.PhyloPsinglent.mean$After_9010 <- NA
# low,low 4565 regions
AtOs.LECIFscore.PhyloPsinglent.mean$After_9010[which(AtOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.1 & AtOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.1)] <- "Low,Low"
# high,low 2985 regions
AtOs.LECIFscore.PhyloPsinglent.mean$After_9010[which(AtOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.9 & AtOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.1)] <- "High,Low"
# high,high 5811 regions
AtOs.LECIFscore.PhyloPsinglent.mean$After_9010[which(AtOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.9 & AtOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.9)] <- "High,High"
# low,high 3721 regions
AtOs.LECIFscore.PhyloPsinglent.mean$After_9010[which(AtOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.1 & AtOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.9)] <- "Low,High"
length(which(is.na(AtOs.LECIFscore.PhyloPsinglent.mean$After_9010) == TRUE)) # check
(4565+2985+5811+3721+431514) == nrow(AtOs.LECIFscore.PhyloPsinglent.mean)

## CS.At and CS.Os columns

AtOs.LECIFscore.PhyloPsinglent.mean$CS1.At <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS2.At <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS3.At <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS4.At <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS5.At <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS6.At <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS7.At <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS8.At <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS9.At <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS10.At <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS11.At <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS12.At <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS13.At <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS14.At <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS15.At <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS16.At <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$Bivalent.At <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$Active.At <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$Divergent.At <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$Heterochromatin.At <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$Quies.At <- NA

AtOs.LECIFscore.PhyloPsinglent.mean$CS1.Os <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS2.Os <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS3.Os <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS4.Os <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS5.Os <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS6.Os <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS7.Os <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS8.Os <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS9.Os <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS10.Os <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS11.Os <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS12.Os <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS13.Os <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS14.Os <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS15.Os <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$CS16.Os <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$Bivalent.Os <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$Active.Os <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$Divergent.Os <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$Heterochromatin.Os <- NA
AtOs.LECIFscore.PhyloPsinglent.mean$Quies.Os <- NA

## loop for annotating my new merged regions with both scores with previously one with CS presence/absence info

for(i in 1:nrow(AtOs.LECIFscore.PhyloPsinglent.mean)){
  print(paste0(i/nrow(AtOs.LECIFscore.PhyloPsinglent.mean)*100, " % covered \n"))
  hit <- which(AtOs.df.F3$Chr.At == AtOs.LECIFscore.PhyloPsinglent.mean$Chr[i] & AtOs.df.F3$Start.At == AtOs.LECIFscore.PhyloPsinglent.mean$Start[i] & AtOs.df.F3$End.At == AtOs.LECIFscore.PhyloPsinglent.mean$End[i])
  if(length(hit) == 1){
    AtOs.LECIFscore.PhyloPsinglent.mean[i,c(13:54)] <- AtOs.df.F3[hit,10:51]
  }else if(length(hit) == 0){
    Chr <- which(AtOs.df.F3$Chr.At == AtOs.LECIFscore.PhyloPsinglent.mean$Chr[i])
    hit2 <- which(AtOs.df.F3$Start.At[Chr] >= AtOs.LECIFscore.PhyloPsinglent.mean$Start[i] & AtOs.df.F3$End.At[Chr] <= AtOs.LECIFscore.PhyloPsinglent.mean$End[i])
    if(length(hit2) == 0){stop("Something go wrong boy ... \n")}
    AtOs.LECIFscore.PhyloPsinglent.mean[i,c(13:54)] <- colSums(AtOs.df.F3[hit2,10:51])
    colto1 <- c(13:54)[which(AtOs.LECIFscore.PhyloPsinglent.mean[i,c(13:54)] > 1)]
    AtOs.LECIFscore.PhyloPsinglent.mean[i, colto1] <- 1
  }else {
    stop("Something go wrong boy ... \n")
  }
}

write.table(AtOs.LECIFscore.PhyloPsinglent.mean, "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/AtOs_CSsimmilarity.DF.LECIFPhyloP.AfterCutOff.txt", col.names = T, row.names = F, sep = "\t", quote = F)

## loop for computing simmilarity/distance between different subsets of LECIF and PhyloP score with different cutoffs (After ...)

CS.Simmilarity.LECIF.PhyloP.AtOs.After <- data.frame(Bin = rep(c("High,High", "Low,Low", "High,Low", "Low,High"),46),
                                        Metric = c(rep(c("Jaccard"),4),rep(c("SMC"),4),rep(c("Sokal"),4), rep(c("Tanimoto"),4),
                                                   rep(c("Dice"),4), rep(c("Hanman"),4), rep(c("Ochiai"),4), rep(c("Sokal2"),4),
                                                   rep(c("Phi"),4), rep(c("S2"),4), 
                                                   rep(c("Eucledian"),4), rep(c("Sorensen"),4), rep(c("Canberra"),4), rep(c("Wavehedges"),4), 
                                                   rep(c("Tanimoto2"),4), rep(c("Fidelity"),4), rep(c("SquaredChord"),4), rep(c("SquaredChi"),4), 
                                                   rep(c("AdditiveSym"),4), rep(c("Topsoe"),4), rep(c("Kumarjhonson"),4), rep(c("Manhattan"),4),
                                                   rep(c("Lorentzian"),4), rep(c("Czekanowski"),4), rep(c("Hassebrook"),4),  rep(c("SquaredEucledian"),4),
                                                   rep(c("ProbSymm"),4), rep(c("Kullbackleiber"),4), rep(c("Jensenshannon"),4), rep(c("Minkowski"),4), 
                                                   rep(c("Soergel"),4), rep(c("Innerproduct"),4), rep(c("Hellinger"),4), rep(c("Divergence"),4), 
                                                   rep(c("Jeffreys"),4), rep(c("Jensendifference"),4), rep(c("Chebysev"),4), rep(c("Kulczynskid"),4), 
                                                   rep(c("Harmonicmean"),4), rep(c("Matusita"),4), rep(c("Neyman"),4), rep(c("Clark"),4), 
                                                   rep(c("Kdivergence"),4), rep(c("Taneja"),4), 
                                                   rep(c("Maximum"),4), rep(c("Hamming"),4)),
                                        CS1 = NA, CS2 = NA, CS3 = NA, CS4 = NA,
                                        CS5 = NA, CS6 = NA, CS7 = NA, CS8 = NA,
                                        CS9 = NA, CS10 = NA, CS11 = NA, CS12 = NA,
                                        CS13 = NA, CS14 = NA, CS15 = NA, CS16 = NA,
                                        Bivalent = NA, Active = NA, Divergent = NA, Heterochromatin = NA, Quies = NA)

List.CS.Simmilarity.LECIF.PhyloP.AtOs <- list(After_5050 = CS.Simmilarity.LECIF.PhyloP.AtOs.After,
                                              After_6040 = CS.Simmilarity.LECIF.PhyloP.AtOs.After,
                                              After_7030 = CS.Simmilarity.LECIF.PhyloP.AtOs.After,
                                              After_8020 = CS.Simmilarity.LECIF.PhyloP.AtOs.After,
                                              After_9010 = CS.Simmilarity.LECIF.PhyloP.AtOs.After)
bincolumns <- c(8:12)
for(z in 1:length(bincolumns)){ # 4 different cutoffs
  cat(paste0(colnames(AtOs.LECIFscore.PhyloPsinglent.mean)[bincolumns[z]], " cutoff processing ... \n"))
  for(i in 1:21){# 16 CS and 5 CS functional groups
    cat(paste0(i / 21 * 100, ' % CS groups started \n'))
    for(j in 1:4){ # lecif score percentile bins
      cat(paste0(j / 4 * 100, ' % cutoffs bins started \n'))
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], CS.Os = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i])), method = 1)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+4,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], CS.Os = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i])), method = 2) 
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+8,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], CS.Os = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i])), method = 3)  
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+12,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], CS.Os = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i])), method = 4) 
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+16,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], CS.Os = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i])), method = 5) 
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+20,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], CS.Os = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i])), method = 6) 
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+24,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], CS.Os = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i])), method = 7) 
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+28,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], CS.Os = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i])), method = 8) 
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+32,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], CS.Os = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i])), method = 9) 
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+36,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], CS.Os = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i])), method = 10) 
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+40,2+i] <- philentropy::euclidean(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F) 
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+44,2+i] <- philentropy::sorensen(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+48,2+i] <- philentropy::canberra(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+52,2+i] <- philentropy::wave_hedges(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+56,2+i] <- philentropy::tanimoto(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+60,2+i] <- philentropy::fidelity(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+64,2+i] <- philentropy::squared_chord(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+68,2+i] <- philentropy::squared_chi_sq(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+72,2+i] <- philentropy::additive_symm_chi_sq(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+76,2+i] <- philentropy::topsoe(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F,unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+80,2+i] <- philentropy::kumar_johnson(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F, epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+84,2+i] <- philentropy::manhattan(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+88,2+i] <- philentropy::lorentzian(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F, unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+92,2+i] <- philentropy::czekanowski(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+96,2+i] <- philentropy::kumar_hassebrook(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+100,2+i] <- philentropy::squared_euclidean(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+104,2+i] <- philentropy::prob_symm_chi_sq(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+108,2+i] <- philentropy::kullback_leibler_distance(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F, unit = "log2", epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+112,2+i] <- philentropy::jensen_shannon(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F, unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+116,2+i] <- philentropy::minkowski(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F, n = 3)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+120,2+i] <- philentropy::soergel(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+124,2+i] <- philentropy::inner_product(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+128,2+i] <- philentropy::hellinger(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+132,2+i] <- philentropy::divergence_sq(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+136,2+i] <- philentropy::jeffreys(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F, unit = "log2", epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+140,2+i] <- philentropy::jensen_difference(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F, unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+144,2+i] <- philentropy::chebyshev(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+148,2+i] <- philentropy::kulczynski_d(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F, epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+152,2+i] <- philentropy::harmonic_mean_dist(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+156,2+i] <- philentropy::matusita(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+160,2+i] <- philentropy::neyman_chi_sq(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F, epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+164,2+i] <- philentropy::clark_sq(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+168,2+i] <- philentropy::k_divergence(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F, unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+172,2+i] <- philentropy::taneja(P = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], Q = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i], testNA = F, unit = "log2", epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+176,2+i] <- proxyC::dist(as.matrix(data.frame(CS.At = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], CS.Os = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i])), method = "maximum", margin = 2)[1,2]
      List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j+180,2+i] <- proxyC::dist(as.matrix(data.frame(CS.At = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),12+i], CS.Os = AtOs.LECIFscore.PhyloPsinglent.mean[which(AtOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtOs[[z]][j,1]),33+i])), method = "hamming", margin = 2)[1,2]
      cat(paste0(j / 4 * 100, " % cutoff bins completed \n"))
    }
    cat(paste0(i / 21 * 100, " % CS groups completed \n"))
  }
  cat(paste0(colnames(AtOs.LECIFscore.PhyloPsinglent.mean)[bincolumns[z]], " cutoff computation completed \n"))
}

write.table(List.CS.Simmilarity.LECIF.PhyloP.AtOs$After_5050, file = "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/AtOs_CSsimmilarity.DF.Distances.After5050.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(List.CS.Simmilarity.LECIF.PhyloP.AtOs$After_6040, file = "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/AtOs_CSsimmilarity.DF.Distances.After6040.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(List.CS.Simmilarity.LECIF.PhyloP.AtOs$After_7030, file = "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/AtOs_CSsimmilarity.DF.Distances.After7030.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(List.CS.Simmilarity.LECIF.PhyloP.AtOs$After_8020, file = "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/AtOs_CSsimmilarity.DF.Distances.After8020.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(List.CS.Simmilarity.LECIF.PhyloP.AtOs$After_9010, file = "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/AtOs_CSsimmilarity.DF.Distances.After9010.txt", col.names = T, row.names = F, sep = "\t", quote = F)

###############################################
########## 1.3.2 Arabidopsis - Zea ########### 
##############################################

AtZm.LECIFscore.PhyloPsinglent.mean <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/At-Zm_LECIFscore_PhyloP_singlent_100overlap_mean.bedgraph", header = F, sep = "\t")
AtZm.df.F3 <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/AtZm_CSsimmilarity.DF3.txt", header = T, sep = "\t")

colnames(AtZm.LECIFscore.PhyloPsinglent.mean) <- c("Chr", "Start", "End", "LECIF", "PhyloP")
AtZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score <- percent_rank(AtZm.LECIFscore.PhyloPsinglent.mean$LECIF)
AtZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score <- percent_rank(AtZm.LECIFscore.PhyloPsinglent.mean$PhyloP)

# After_50;50
AtZm.LECIFscore.PhyloPsinglent.mean$After_5050 <- NA
# low,low 88184 regions
AtZm.LECIFscore.PhyloPsinglent.mean$After_5050[which(AtZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.5 & AtZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.5)] <- "Low,Low"
# high,low 85088 regions
AtZm.LECIFscore.PhyloPsinglent.mean$After_5050[which(AtZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.5 & AtZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.5)] <- "High,Low"
# high,high 88180 regions
AtZm.LECIFscore.PhyloPsinglent.mean$After_5050[which(AtZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.5 & AtZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.5)] <- "High,High"
# low,high 85092 regions
AtZm.LECIFscore.PhyloPsinglent.mean$After_5050[which(AtZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.5 & AtZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.5)] <- "Low,High"
(88184 + 85088 + 88180 + 85092) == nrow(AtZm.LECIFscore.PhyloPsinglent.mean) # check
length(which(is.na(AtZm.LECIFscore.PhyloPsinglent.mean$After_5050) == TRUE)) == 0 # check

# After_60;40
AtZm.LECIFscore.PhyloPsinglent.mean$After_6040 <- NA
# low,low 58097 regions
AtZm.LECIFscore.PhyloPsinglent.mean$After_6040[which(AtZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.4 & AtZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.4)] <- "Low,Low"
# high,low 54146 regions
AtZm.LECIFscore.PhyloPsinglent.mean$After_6040[which(AtZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.6 & AtZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.4)] <- "High,Low"
# high,high 56295 regions
AtZm.LECIFscore.PhyloPsinglent.mean$After_6040[which(AtZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.6 & AtZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.6)] <- "High,High"
# low,high 53443 regions
AtZm.LECIFscore.PhyloPsinglent.mean$After_6040[which(AtZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.4 & AtZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.6)] <- "Low,High"
length(which(is.na(AtZm.LECIFscore.PhyloPsinglent.mean$After_6040) == TRUE)) # check
(58097+54146+56295+53443+124563) == nrow(AtZm.LECIFscore.PhyloPsinglent.mean)

# After_70;30
AtZm.LECIFscore.PhyloPsinglent.mean$After_7030 <- NA
# low,low 33984 regions
AtZm.LECIFscore.PhyloPsinglent.mean$After_7030[which(AtZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.3 & AtZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.3)] <- "Low,Low"
# high,low 30811 regions
AtZm.LECIFscore.PhyloPsinglent.mean$After_7030[which(AtZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.7 & AtZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.3)] <- "High,Low"
# high,high 31933 regions
AtZm.LECIFscore.PhyloPsinglent.mean$After_7030[which(AtZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.7 & AtZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.7)] <- "High,High"
# low,high 28806 regions
AtZm.LECIFscore.PhyloPsinglent.mean$After_7030[which(AtZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.3 & AtZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.7)] <- "Low,High"
length(which(is.na(AtZm.LECIFscore.PhyloPsinglent.mean$After_7030) == TRUE)) # check
(33984+30811+31933+28806+221010) == nrow(AtZm.LECIFscore.PhyloPsinglent.mean)

# After_80;20
AtZm.LECIFscore.PhyloPsinglent.mean$After_8020 <- NA
# low,low 16010 regions
AtZm.LECIFscore.PhyloPsinglent.mean$After_8020[which(AtZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.2 & AtZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.2)] <- "Low,Low"
# high,low 15530 regions
AtZm.LECIFscore.PhyloPsinglent.mean$After_8020[which(AtZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.8 & AtZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.2)] <- "High,Low"
# high,high 14106 regions
AtZm.LECIFscore.PhyloPsinglent.mean$After_8020[which(AtZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.8 & AtZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.8)] <- "High,High"
# low,high 12115 regions
AtZm.LECIFscore.PhyloPsinglent.mean$After_8020[which(AtZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.2 & AtZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.8)] <- "Low,High"
length(which(is.na(AtZm.LECIFscore.PhyloPsinglent.mean$After_8020) == TRUE)) # check
(16010+15530+14106+12115+288783) == nrow(AtZm.LECIFscore.PhyloPsinglent.mean)

# After_90;10
AtZm.LECIFscore.PhyloPsinglent.mean$After_9010 <- NA
# low,low 4218 regions
AtZm.LECIFscore.PhyloPsinglent.mean$After_9010[which(AtZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.1 & AtZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.1)] <- "Low,Low"
# high,low 3834 regions
AtZm.LECIFscore.PhyloPsinglent.mean$After_9010[which(AtZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.9 & AtZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.1)] <- "High,Low"
# high,high 3803 regions
AtZm.LECIFscore.PhyloPsinglent.mean$After_9010[which(AtZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.9 & AtZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.9)] <- "High,High"
# low,high 2887 regions
AtZm.LECIFscore.PhyloPsinglent.mean$After_9010[which(AtZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.1 & AtZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.9)] <- "Low,High"
length(which(is.na(AtZm.LECIFscore.PhyloPsinglent.mean$After_9010) == TRUE)) # check
(4218+3834+3803+2887+331802) == nrow(AtZm.LECIFscore.PhyloPsinglent.mean)

## CS.At and CS.Zm columns

AtZm.LECIFscore.PhyloPsinglent.mean$CS1.At <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS2.At <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS3.At <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS4.At <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS5.At <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS6.At <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS7.At <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS8.At <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS9.At <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS10.At <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS11.At <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS12.At <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS13.At <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS14.At <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS15.At <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS16.At <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$Bivalent.At <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$Active.At <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$Divergent.At <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$Heterochromatin.At <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$Quies.At <- NA

AtZm.LECIFscore.PhyloPsinglent.mean$CS1.Zm <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS2.Zm <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS3.Zm <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS4.Zm <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS5.Zm <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS6.Zm <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS7.Zm <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS8.Zm <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS9.Zm <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS10.Zm <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS11.Zm <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS12.Zm <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS13.Zm <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS14.Zm <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS15.Zm <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$CS16.Zm <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$Bivalent.Zm <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$Active.Zm <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$Divergent.Zm <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$Heterochromatin.Zm <- NA
AtZm.LECIFscore.PhyloPsinglent.mean$Quies.Zm <- NA

## loop for annotating my new merged regions with both scores with previously one with CS presence/absence info

for(i in 1:nrow(AtZm.LECIFscore.PhyloPsinglent.mean)){
  print(paste0(i/nrow(AtZm.LECIFscore.PhyloPsinglent.mean)*100, " % covered \n"))
  hit <- which(AtZm.df.F3$Chr.At == AtZm.LECIFscore.PhyloPsinglent.mean$Chr[i] & AtZm.df.F3$Start.At == AtZm.LECIFscore.PhyloPsinglent.mean$Start[i] & AtZm.df.F3$End.At == AtZm.LECIFscore.PhyloPsinglent.mean$End[i])
  if(length(hit) == 1){
    AtZm.LECIFscore.PhyloPsinglent.mean[i,c(13:54)] <- AtZm.df.F3[hit,10:51]
  }else if(length(hit) == 0){
    Chr <- which(AtZm.df.F3$Chr.At == AtZm.LECIFscore.PhyloPsinglent.mean$Chr[i])
    hit2 <- which(AtZm.df.F3$Start.At[Chr] >= AtZm.LECIFscore.PhyloPsinglent.mean$Start[i] & AtZm.df.F3$End.At[Chr] <= AtZm.LECIFscore.PhyloPsinglent.mean$End[i])
    if(length(hit2) == 0){stop("Something go wrong boy ... \n")}
    AtZm.LECIFscore.PhyloPsinglent.mean[i,c(13:54)] <- colSums(AtZm.df.F3[hit2,10:51])
    colto1 <- c(13:54)[which(AtZm.LECIFscore.PhyloPsinglent.mean[i,c(13:54)] > 1)]
    AtZm.LECIFscore.PhyloPsinglent.mean[i, colto1] <- 1
  }else {
    stop("Something go wrong boy ... \n")
  }
}

write.table(AtZm.LECIFscore.PhyloPsinglent.mean, "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/AtZm_CSsimmilarity.DF.LECIFPhyloP.AfterCutOff.txt", col.names = T, row.names = F, sep = "\t", quote = F)

## loop for computing simmilarity/distance between different subsets of LECIF and PhyloP score with different cutoffs (After ...)

CS.Simmilarity.LECIF.PhyloP.AtZm.After <- data.frame(Bin = rep(c("High,High", "Low,Low", "High,Low", "Low,High"),46),
                                                     Metric = c(rep(c("Jaccard"),4),rep(c("SMC"),4),rep(c("Sokal"),4), rep(c("Tanimoto"),4),
                                                                rep(c("Dice"),4), rep(c("Hanman"),4), rep(c("Ochiai"),4), rep(c("Sokal2"),4),
                                                                rep(c("Phi"),4), rep(c("S2"),4), 
                                                                rep(c("Eucledian"),4), rep(c("Sorensen"),4), rep(c("Canberra"),4), rep(c("Wavehedges"),4), 
                                                                rep(c("Tanimoto2"),4), rep(c("Fidelity"),4), rep(c("SquaredChord"),4), rep(c("SquaredChi"),4), 
                                                                rep(c("AdditiveSym"),4), rep(c("Topsoe"),4), rep(c("Kumarjhonson"),4), rep(c("Manhattan"),4),
                                                                rep(c("Lorentzian"),4), rep(c("Czekanowski"),4), rep(c("Hassebrook"),4),  rep(c("SquaredEucledian"),4),
                                                                rep(c("ProbSymm"),4), rep(c("Kullbackleiber"),4), rep(c("Jensenshannon"),4), rep(c("Minkowski"),4), 
                                                                rep(c("Soergel"),4), rep(c("Innerproduct"),4), rep(c("Hellinger"),4), rep(c("Divergence"),4), 
                                                                rep(c("Jeffreys"),4), rep(c("Jensendifference"),4), rep(c("Chebysev"),4), rep(c("Kulczynskid"),4), 
                                                                rep(c("Harmonicmean"),4), rep(c("Matusita"),4), rep(c("Neyman"),4), rep(c("Clark"),4), 
                                                                rep(c("Kdivergence"),4), rep(c("Taneja"),4), 
                                                                rep(c("Maximum"),4), rep(c("Hamming"),4)),
                                                     CS1 = NA, CS2 = NA, CS3 = NA, CS4 = NA,
                                                     CS5 = NA, CS6 = NA, CS7 = NA, CS8 = NA,
                                                     CS9 = NA, CS10 = NA, CS11 = NA, CS12 = NA,
                                                     CS13 = NA, CS14 = NA, CS15 = NA, CS16 = NA,
                                                     Bivalent = NA, Active = NA, Divergent = NA, Heterochromatin = NA, Quies = NA)

List.CS.Simmilarity.LECIF.PhyloP.AtZm <- list(After_5050 = CS.Simmilarity.LECIF.PhyloP.AtZm.After,
                                              After_6040 = CS.Simmilarity.LECIF.PhyloP.AtZm.After,
                                              After_7030 = CS.Simmilarity.LECIF.PhyloP.AtZm.After,
                                              After_8020 = CS.Simmilarity.LECIF.PhyloP.AtZm.After,
                                              After_9010 = CS.Simmilarity.LECIF.PhyloP.AtZm.After)
bincolumns <- c(8:12)
for(z in 1:length(bincolumns)){ # 4 different cutoffs
  cat(paste0(colnames(AtZm.LECIFscore.PhyloPsinglent.mean)[bincolumns[z]], " cutoff processing ... \n"))
  for(i in 1:21){# 16 CS and 5 CS functional groups
    cat(paste0(i / 21 * 100, ' % CS groups started \n'))
    for(j in 1:4){ # lecif score percentile bins
      cat(paste0(j / 4 * 100, ' % cutoffs bins started \n'))
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], CS.Zm = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i])), method = 1)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+4,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], CS.Zm = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i])), method = 2) 
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+8,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], CS.Zm = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i])), method = 3)  
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+12,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], CS.Zm = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i])), method = 4) 
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+16,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], CS.Zm = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i])), method = 5) 
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+20,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], CS.Zm = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i])), method = 6) 
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+24,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], CS.Zm = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i])), method = 7) 
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+28,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], CS.Zm = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i])), method = 8) 
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+32,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], CS.Zm = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i])), method = 9) 
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+36,2+i] <- ade4::dist.binary(t(data.frame(CS.At = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], CS.Zm = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i])), method = 10) 
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+40,2+i] <- philentropy::euclidean(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F) 
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+44,2+i] <- philentropy::sorensen(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+48,2+i] <- philentropy::canberra(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+52,2+i] <- philentropy::wave_hedges(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+56,2+i] <- philentropy::tanimoto(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+60,2+i] <- philentropy::fidelity(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+64,2+i] <- philentropy::squared_chord(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+68,2+i] <- philentropy::squared_chi_sq(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+72,2+i] <- philentropy::additive_symm_chi_sq(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+76,2+i] <- philentropy::topsoe(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F,unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+80,2+i] <- philentropy::kumar_johnson(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F, epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+84,2+i] <- philentropy::manhattan(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+88,2+i] <- philentropy::lorentzian(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F, unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+92,2+i] <- philentropy::czekanowski(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+96,2+i] <- philentropy::kumar_hassebrook(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+100,2+i] <- philentropy::squared_euclidean(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+104,2+i] <- philentropy::prob_symm_chi_sq(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+108,2+i] <- philentropy::kullback_leibler_distance(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F, unit = "log2", epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+112,2+i] <- philentropy::jensen_shannon(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F, unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+116,2+i] <- philentropy::minkowski(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F, n = 3)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+120,2+i] <- philentropy::soergel(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+124,2+i] <- philentropy::inner_product(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+128,2+i] <- philentropy::hellinger(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+132,2+i] <- philentropy::divergence_sq(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+136,2+i] <- philentropy::jeffreys(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F, unit = "log2", epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+140,2+i] <- philentropy::jensen_difference(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F, unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+144,2+i] <- philentropy::chebyshev(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+148,2+i] <- philentropy::kulczynski_d(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F, epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+152,2+i] <- philentropy::harmonic_mean_dist(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+156,2+i] <- philentropy::matusita(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+160,2+i] <- philentropy::neyman_chi_sq(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F, epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+164,2+i] <- philentropy::clark_sq(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+168,2+i] <- philentropy::k_divergence(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F, unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+172,2+i] <- philentropy::taneja(P = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], Q = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i], testNA = F, unit = "log2", epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+176,2+i] <- proxyC::dist(as.matrix(data.frame(CS.At = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], CS.Zm = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i])), method = "maximum", margin = 2)[1,2]
      List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j+180,2+i] <- proxyC::dist(as.matrix(data.frame(CS.At = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),12+i], CS.Zm = AtZm.LECIFscore.PhyloPsinglent.mean[which(AtZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.AtZm[[z]][j,1]),33+i])), method = "hamming", margin = 2)[1,2]
      cat(paste0(j / 4 * 100, " % cutoff bins completed \n"))
    }
    cat(paste0(i / 21 * 100, " % CS groups completed \n"))
  }
  cat(paste0(colnames(AtZm.LECIFscore.PhyloPsinglent.mean)[bincolumns[z]], " cutoff computation completed \n"))
}

write.table(List.CS.Simmilarity.LECIF.PhyloP.AtZm$After_5050, file = "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/AtZm_CSsimmilarity.DF.Distances.After5050.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(List.CS.Simmilarity.LECIF.PhyloP.AtZm$After_6040, file = "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/AtZm_CSsimmilarity.DF.Distances.After6040.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(List.CS.Simmilarity.LECIF.PhyloP.AtZm$After_7030, file = "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/AtZm_CSsimmilarity.DF.Distances.After7030.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(List.CS.Simmilarity.LECIF.PhyloP.AtZm$After_8020, file = "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/AtZm_CSsimmilarity.DF.Distances.After8020.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(List.CS.Simmilarity.LECIF.PhyloP.AtZm$After_9010, file = "H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/AtZm_CSsimmilarity.DF.Distances.After9010.txt", col.names = T, row.names = F, sep = "\t", quote = F)

###############################################
########## 1.3.3 Oryza - Arabidopsis ######### 
##############################################

## Condition PhyloP scores for all 100% bases in LECIF region

OsAt.LECIFscore.PhyloPsinglent.mean <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/Os-At_LECIFscore_PhyloP_singlent_100overlap_mean.bedgraph", header = F, sep = "\t")
OsAt.df.F3 <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/OsAt_CSsimmilarity.DF3.txt", header = T, sep = "\t")

colnames(OsAt.LECIFscore.PhyloPsinglent.mean) <- c("Chr", "Start", "End", "LECIF", "PhyloP")
OsAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score <- percent_rank(OsAt.LECIFscore.PhyloPsinglent.mean$LECIF)
OsAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score <- percent_rank(OsAt.LECIFscore.PhyloPsinglent.mean$PhyloP)

# After_50;50
OsAt.LECIFscore.PhyloPsinglent.mean$After_5050 <- NA
# low,low 149450 regions
OsAt.LECIFscore.PhyloPsinglent.mean$After_5050[which(OsAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.5 & OsAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.5)] <- "Low,Low"
# high,low 149761 regions
OsAt.LECIFscore.PhyloPsinglent.mean$After_5050[which(OsAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.5 & OsAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.5)] <- "High,Low"
# high,high 149182 regions
OsAt.LECIFscore.PhyloPsinglent.mean$After_5050[which(OsAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.5 & OsAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.5)] <- "High,High"
# low,high 150028 regions
OsAt.LECIFscore.PhyloPsinglent.mean$After_5050[which(OsAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.5 & OsAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.5)] <- "Low,High"
(149450 + 149761 + 149182 + 150028) == nrow(OsAt.LECIFscore.PhyloPsinglent.mean) # check
length(which(is.na(OsAt.LECIFscore.PhyloPsinglent.mean$After_5050) == TRUE)) == 0 # check

# After_60;40
OsAt.LECIFscore.PhyloPsinglent.mean$After_6040 <- NA
# low,low 94410 regions
OsAt.LECIFscore.PhyloPsinglent.mean$After_6040[which(OsAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.4 & OsAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.4)] <- "Low,Low"
# high,low 94429 regions
OsAt.LECIFscore.PhyloPsinglent.mean$After_6040[which(OsAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.6 & OsAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.4)] <- "High,Low"
# high,high 97735 regions
OsAt.LECIFscore.PhyloPsinglent.mean$After_6040[which(OsAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.6 & OsAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.6)] <- "High,High"
# low,high 94922 regions
OsAt.LECIFscore.PhyloPsinglent.mean$After_6040[which(OsAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.4 & OsAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.6)] <- "Low,High"
length(which(is.na(OsAt.LECIFscore.PhyloPsinglent.mean$After_6040) == TRUE)) # check
(94410+94429+97735+94922+216925) == nrow(OsAt.LECIFscore.PhyloPsinglent.mean)

# After_70;30
OsAt.LECIFscore.PhyloPsinglent.mean$After_7030 <- NA
# low,low 52189 regions
OsAt.LECIFscore.PhyloPsinglent.mean$After_7030[which(OsAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.3 & OsAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.3)] <- "Low,Low"
# high,low 49486 regions
OsAt.LECIFscore.PhyloPsinglent.mean$After_7030[which(OsAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.7 & OsAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.3)] <- "High,Low"
# high,high 57600 regions
OsAt.LECIFscore.PhyloPsinglent.mean$After_7030[which(OsAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.7 & OsAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.7)] <- "High,High"
# low,high 51445 regions
OsAt.LECIFscore.PhyloPsinglent.mean$After_7030[which(OsAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.3 & OsAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.7)] <- "Low,High"
length(which(is.na(OsAt.LECIFscore.PhyloPsinglent.mean$After_7030) == TRUE)) # check
(52189+49486+57600+51445+387701) == nrow(OsAt.LECIFscore.PhyloPsinglent.mean)

# After_80;20
OsAt.LECIFscore.PhyloPsinglent.mean$After_8020 <- NA
# low,low 23824 regions
OsAt.LECIFscore.PhyloPsinglent.mean$After_8020[which(OsAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.2 & OsAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.2)] <- "Low,Low"
# high,low 23440 regions
OsAt.LECIFscore.PhyloPsinglent.mean$After_8020[which(OsAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.8 & OsAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.2)] <- "High,Low"
# high,high 26515 regions
OsAt.LECIFscore.PhyloPsinglent.mean$After_8020[which(OsAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.8 & OsAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.8)] <- "High,High"
# low,high 21140 regions
OsAt.LECIFscore.PhyloPsinglent.mean$After_8020[which(OsAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.2 & OsAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.8)] <- "Low,High"
length(which(is.na(OsAt.LECIFscore.PhyloPsinglent.mean$After_8020) == TRUE)) # check
(23824+23440+26515+21140+503502) == nrow(OsAt.LECIFscore.PhyloPsinglent.mean)

# After_90;10
OsAt.LECIFscore.PhyloPsinglent.mean$After_9010 <- NA
# low,low 6570 regions
OsAt.LECIFscore.PhyloPsinglent.mean$After_9010[which(OsAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.1 & OsAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.1)] <- "Low,Low"
# high,low 5663 regions
OsAt.LECIFscore.PhyloPsinglent.mean$After_9010[which(OsAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.9 & OsAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.1)] <- "High,Low"
# high,high 7778 regions
OsAt.LECIFscore.PhyloPsinglent.mean$After_9010[which(OsAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.9 & OsAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.9)] <- "High,High"
# low,high 4421 regions
OsAt.LECIFscore.PhyloPsinglent.mean$After_9010[which(OsAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.1 & OsAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.9)] <- "Low,High"
length(which(is.na(OsAt.LECIFscore.PhyloPsinglent.mean$After_9010) == TRUE)) # check
(6570+5663+7778+4421+573989) == nrow(OsAt.LECIFscore.PhyloPsinglent.mean)

## CS.At and CS.Os columns

OsAt.LECIFscore.PhyloPsinglent.mean$CS1.Os <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS2.Os <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS3.Os <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS4.Os <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS5.Os <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS6.Os <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS7.Os <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS8.Os <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS9.Os <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS10.Os <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS11.Os <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS12.Os <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS13.Os <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS14.Os <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS15.Os <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS16.Os <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$Bivalent.Os <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$Active.Os <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$Divergent.Os <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$Heterochromatin.Os <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$Quies.Os <- NA

OsAt.LECIFscore.PhyloPsinglent.mean$CS1.At <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS2.At <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS3.At <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS4.At <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS5.At <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS6.At <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS7.At <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS8.At <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS9.At <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS10.At <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS11.At <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS12.At <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS13.At <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS14.At <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS15.At <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$CS16.At <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$Bivalent.At <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$Active.At <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$Divergent.At <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$Heterochromatin.At <- NA
OsAt.LECIFscore.PhyloPsinglent.mean$Quies.At <- NA

## loop for annotating my new merged regions with both scores with previously one with CS presence/absence info

for(i in 1:nrow(OsAt.LECIFscore.PhyloPsinglent.mean)){
  print(paste0(i/nrow(OsAt.LECIFscore.PhyloPsinglent.mean)*100, " % covered \n"))
  hit <- which(OsAt.df.F3$Chr.Os == OsAt.LECIFscore.PhyloPsinglent.mean$Chr[i] & OsAt.df.F3$Start.Os == OsAt.LECIFscore.PhyloPsinglent.mean$Start[i] & OsAt.df.F3$End.Os == OsAt.LECIFscore.PhyloPsinglent.mean$End[i])
  if(length(hit) == 1){
    OsAt.LECIFscore.PhyloPsinglent.mean[i,c(13:54)] <- OsAt.df.F3[hit,10:51]
  }else if(length(hit) == 0){
    Chr <- which(OsAt.df.F3$Chr.Os == OsAt.LECIFscore.PhyloPsinglent.mean$Chr[i])
    hit2 <- which(OsAt.df.F3$Start.Os[Chr] >= OsAt.LECIFscore.PhyloPsinglent.mean$Start[i] & OsAt.df.F3$End.Os[Chr] <= OsAt.LECIFscore.PhyloPsinglent.mean$End[i])
    if(length(hit2) == 0){stop("Something go wrong boy ... \n")}
    OsAt.LECIFscore.PhyloPsinglent.mean[i,c(13:54)] <- colSums(OsAt.df.F3[hit2,10:51])
    colto1 <- c(13:54)[which(OsAt.LECIFscore.PhyloPsinglent.mean[i,c(13:54)] > 1)]
    OsAt.LECIFscore.PhyloPsinglent.mean[i, colto1] <- 1
  }else {
    stop("Something go wrong boy ... \n")
  }
}

write.table(OsAt.LECIFscore.PhyloPsinglent.mean, "H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/OsAt_CSsimmilarity.DF.LECIFPhyloP.AfterCutOff.txt", col.names = T, row.names = F, sep = "\t", quote = F)

## loop for computing simmilarity/distance between different subsets of LECIF and PhyloP score with different cutoffs (After ...)

CS.Simmilarity.LECIF.PhyloP.OsAt.After <- data.frame(Bin = rep(c("High,High", "Low,Low", "High,Low", "Low,High"),46),
                                                     Metric = c(rep(c("Jaccard"),4),rep(c("SMC"),4),rep(c("Sokal"),4), rep(c("Tanimoto"),4),
                                                                rep(c("Dice"),4), rep(c("Hanman"),4), rep(c("Ochiai"),4), rep(c("Sokal2"),4),
                                                                rep(c("Phi"),4), rep(c("S2"),4), 
                                                                rep(c("Eucledian"),4), rep(c("Sorensen"),4), rep(c("Canberra"),4), rep(c("Wavehedges"),4), 
                                                                rep(c("Tanimoto2"),4), rep(c("Fidelity"),4), rep(c("SquaredChord"),4), rep(c("SquaredChi"),4), 
                                                                rep(c("AdditiveSym"),4), rep(c("Topsoe"),4), rep(c("Kumarjhonson"),4), rep(c("Manhattan"),4),
                                                                rep(c("Lorentzian"),4), rep(c("Czekanowski"),4), rep(c("Hassebrook"),4),  rep(c("SquaredEucledian"),4),
                                                                rep(c("ProbSymm"),4), rep(c("Kullbackleiber"),4), rep(c("Jensenshannon"),4), rep(c("Minkowski"),4), 
                                                                rep(c("Soergel"),4), rep(c("Innerproduct"),4), rep(c("Hellinger"),4), rep(c("Divergence"),4), 
                                                                rep(c("Jeffreys"),4), rep(c("Jensendifference"),4), rep(c("Chebysev"),4), rep(c("Kulczynskid"),4), 
                                                                rep(c("Harmonicmean"),4), rep(c("Matusita"),4), rep(c("Neyman"),4), rep(c("Clark"),4), 
                                                                rep(c("Kdivergence"),4), rep(c("Taneja"),4), 
                                                                rep(c("Maximum"),4), rep(c("Hamming"),4)),
                                                     CS1 = NA, CS2 = NA, CS3 = NA, CS4 = NA,
                                                     CS5 = NA, CS6 = NA, CS7 = NA, CS8 = NA,
                                                     CS9 = NA, CS10 = NA, CS11 = NA, CS12 = NA,
                                                     CS13 = NA, CS14 = NA, CS15 = NA, CS16 = NA,
                                                     Bivalent = NA, Active = NA, Divergent = NA, Heterochromatin = NA, Quies = NA)

List.CS.Simmilarity.LECIF.PhyloP.OsAt <- list(After_5050 = CS.Simmilarity.LECIF.PhyloP.OsAt.After,
                                              After_6040 = CS.Simmilarity.LECIF.PhyloP.OsAt.After,
                                              After_7030 = CS.Simmilarity.LECIF.PhyloP.OsAt.After,
                                              After_8020 = CS.Simmilarity.LECIF.PhyloP.OsAt.After,
                                              After_9010 = CS.Simmilarity.LECIF.PhyloP.OsAt.After)
bincolumns <- c(8:12)
for(z in 1:length(bincolumns)){ # 4 different cutoffs
  cat(paste0(colnames(OsAt.LECIFscore.PhyloPsinglent.mean)[bincolumns[z]], " cutoff processing ... \n"))
  for(i in 1:21){# 16 CS and 5 CS functional groups
    cat(paste0(i / 21 * 100, ' % CS groups started \n'))
    for(j in 1:4){ # lecif score percentile bins
      cat(paste0(j / 4 * 100, ' % cutoffs bins started \n'))
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], CS.At = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i])), method = 1)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+4,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], CS.At = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i])), method = 2) 
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+8,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], CS.At = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i])), method = 3)  
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+12,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], CS.At = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i])), method = 4) 
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+16,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], CS.At = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i])), method = 5) 
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+20,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], CS.At = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i])), method = 6) 
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+24,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], CS.At = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i])), method = 7) 
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+28,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], CS.At = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i])), method = 8) 
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+32,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], CS.At = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i])), method = 9) 
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+36,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], CS.At = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i])), method = 10) 
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+40,2+i] <- philentropy::euclidean(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F) 
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+44,2+i] <- philentropy::sorensen(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+48,2+i] <- philentropy::canberra(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+52,2+i] <- philentropy::wave_hedges(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+56,2+i] <- philentropy::tanimoto(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+60,2+i] <- philentropy::fidelity(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+64,2+i] <- philentropy::squared_chord(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+68,2+i] <- philentropy::squared_chi_sq(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+72,2+i] <- philentropy::additive_symm_chi_sq(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+76,2+i] <- philentropy::topsoe(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F,unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+80,2+i] <- philentropy::kumar_johnson(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F, epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+84,2+i] <- philentropy::manhattan(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+88,2+i] <- philentropy::lorentzian(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F, unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+92,2+i] <- philentropy::czekanowski(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+96,2+i] <- philentropy::kumar_hassebrook(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+100,2+i] <- philentropy::squared_euclidean(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+104,2+i] <- philentropy::prob_symm_chi_sq(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+108,2+i] <- philentropy::kullback_leibler_distance(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F, unit = "log2", epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+112,2+i] <- philentropy::jensen_shannon(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F, unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+116,2+i] <- philentropy::minkowski(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F, n = 3)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+120,2+i] <- philentropy::soergel(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+124,2+i] <- philentropy::inner_product(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+128,2+i] <- philentropy::hellinger(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+132,2+i] <- philentropy::divergence_sq(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+136,2+i] <- philentropy::jeffreys(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F, unit = "log2", epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+140,2+i] <- philentropy::jensen_difference(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F, unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+144,2+i] <- philentropy::chebyshev(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+148,2+i] <- philentropy::kulczynski_d(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F, epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+152,2+i] <- philentropy::harmonic_mean_dist(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+156,2+i] <- philentropy::matusita(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+160,2+i] <- philentropy::neyman_chi_sq(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F, epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+164,2+i] <- philentropy::clark_sq(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+168,2+i] <- philentropy::k_divergence(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F, unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+172,2+i] <- philentropy::taneja(P = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], Q = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i], testNA = F, unit = "log2", epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+176,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Os = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], CS.At = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i])), method = "maximum", margin = 2)[1,2]
      List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j+180,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Os = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),12+i], CS.At = OsAt.LECIFscore.PhyloPsinglent.mean[which(OsAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsAt[[z]][j,1]),33+i])), method = "hamming", margin = 2)[1,2]
      cat(paste0(j / 4 * 100, " % cutoff bins completed \n"))
    }
    cat(paste0(i / 21 * 100, " % CS groups completed \n"))
  }
  cat(paste0(colnames(OsAt.LECIFscore.PhyloPsinglent.mean)[bincolumns[z]], " cutoff computation completed \n"))
}

write.table(List.CS.Simmilarity.LECIF.PhyloP.OsAt$After_5050, file = "H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/OsAt_CSsimmilarity.DF.Distances.After5050.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(List.CS.Simmilarity.LECIF.PhyloP.OsAt$After_6040, file = "H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/OsAt_CSsimmilarity.DF.Distances.After6040.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(List.CS.Simmilarity.LECIF.PhyloP.OsAt$After_7030, file = "H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/OsAt_CSsimmilarity.DF.Distances.After7030.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(List.CS.Simmilarity.LECIF.PhyloP.OsAt$After_8020, file = "H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/OsAt_CSsimmilarity.DF.Distances.After8020.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(List.CS.Simmilarity.LECIF.PhyloP.OsAt$After_9010, file = "H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/OsAt_CSsimmilarity.DF.Distances.After9010.txt", col.names = T, row.names = F, sep = "\t", quote = F)


###############################################
############# 1.3.4 Oryza - Zea  ############# 
##############################################

OsZm.LECIFscore.PhyloPsinglent.mean <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/Os-Zm_LECIFscore_PhyloP_singlent_100overlap_mean.bedgraph", header = F, sep = "\t")
OsZm.df.F3 <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/OsZm_CSsimmilarity.DF3.txt", header = T, sep = "\t")

colnames(OsZm.LECIFscore.PhyloPsinglent.mean) <- c("Chr", "Start", "End", "LECIF", "PhyloP")
OsZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score <- percent_rank(OsZm.LECIFscore.PhyloPsinglent.mean$LECIF)
OsZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score <- percent_rank(OsZm.LECIFscore.PhyloPsinglent.mean$PhyloP)

# After_50;50
OsZm.LECIFscore.PhyloPsinglent.mean$After_5050 <- NA
# low,low 251031 regions
OsZm.LECIFscore.PhyloPsinglent.mean$After_5050[which(OsZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.5 & OsZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.5)] <- "Low,Low"
# high,low 206911 regions
OsZm.LECIFscore.PhyloPsinglent.mean$After_5050[which(OsZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.5 & OsZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.5)] <- "High,Low"
# high,high 251018 regions
OsZm.LECIFscore.PhyloPsinglent.mean$After_5050[which(OsZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.5 & OsZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.5)] <- "High,High"
# low,high 206924 regions
OsZm.LECIFscore.PhyloPsinglent.mean$After_5050[which(OsZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.5 & OsZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.5)] <- "Low,High"
(251031 + 206911 + 251018 + 206924) == nrow(OsZm.LECIFscore.PhyloPsinglent.mean) # check
length(which(is.na(OsZm.LECIFscore.PhyloPsinglent.mean$After_5050) == TRUE)) == 0 # check

# After_60;40
OsZm.LECIFscore.PhyloPsinglent.mean$After_6040 <- NA
# low,low 171623 regions
OsZm.LECIFscore.PhyloPsinglent.mean$After_6040[which(OsZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.4 & OsZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.4)] <- "Low,Low"
# high,low 125205 regions
OsZm.LECIFscore.PhyloPsinglent.mean$After_6040[which(OsZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.6 & OsZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.4)] <- "High,Low"
# high,high 161670 regions
OsZm.LECIFscore.PhyloPsinglent.mean$After_6040[which(OsZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.6 & OsZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.6)] <- "High,High"
# low,high 127028 regions
OsZm.LECIFscore.PhyloPsinglent.mean$After_6040[which(OsZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.4 & OsZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.6)] <- "Low,High"
length(which(is.na(OsZm.LECIFscore.PhyloPsinglent.mean$After_6040) == TRUE)) # check
(171623+125205+161670+127028+330358) == nrow(OsZm.LECIFscore.PhyloPsinglent.mean)

# After_70;30
OsZm.LECIFscore.PhyloPsinglent.mean$After_7030 <- NA
# low,low 103207 regions
OsZm.LECIFscore.PhyloPsinglent.mean$After_7030[which(OsZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.3 & OsZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.3)] <- "Low,Low"
# high,low 69848 regions
OsZm.LECIFscore.PhyloPsinglent.mean$After_7030[which(OsZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.7 & OsZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.3)] <- "High,Low"
# high,high 89650 regions
OsZm.LECIFscore.PhyloPsinglent.mean$After_7030[which(OsZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.7 & OsZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.7)] <- "High,High"
# low,high 67280 regions
OsZm.LECIFscore.PhyloPsinglent.mean$After_7030[which(OsZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.3 & OsZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.7)] <- "Low,High"
length(which(is.na(OsZm.LECIFscore.PhyloPsinglent.mean$After_7030) == TRUE)) # check
(103207+69848+89650+67280+585899) == nrow(OsZm.LECIFscore.PhyloPsinglent.mean)

# After_80;20
OsZm.LECIFscore.PhyloPsinglent.mean$After_8020 <- NA
# low,low 49854 regions
OsZm.LECIFscore.PhyloPsinglent.mean$After_8020[which(OsZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.2 & OsZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.2)] <- "Low,Low"
# high,low 32793 regions
OsZm.LECIFscore.PhyloPsinglent.mean$After_8020[which(OsZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.8 & OsZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.2)] <- "High,Low"
# high,high 38945 regions
OsZm.LECIFscore.PhyloPsinglent.mean$After_8020[which(OsZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.8 & OsZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.8)] <- "High,High"
# low,high 27647 regions
OsZm.LECIFscore.PhyloPsinglent.mean$After_8020[which(OsZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.2 & OsZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.8)] <- "Low,High"
length(which(is.na(OsZm.LECIFscore.PhyloPsinglent.mean$After_8020) == TRUE)) # check
(49854+32793+38945+27647+766645) == nrow(OsZm.LECIFscore.PhyloPsinglent.mean)

# After_90;10
OsZm.LECIFscore.PhyloPsinglent.mean$After_9010 <- NA
# low,low 14737 regions
OsZm.LECIFscore.PhyloPsinglent.mean$After_9010[which(OsZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.1 & OsZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.1)] <- "Low,Low"
# high,low 8508 regions
OsZm.LECIFscore.PhyloPsinglent.mean$After_9010[which(OsZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.9 & OsZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.1)] <- "High,Low"
# high,high 9636 regions
OsZm.LECIFscore.PhyloPsinglent.mean$After_9010[which(OsZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.9 & OsZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.9)] <- "High,High"
# low,high 5855 regions
OsZm.LECIFscore.PhyloPsinglent.mean$After_9010[which(OsZm.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.1 & OsZm.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.9)] <- "Low,High"
length(which(is.na(OsZm.LECIFscore.PhyloPsinglent.mean$After_9010) == TRUE)) # check
(14737+8508+9636+5855+877148) == nrow(OsZm.LECIFscore.PhyloPsinglent.mean)

## CS.At and CS.Zm columns

OsZm.LECIFscore.PhyloPsinglent.mean$CS1.Os <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS2.Os <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS3.Os <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS4.Os <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS5.Os <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS6.Os <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS7.Os <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS8.Os <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS9.Os <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS10.Os <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS11.Os <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS12.Os <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS13.Os <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS14.Os <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS15.Os <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS16.Os <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$Bivalent.Os <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$Active.Os <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$Divergent.Os <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$Heterochromatin.Os <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$Quies.Os <- NA

OsZm.LECIFscore.PhyloPsinglent.mean$CS1.Zm <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS2.Zm <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS3.Zm <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS4.Zm <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS5.Zm <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS6.Zm <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS7.Zm <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS8.Zm <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS9.Zm <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS10.Zm <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS11.Zm <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS12.Zm <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS13.Zm <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS14.Zm <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS15.Zm <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$CS16.Zm <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$Bivalent.Zm <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$Active.Zm <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$Divergent.Zm <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$Heterochromatin.Zm <- NA
OsZm.LECIFscore.PhyloPsinglent.mean$Quies.Zm <- NA

## loop for annotating my new merged regions with both scores with previously one with CS presence/absence info

for(i in 1:nrow(OsZm.LECIFscore.PhyloPsinglent.mean)){
  print(paste0(i/nrow(OsZm.LECIFscore.PhyloPsinglent.mean)*100, " % covered \n"))
  hit <- which(OsZm.df.F3$Chr.Os == OsZm.LECIFscore.PhyloPsinglent.mean$Chr[i] & OsZm.df.F3$Start.Os == OsZm.LECIFscore.PhyloPsinglent.mean$Start[i] & OsZm.df.F3$End.Os == OsZm.LECIFscore.PhyloPsinglent.mean$End[i])
  if(length(hit) == 1){
    OsZm.LECIFscore.PhyloPsinglent.mean[i,c(13:54)] <- OsZm.df.F3[hit,10:51]
  }else if(length(hit) == 0){
    Chr <- which(OsZm.df.F3$Chr.Os == OsZm.LECIFscore.PhyloPsinglent.mean$Chr[i])
    hit2 <- which(OsZm.df.F3$Start.Os[Chr] >= OsZm.LECIFscore.PhyloPsinglent.mean$Start[i] & OsZm.df.F3$End.Os[Chr] <= OsZm.LECIFscore.PhyloPsinglent.mean$End[i])
    if(length(hit2) == 0){stop("Something go wrong boy ... \n")}
    OsZm.LECIFscore.PhyloPsinglent.mean[i,c(13:54)] <- colSums(OsZm.df.F3[hit2,10:51])
    colto1 <- c(13:54)[which(OsZm.LECIFscore.PhyloPsinglent.mean[i,c(13:54)] > 1)]
    OsZm.LECIFscore.PhyloPsinglent.mean[i, colto1] <- 1
  }else {
    stop("Something go wrong boy ... \n")
  }
}

write.table(OsZm.LECIFscore.PhyloPsinglent.mean, "H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/OsZm_CSsimmilarity.DF.LECIFPhyloP.AfterCutOff.txt", col.names = T, row.names = F, sep = "\t", quote = F)

## loop for computing simmilarity/distance between different subsets of LECIF and PhyloP score with different cutoffs (After ...)

CS.Simmilarity.LECIF.PhyloP.OsZm.After <- data.frame(Bin = rep(c("High,High", "Low,Low", "High,Low", "Low,High"),46),
                                                     Metric = c(rep(c("Jaccard"),4),rep(c("SMC"),4),rep(c("Sokal"),4), rep(c("Tanimoto"),4),
                                                                rep(c("Dice"),4), rep(c("Hanman"),4), rep(c("Ochiai"),4), rep(c("Sokal2"),4),
                                                                rep(c("Phi"),4), rep(c("S2"),4), 
                                                                rep(c("Eucledian"),4), rep(c("Sorensen"),4), rep(c("Canberra"),4), rep(c("Wavehedges"),4), 
                                                                rep(c("Tanimoto2"),4), rep(c("Fidelity"),4), rep(c("SquaredChord"),4), rep(c("SquaredChi"),4), 
                                                                rep(c("AdditiveSym"),4), rep(c("Topsoe"),4), rep(c("Kumarjhonson"),4), rep(c("Manhattan"),4),
                                                                rep(c("Lorentzian"),4), rep(c("Czekanowski"),4), rep(c("Hassebrook"),4),  rep(c("SquaredEucledian"),4),
                                                                rep(c("ProbSymm"),4), rep(c("Kullbackleiber"),4), rep(c("Jensenshannon"),4), rep(c("Minkowski"),4), 
                                                                rep(c("Soergel"),4), rep(c("Innerproduct"),4), rep(c("Hellinger"),4), rep(c("Divergence"),4), 
                                                                rep(c("Jeffreys"),4), rep(c("Jensendifference"),4), rep(c("Chebysev"),4), rep(c("Kulczynskid"),4), 
                                                                rep(c("Harmonicmean"),4), rep(c("Matusita"),4), rep(c("Neyman"),4), rep(c("Clark"),4), 
                                                                rep(c("Kdivergence"),4), rep(c("Taneja"),4), 
                                                                rep(c("Maximum"),4), rep(c("Hamming"),4)),
                                                     CS1 = NA, CS2 = NA, CS3 = NA, CS4 = NA,
                                                     CS5 = NA, CS6 = NA, CS7 = NA, CS8 = NA,
                                                     CS9 = NA, CS10 = NA, CS11 = NA, CS12 = NA,
                                                     CS13 = NA, CS14 = NA, CS15 = NA, CS16 = NA,
                                                     Bivalent = NA, Active = NA, Divergent = NA, Heterochromatin = NA, Quies = NA)

List.CS.Simmilarity.LECIF.PhyloP.OsZm <- list(After_5050 = CS.Simmilarity.LECIF.PhyloP.OsZm.After,
                                              After_6040 = CS.Simmilarity.LECIF.PhyloP.OsZm.After,
                                              After_7030 = CS.Simmilarity.LECIF.PhyloP.OsZm.After,
                                              After_8020 = CS.Simmilarity.LECIF.PhyloP.OsZm.After,
                                              After_9010 = CS.Simmilarity.LECIF.PhyloP.OsZm.After)
bincolumns <- c(8:12)
for(z in 1:length(bincolumns)){ # 4 different cutoffs
  cat(paste0(colnames(OsZm.LECIFscore.PhyloPsinglent.mean)[bincolumns[z]], " cutoff processing ... \n"))
  for(i in 1:21){# 16 CS and 5 CS functional groups
    cat(paste0(i / 21 * 100, ' % CS groups started \n'))
    for(j in 1:4){ # lecif score percentile bins
      cat(paste0(j / 4 * 100, ' % cutoffs bins started \n'))
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], CS.Zm = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i])), method = 1)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+4,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], CS.Zm = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i])), method = 2) 
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+8,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], CS.Zm = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i])), method = 3)  
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+12,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], CS.Zm = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i])), method = 4) 
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+16,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], CS.Zm = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i])), method = 5) 
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+20,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], CS.Zm = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i])), method = 6) 
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+24,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], CS.Zm = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i])), method = 7) 
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+28,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], CS.Zm = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i])), method = 8) 
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+32,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], CS.Zm = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i])), method = 9) 
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+36,2+i] <- ade4::dist.binary(t(data.frame(CS.Os = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], CS.Zm = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i])), method = 10) 
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+40,2+i] <- philentropy::euclidean(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F) 
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+44,2+i] <- philentropy::sorensen(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+48,2+i] <- philentropy::canberra(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+52,2+i] <- philentropy::wave_hedges(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+56,2+i] <- philentropy::tanimoto(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+60,2+i] <- philentropy::fidelity(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+64,2+i] <- philentropy::squared_chord(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+68,2+i] <- philentropy::squared_chi_sq(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+72,2+i] <- philentropy::additive_symm_chi_sq(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+76,2+i] <- philentropy::topsoe(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F,unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+80,2+i] <- philentropy::kumar_johnson(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F, epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+84,2+i] <- philentropy::manhattan(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+88,2+i] <- philentropy::lorentzian(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F, unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+92,2+i] <- philentropy::czekanowski(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+96,2+i] <- philentropy::kumar_hassebrook(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+100,2+i] <- philentropy::squared_euclidean(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+104,2+i] <- philentropy::prob_symm_chi_sq(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+108,2+i] <- philentropy::kullback_leibler_distance(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F, unit = "log2", epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+112,2+i] <- philentropy::jensen_shannon(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F, unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+116,2+i] <- philentropy::minkowski(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F, n = 3)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+120,2+i] <- philentropy::soergel(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+124,2+i] <- philentropy::inner_product(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+128,2+i] <- philentropy::hellinger(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+132,2+i] <- philentropy::divergence_sq(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+136,2+i] <- philentropy::jeffreys(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F, unit = "log2", epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+140,2+i] <- philentropy::jensen_difference(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F, unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+144,2+i] <- philentropy::chebyshev(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+148,2+i] <- philentropy::kulczynski_d(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F, epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+152,2+i] <- philentropy::harmonic_mean_dist(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+156,2+i] <- philentropy::matusita(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+160,2+i] <- philentropy::neyman_chi_sq(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F, epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+164,2+i] <- philentropy::clark_sq(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+168,2+i] <- philentropy::k_divergence(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F, unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+172,2+i] <- philentropy::taneja(P = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], Q = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i], testNA = F, unit = "log2", epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+176,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Os = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], CS.Zm = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i])), method = "maximum", margin = 2)[1,2]
      List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j+180,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Os = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),12+i], CS.Zm = OsZm.LECIFscore.PhyloPsinglent.mean[which(OsZm.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.OsZm[[z]][j,1]),33+i])), method = "hamming", margin = 2)[1,2]
      cat(paste0(j / 4 * 100, " % cutoff bins completed \n"))
    }
    cat(paste0(i / 21 * 100, " % CS groups completed \n"))
  }
  cat(paste0(colnames(OsZm.LECIFscore.PhyloPsinglent.mean)[bincolumns[z]], " cutoff computation completed \n"))
}

write.table(List.CS.Simmilarity.LECIF.PhyloP.OsZm$After_5050, file = "H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/OsZm_CSsimmilarity.DF.Distances.After5050.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(List.CS.Simmilarity.LECIF.PhyloP.OsZm$After_6040, file = "H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/OsZm_CSsimmilarity.DF.Distances.After6040.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(List.CS.Simmilarity.LECIF.PhyloP.OsZm$After_7030, file = "H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/OsZm_CSsimmilarity.DF.Distances.After7030.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(List.CS.Simmilarity.LECIF.PhyloP.OsZm$After_8020, file = "H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/OsZm_CSsimmilarity.DF.Distances.After8020.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(List.CS.Simmilarity.LECIF.PhyloP.OsZm$After_9010, file = "H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/OsZm_CSsimmilarity.DF.Distances.After9010.txt", col.names = T, row.names = F, sep = "\t", quote = F)

###############################################
########## 1.3.5 Zea - Arabidopsis  ###########
##############################################

## Condition PhyloP scores for all 100% bases in LECIF region

ZmAt.LECIFscore.PhyloPsinglent.mean <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/Zm-At_LECIFscore_PhyloP_singlent_100overlap_mean.bedgraph", header = F, sep = "\t")
ZmAt.df.F3 <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/ZmAt_CSsimmilarity.DF3.txt", header = T, sep = "\t")

colnames(ZmAt.LECIFscore.PhyloPsinglent.mean) <- c("Chr", "Start", "End", "LECIF", "PhyloP")
ZmAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score <- percent_rank(ZmAt.LECIFscore.PhyloPsinglent.mean$LECIF)
ZmAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score <- percent_rank(ZmAt.LECIFscore.PhyloPsinglent.mean$PhyloP)

# After_50;50
ZmAt.LECIFscore.PhyloPsinglent.mean$After_5050 <- NA
# low,low 36617 regions
ZmAt.LECIFscore.PhyloPsinglent.mean$After_5050[which(ZmAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.5 & ZmAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.5)] <- "Low,Low"
# high,low 36439 regions
ZmAt.LECIFscore.PhyloPsinglent.mean$After_5050[which(ZmAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.5 & ZmAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.5)] <- "High,Low"
# high,high 36611 regions
ZmAt.LECIFscore.PhyloPsinglent.mean$After_5050[which(ZmAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.5 & ZmAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.5)] <- "High,High"
# low,high 36444 regions
ZmAt.LECIFscore.PhyloPsinglent.mean$After_5050[which(ZmAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.5 & ZmAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.5)] <- "Low,High"
(36617 + 36439 + 36611 + 36444) == nrow(ZmAt.LECIFscore.PhyloPsinglent.mean) # check
length(which(is.na(ZmAt.LECIFscore.PhyloPsinglent.mean$After_5050) == TRUE)) == 0 # check

# After_60;40
ZmAt.LECIFscore.PhyloPsinglent.mean$After_6040 <- NA
# low,low 23695 regions
ZmAt.LECIFscore.PhyloPsinglent.mean$After_6040[which(ZmAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.4 & ZmAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.4)] <- "Low,Low"
# high,low 23123 regions
ZmAt.LECIFscore.PhyloPsinglent.mean$After_6040[which(ZmAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.6 & ZmAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.4)] <- "High,Low"
# high,high 23471 regions
ZmAt.LECIFscore.PhyloPsinglent.mean$After_6040[which(ZmAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.6 & ZmAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.6)] <- "High,High"
# low,high 23170 regions
ZmAt.LECIFscore.PhyloPsinglent.mean$After_6040[which(ZmAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.4 & ZmAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.6)] <- "Low,High"
length(which(is.na(ZmAt.LECIFscore.PhyloPsinglent.mean$After_6040) == TRUE)) # check
(23695+23123+23471+23170+52652) == nrow(ZmAt.LECIFscore.PhyloPsinglent.mean)

# After_70;30
ZmAt.LECIFscore.PhyloPsinglent.mean$After_7030 <- NA
# low,low 13440 regions
ZmAt.LECIFscore.PhyloPsinglent.mean$After_7030[which(ZmAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.3 & ZmAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.3)] <- "Low,Low"
# high,low 12911 regions
ZmAt.LECIFscore.PhyloPsinglent.mean$After_7030[which(ZmAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.7 & ZmAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.3)] <- "High,Low"
# high,high 13158 regions
ZmAt.LECIFscore.PhyloPsinglent.mean$After_7030[which(ZmAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.7 & ZmAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.7)] <- "High,High"
# low,high 13089 regions
ZmAt.LECIFscore.PhyloPsinglent.mean$After_7030[which(ZmAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.3 & ZmAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.7)] <- "Low,High"
length(which(is.na(ZmAt.LECIFscore.PhyloPsinglent.mean$After_7030) == TRUE)) # check
(13440+12911+13158+13089+93513) == nrow(ZmAt.LECIFscore.PhyloPsinglent.mean)

# After_80;20
ZmAt.LECIFscore.PhyloPsinglent.mean$After_8020 <- NA
# low,low 6000 regions
ZmAt.LECIFscore.PhyloPsinglent.mean$After_8020[which(ZmAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.2 & ZmAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.2)] <- "Low,Low"
# high,low 5641 regions
ZmAt.LECIFscore.PhyloPsinglent.mean$After_8020[which(ZmAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.8 & ZmAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.2)] <- "High,Low"
# high,high 5722 regions
ZmAt.LECIFscore.PhyloPsinglent.mean$After_8020[which(ZmAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.8 & ZmAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.8)] <- "High,High"
# low,high 5834 regions
ZmAt.LECIFscore.PhyloPsinglent.mean$After_8020[which(ZmAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.2 & ZmAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.8)] <- "Low,High"
length(which(is.na(ZmAt.LECIFscore.PhyloPsinglent.mean$After_8020) == TRUE)) # check
(6000+5641+5722+5834+122914) == nrow(ZmAt.LECIFscore.PhyloPsinglent.mean)

# After_90;10
ZmAt.LECIFscore.PhyloPsinglent.mean$After_9010 <- NA
# low,low 1599 regions
ZmAt.LECIFscore.PhyloPsinglent.mean$After_9010[which(ZmAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.1 & ZmAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.1)] <- "Low,Low"
# high,low 1355 regions
ZmAt.LECIFscore.PhyloPsinglent.mean$After_9010[which(ZmAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.9 & ZmAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.1)] <- "High,Low"
# high,high 1363 regions
ZmAt.LECIFscore.PhyloPsinglent.mean$After_9010[which(ZmAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.9 & ZmAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.9)] <- "High,High"
# low,high 1376 regions
ZmAt.LECIFscore.PhyloPsinglent.mean$After_9010[which(ZmAt.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.1 & ZmAt.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.9)] <- "Low,High"
length(which(is.na(ZmAt.LECIFscore.PhyloPsinglent.mean$After_9010) == TRUE)) # check
(1599+1355+1363+1376+140418) == nrow(ZmAt.LECIFscore.PhyloPsinglent.mean)

## CS.At and CS.Zm columns

ZmAt.LECIFscore.PhyloPsinglent.mean$CS1.Zm <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS2.Zm <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS3.Zm <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS4.Zm <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS5.Zm <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS6.Zm <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS7.Zm <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS8.Zm <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS9.Zm <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS10.Zm <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS11.Zm <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS12.Zm <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS13.Zm <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS14.Zm <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS15.Zm <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS16.Zm <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$Bivalent.Zm <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$Active.Zm <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$Divergent.Zm <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$Heterochromatin.Zm <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$Quies.Zm <- NA

ZmAt.LECIFscore.PhyloPsinglent.mean$CS1.At <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS2.At <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS3.At <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS4.At <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS5.At <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS6.At <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS7.At <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS8.At <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS9.At <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS10.At <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS11.At <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS12.At <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS13.At <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS14.At <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS15.At <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$CS16.At <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$Bivalent.At <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$Active.At <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$Divergent.At <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$Heterochromatin.At <- NA
ZmAt.LECIFscore.PhyloPsinglent.mean$Quies.At <- NA

## loop for annotating my new merged regions with both scores with previously one with CS presence/absence info

for(i in 1:nrow(ZmAt.LECIFscore.PhyloPsinglent.mean)){
  print(paste0(i/nrow(ZmAt.LECIFscore.PhyloPsinglent.mean)*100, " % covered \n"))
  hit <- which(ZmAt.df.F3$Chr.Zm == ZmAt.LECIFscore.PhyloPsinglent.mean$Chr[i] & ZmAt.df.F3$Start.Zm == ZmAt.LECIFscore.PhyloPsinglent.mean$Start[i] & ZmAt.df.F3$End.Zm == ZmAt.LECIFscore.PhyloPsinglent.mean$End[i])
  if(length(hit) == 1){
    ZmAt.LECIFscore.PhyloPsinglent.mean[i,c(13:54)] <- ZmAt.df.F3[hit,10:51]
  }else if(length(hit) == 0){
    Chr <- which(ZmAt.df.F3$Chr.Zm == ZmAt.LECIFscore.PhyloPsinglent.mean$Chr[i])
    hit2 <- which(ZmAt.df.F3$Start.Zm[Chr] >= ZmAt.LECIFscore.PhyloPsinglent.mean$Start[i] & ZmAt.df.F3$End.Zm[Chr] <= ZmAt.LECIFscore.PhyloPsinglent.mean$End[i])
    if(length(hit2) == 0){stop("Something go wrong boy ... \n")}
    ZmAt.LECIFscore.PhyloPsinglent.mean[i,c(13:54)] <- colSums(ZmAt.df.F3[hit2,10:51])
    colto1 <- c(13:54)[which(ZmAt.LECIFscore.PhyloPsinglent.mean[i,c(13:54)] > 1)]
    ZmAt.LECIFscore.PhyloPsinglent.mean[i, colto1] <- 1
  }else {
    stop("Something go wrong boy ... \n")
  }
}

write.table(ZmAt.LECIFscore.PhyloPsinglent.mean, "H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/ZmAt_CSsimmilarity.DF.LECIFPhyloP.AfterCutOff.txt", col.names = T, row.names = F, sep = "\t", quote = F)

## loop for computing simmilarity/distance between different subsets of LECIF and PhyloP score with different cutoffs (After ...)

CS.Simmilarity.LECIF.PhyloP.ZmAt.After <- data.frame(Bin = rep(c("High,High", "Low,Low", "High,Low", "Low,High"),46),
                                                     Metric = c(rep(c("Jaccard"),4),rep(c("SMC"),4),rep(c("Sokal"),4), rep(c("Tanimoto"),4),
                                                                rep(c("Dice"),4), rep(c("Hanman"),4), rep(c("Ochiai"),4), rep(c("Sokal2"),4),
                                                                rep(c("Phi"),4), rep(c("S2"),4), 
                                                                rep(c("Eucledian"),4), rep(c("Sorensen"),4), rep(c("Canberra"),4), rep(c("Wavehedges"),4), 
                                                                rep(c("Tanimoto2"),4), rep(c("Fidelity"),4), rep(c("SquaredChord"),4), rep(c("SquaredChi"),4), 
                                                                rep(c("AdditiveSym"),4), rep(c("Topsoe"),4), rep(c("Kumarjhonson"),4), rep(c("Manhattan"),4),
                                                                rep(c("Lorentzian"),4), rep(c("Czekanowski"),4), rep(c("Hassebrook"),4),  rep(c("SquaredEucledian"),4),
                                                                rep(c("ProbSymm"),4), rep(c("Kullbackleiber"),4), rep(c("Jensenshannon"),4), rep(c("Minkowski"),4), 
                                                                rep(c("Soergel"),4), rep(c("Innerproduct"),4), rep(c("Hellinger"),4), rep(c("Divergence"),4), 
                                                                rep(c("Jeffreys"),4), rep(c("Jensendifference"),4), rep(c("Chebysev"),4), rep(c("Kulczynskid"),4), 
                                                                rep(c("Harmonicmean"),4), rep(c("Matusita"),4), rep(c("Neyman"),4), rep(c("Clark"),4), 
                                                                rep(c("Kdivergence"),4), rep(c("Taneja"),4), 
                                                                rep(c("Maximum"),4), rep(c("Hamming"),4)),
                                                     CS1 = NA, CS2 = NA, CS3 = NA, CS4 = NA,
                                                     CS5 = NA, CS6 = NA, CS7 = NA, CS8 = NA,
                                                     CS9 = NA, CS10 = NA, CS11 = NA, CS12 = NA,
                                                     CS13 = NA, CS14 = NA, CS15 = NA, CS16 = NA,
                                                     Bivalent = NA, Active = NA, Divergent = NA, Heterochromatin = NA, Quies = NA)

List.CS.Simmilarity.LECIF.PhyloP.ZmAt <- list(After_5050 = CS.Simmilarity.LECIF.PhyloP.ZmAt.After,
                                              After_6040 = CS.Simmilarity.LECIF.PhyloP.ZmAt.After,
                                              After_7030 = CS.Simmilarity.LECIF.PhyloP.ZmAt.After,
                                              After_8020 = CS.Simmilarity.LECIF.PhyloP.ZmAt.After,
                                              After_9010 = CS.Simmilarity.LECIF.PhyloP.ZmAt.After)
bincolumns <- c(8:12)
for(z in 1:length(bincolumns)){ # 4 different cutoffs
  cat(paste0(colnames(ZmAt.LECIFscore.PhyloPsinglent.mean)[bincolumns[z]], " cutoff processing ... \n"))
  for(i in 1:21){# 16 CS and 5 CS functional groups
    cat(paste0(i / 21 * 100, ' % CS groups started \n'))
    for(j in 1:4){ # lecif score percentile bins
      cat(paste0(j / 4 * 100, ' % cutoffs bins started \n'))
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], CS.At = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i])), method = 1)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+4,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], CS.At = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i])), method = 2) 
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+8,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], CS.At = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i])), method = 3)  
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+12,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], CS.At = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i])), method = 4) 
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+16,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], CS.At = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i])), method = 5) 
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+20,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], CS.At = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i])), method = 6) 
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+24,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], CS.At = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i])), method = 7) 
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+28,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], CS.At = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i])), method = 8) 
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+32,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], CS.At = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i])), method = 9) 
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+36,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], CS.At = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i])), method = 10) 
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+40,2+i] <- philentropy::euclidean(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F) 
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+44,2+i] <- philentropy::sorensen(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+48,2+i] <- philentropy::canberra(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+52,2+i] <- philentropy::wave_hedges(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+56,2+i] <- philentropy::tanimoto(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+60,2+i] <- philentropy::fidelity(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+64,2+i] <- philentropy::squared_chord(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+68,2+i] <- philentropy::squared_chi_sq(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+72,2+i] <- philentropy::additive_symm_chi_sq(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+76,2+i] <- philentropy::topsoe(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F,unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+80,2+i] <- philentropy::kumar_johnson(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F, epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+84,2+i] <- philentropy::manhattan(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+88,2+i] <- philentropy::lorentzian(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F, unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+92,2+i] <- philentropy::czekanowski(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+96,2+i] <- philentropy::kumar_hassebrook(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+100,2+i] <- philentropy::squared_euclidean(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+104,2+i] <- philentropy::prob_symm_chi_sq(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+108,2+i] <- philentropy::kullback_leibler_distance(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F, unit = "log2", epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+112,2+i] <- philentropy::jensen_shannon(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F, unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+116,2+i] <- philentropy::minkowski(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F, n = 3)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+120,2+i] <- philentropy::soergel(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+124,2+i] <- philentropy::inner_product(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+128,2+i] <- philentropy::hellinger(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+132,2+i] <- philentropy::divergence_sq(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+136,2+i] <- philentropy::jeffreys(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F, unit = "log2", epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+140,2+i] <- philentropy::jensen_difference(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F, unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+144,2+i] <- philentropy::chebyshev(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+148,2+i] <- philentropy::kulczynski_d(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F, epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+152,2+i] <- philentropy::harmonic_mean_dist(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+156,2+i] <- philentropy::matusita(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+160,2+i] <- philentropy::neyman_chi_sq(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F, epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+164,2+i] <- philentropy::clark_sq(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+168,2+i] <- philentropy::k_divergence(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F, unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+172,2+i] <- philentropy::taneja(P = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], Q = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i], testNA = F, unit = "log2", epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+176,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Zm = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], CS.At = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i])), method = "maximum", margin = 2)[1,2]
      List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j+180,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Zm = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),12+i], CS.At = ZmAt.LECIFscore.PhyloPsinglent.mean[which(ZmAt.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmAt[[z]][j,1]),33+i])), method = "hamming", margin = 2)[1,2]
      cat(paste0(j / 4 * 100, " % cutoff bins completed \n"))
    }
    cat(paste0(i / 21 * 100, " % CS groups completed \n"))
  }
  cat(paste0(colnames(ZmAt.LECIFscore.PhyloPsinglent.mean)[bincolumns[z]], " cutoff computation completed \n"))
}

write.table(List.CS.Simmilarity.LECIF.PhyloP.ZmAt$After_5050, file = "H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/ZmAt_CSsimmilarity.DF.Distances.After5050.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(List.CS.Simmilarity.LECIF.PhyloP.ZmAt$After_6040, file = "H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/ZmAt_CSsimmilarity.DF.Distances.After6040.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(List.CS.Simmilarity.LECIF.PhyloP.ZmAt$After_7030, file = "H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/ZmAt_CSsimmilarity.DF.Distances.After7030.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(List.CS.Simmilarity.LECIF.PhyloP.ZmAt$After_8020, file = "H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/ZmAt_CSsimmilarity.DF.Distances.After8020.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(List.CS.Simmilarity.LECIF.PhyloP.ZmAt$After_9010, file = "H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/ZmAt_CSsimmilarity.DF.Distances.After9010.txt", col.names = T, row.names = F, sep = "\t", quote = F)

###############################################
############# 1.3.6 Zea - Oryza  ##############
##############################################

## Condition PhyloP scores for all 100% bases in LECIF region

ZmOs.LECIFscore.PhyloPsinglent.mean <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/Zm-Os_LECIFscore_PhyloP_singlent_100overlap_mean.bedgraph", header = F, sep = "\t")
ZmOs.df.F3 <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/ZmOs_CSsimmilarity.DF3.txt", header = T, sep = "\t")

colnames(ZmOs.LECIFscore.PhyloPsinglent.mean) <- c("Chr", "Start", "End", "LECIF", "PhyloP")
ZmOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score <- percent_rank(ZmOs.LECIFscore.PhyloPsinglent.mean$LECIF)
ZmOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score <- percent_rank(ZmOs.LECIFscore.PhyloPsinglent.mean$PhyloP)

# After_50;50
ZmOs.LECIFscore.PhyloPsinglent.mean$After_5050 <- NA
# low,low 87767 regions
ZmOs.LECIFscore.PhyloPsinglent.mean$After_5050[which(ZmOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.5 & ZmOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.5)] <- "Low,Low"
# high,low 83852 regions
ZmOs.LECIFscore.PhyloPsinglent.mean$After_5050[which(ZmOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.5 & ZmOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.5)] <- "High,Low"
# high,high 87767 regions
ZmOs.LECIFscore.PhyloPsinglent.mean$After_5050[which(ZmOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.5 & ZmOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.5)] <- "High,High"
# low,high 83852 regions
ZmOs.LECIFscore.PhyloPsinglent.mean$After_5050[which(ZmOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.5 & ZmOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.5)] <- "Low,High"
(87767 + 83852 + 87767 + 83852) == nrow(ZmOs.LECIFscore.PhyloPsinglent.mean) # check
length(which(is.na(ZmOs.LECIFscore.PhyloPsinglent.mean$After_5050) == TRUE)) == 0 # check

# After_60;40
ZmOs.LECIFscore.PhyloPsinglent.mean$After_6040 <- NA
# low,low 56131 regions
ZmOs.LECIFscore.PhyloPsinglent.mean$After_6040[which(ZmOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.4 & ZmOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.4)] <- "Low,Low"
# high,low 53189 regions
ZmOs.LECIFscore.PhyloPsinglent.mean$After_6040[which(ZmOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.6 & ZmOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.4)] <- "High,Low"
# high,high 57328 regions
ZmOs.LECIFscore.PhyloPsinglent.mean$After_6040[which(ZmOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.6 & ZmOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.6)] <- "High,High"
# low,high 53029 regions
ZmOs.LECIFscore.PhyloPsinglent.mean$After_6040[which(ZmOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.4 & ZmOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.6)] <- "Low,High"
length(which(is.na(ZmOs.LECIFscore.PhyloPsinglent.mean$After_6040) == TRUE)) # check
(56131+53189+57328+53029+123561) == nrow(ZmOs.LECIFscore.PhyloPsinglent.mean)

# After_70;30
ZmOs.LECIFscore.PhyloPsinglent.mean$After_7030 <- NA
# low,low 31551 regions
ZmOs.LECIFscore.PhyloPsinglent.mean$After_7030[which(ZmOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.3 & ZmOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.3)] <- "Low,Low"
# high,low 29634 regions
ZmOs.LECIFscore.PhyloPsinglent.mean$After_7030[which(ZmOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.7 & ZmOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.3)] <- "High,Low"
# high,high 32860 regions
ZmOs.LECIFscore.PhyloPsinglent.mean$After_7030[which(ZmOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.7 & ZmOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.7)] <- "High,High"
# low,high 29409 regions
ZmOs.LECIFscore.PhyloPsinglent.mean$After_7030[which(ZmOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.3 & ZmOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.7)] <- "Low,High"
length(which(is.na(ZmOs.LECIFscore.PhyloPsinglent.mean$After_7030) == TRUE)) # check
(31551+29634+32860+29409+219784) == nrow(ZmOs.LECIFscore.PhyloPsinglent.mean)

# After_80;20
ZmOs.LECIFscore.PhyloPsinglent.mean$After_8020 <- NA
# low,low 14054 regions
ZmOs.LECIFscore.PhyloPsinglent.mean$After_8020[which(ZmOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.2 & ZmOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.2)] <- "Low,Low"
# high,low 13209 regions
ZmOs.LECIFscore.PhyloPsinglent.mean$After_8020[which(ZmOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.8 & ZmOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.2)] <- "High,Low"
# high,high 14718 regions
ZmOs.LECIFscore.PhyloPsinglent.mean$After_8020[which(ZmOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.8 & ZmOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.8)] <- "High,High"
# low,high 13003 regions
ZmOs.LECIFscore.PhyloPsinglent.mean$After_8020[which(ZmOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.2 & ZmOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.8)] <- "Low,High"
length(which(is.na(ZmOs.LECIFscore.PhyloPsinglent.mean$After_8020) == TRUE)) # check
(14054+13209+14718+13003+288254) == nrow(ZmOs.LECIFscore.PhyloPsinglent.mean)

# After_90;10
ZmOs.LECIFscore.PhyloPsinglent.mean$After_9010 <- NA
# low,low 3546 regions
ZmOs.LECIFscore.PhyloPsinglent.mean$After_9010[which(ZmOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.1 & ZmOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.1)] <- "Low,Low"
# high,low 3225 regions
ZmOs.LECIFscore.PhyloPsinglent.mean$After_9010[which(ZmOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.9 & ZmOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score < 0.1)] <- "High,Low"
# high,high 3771 regions
ZmOs.LECIFscore.PhyloPsinglent.mean$After_9010[which(ZmOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score > 0.9 & ZmOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.9)] <- "High,High"
# low,high 3293 regions
ZmOs.LECIFscore.PhyloPsinglent.mean$After_9010[which(ZmOs.LECIFscore.PhyloPsinglent.mean$LECIF.Percentile.Rank.Score < 0.1 & ZmOs.LECIFscore.PhyloPsinglent.mean$PhyloP.Percentile.Rank.Score > 0.9)] <- "Low,High"
length(which(is.na(ZmOs.LECIFscore.PhyloPsinglent.mean$After_9010) == TRUE)) # check
(3546+3225+3771+3293+329403) == nrow(ZmOs.LECIFscore.PhyloPsinglent.mean)

## CS.At and CS.Os columns

ZmOs.LECIFscore.PhyloPsinglent.mean$CS1.Zm <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS2.Zm <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS3.Zm <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS4.Zm <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS5.Zm <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS6.Zm <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS7.Zm <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS8.Zm <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS9.Zm <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS10.Zm <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS11.Zm <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS12.Zm <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS13.Zm <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS14.Zm <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS15.Zm <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS16.Zm <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$Bivalent.Zm <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$Active.Zm <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$Divergent.Zm <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$Heterochromatin.Zm <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$Quies.Zm <- NA

ZmOs.LECIFscore.PhyloPsinglent.mean$CS1.Os <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS2.Os <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS3.Os <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS4.Os <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS5.Os <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS6.Os <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS7.Os <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS8.Os <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS9.Os <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS10.Os <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS11.Os <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS12.Os <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS13.Os <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS14.Os <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS15.Os <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$CS16.Os <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$Bivalent.Os <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$Active.Os <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$Divergent.Os <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$Heterochromatin.Os <- NA
ZmOs.LECIFscore.PhyloPsinglent.mean$Quies.Os <- NA

## loop for annotating my new merged regions with both scores with previously one with CS presence/absence info

for(i in 1:nrow(ZmOs.LECIFscore.PhyloPsinglent.mean)){
  print(paste0(i/nrow(ZmOs.LECIFscore.PhyloPsinglent.mean)*100, " % covered \n"))
  hit <- which(ZmOs.df.F3$Chr.Zm == ZmOs.LECIFscore.PhyloPsinglent.mean$Chr[i] & ZmOs.df.F3$Start.Zm == ZmOs.LECIFscore.PhyloPsinglent.mean$Start[i] & ZmOs.df.F3$End.Zm == ZmOs.LECIFscore.PhyloPsinglent.mean$End[i])
  if(length(hit) == 1){
    ZmOs.LECIFscore.PhyloPsinglent.mean[i,c(13:54)] <- ZmOs.df.F3[hit,10:51]
  }else if(length(hit) == 0){
    Chr <- which(ZmOs.df.F3$Chr.Zm == ZmOs.LECIFscore.PhyloPsinglent.mean$Chr[i])
    hit2 <- which(ZmOs.df.F3$Start.Zm[Chr] >= ZmOs.LECIFscore.PhyloPsinglent.mean$Start[i] & ZmOs.df.F3$End.Zm[Chr] <= ZmOs.LECIFscore.PhyloPsinglent.mean$End[i])
    if(length(hit2) == 0){stop("Something go wrong boy ... \n")}
    ZmOs.LECIFscore.PhyloPsinglent.mean[i,c(13:54)] <- colSums(ZmOs.df.F3[hit2,10:51])
    colto1 <- c(13:54)[which(ZmOs.LECIFscore.PhyloPsinglent.mean[i,c(13:54)] > 1)]
    ZmOs.LECIFscore.PhyloPsinglent.mean[i, colto1] <- 1
  }else {
    stop("Something go wrong boy ... \n")
  }
}

write.table(ZmOs.LECIFscore.PhyloPsinglent.mean, "H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/ZmOs_CSsimmilarity.DF.LECIFPhyloP.AfterCutOff.txt", col.names = T, row.names = F, sep = "\t", quote = F)

## loop for computing simmilarity/distance between different subsets of LECIF and PhyloP score with different cutoffs (After ...)

CS.Simmilarity.LECIF.PhyloP.ZmOs.After <- data.frame(Bin = rep(c("High,High", "Low,Low", "High,Low", "Low,High"),46),
                                                     Metric = c(rep(c("Jaccard"),4),rep(c("SMC"),4),rep(c("Sokal"),4), rep(c("Tanimoto"),4),
                                                                rep(c("Dice"),4), rep(c("Hanman"),4), rep(c("Ochiai"),4), rep(c("Sokal2"),4),
                                                                rep(c("Phi"),4), rep(c("S2"),4), 
                                                                rep(c("Eucledian"),4), rep(c("Sorensen"),4), rep(c("Canberra"),4), rep(c("Wavehedges"),4), 
                                                                rep(c("Tanimoto2"),4), rep(c("Fidelity"),4), rep(c("SquaredChord"),4), rep(c("SquaredChi"),4), 
                                                                rep(c("AdditiveSym"),4), rep(c("Topsoe"),4), rep(c("Kumarjhonson"),4), rep(c("Manhattan"),4),
                                                                rep(c("Lorentzian"),4), rep(c("Czekanowski"),4), rep(c("Hassebrook"),4),  rep(c("SquaredEucledian"),4),
                                                                rep(c("ProbSymm"),4), rep(c("Kullbackleiber"),4), rep(c("Jensenshannon"),4), rep(c("Minkowski"),4), 
                                                                rep(c("Soergel"),4), rep(c("Innerproduct"),4), rep(c("Hellinger"),4), rep(c("Divergence"),4), 
                                                                rep(c("Jeffreys"),4), rep(c("Jensendifference"),4), rep(c("Chebysev"),4), rep(c("Kulczynskid"),4), 
                                                                rep(c("Harmonicmean"),4), rep(c("Matusita"),4), rep(c("Neyman"),4), rep(c("Clark"),4), 
                                                                rep(c("Kdivergence"),4), rep(c("Taneja"),4), 
                                                                rep(c("Maximum"),4), rep(c("Hamming"),4)),
                                                     CS1 = NA, CS2 = NA, CS3 = NA, CS4 = NA,
                                                     CS5 = NA, CS6 = NA, CS7 = NA, CS8 = NA,
                                                     CS9 = NA, CS10 = NA, CS11 = NA, CS12 = NA,
                                                     CS13 = NA, CS14 = NA, CS15 = NA, CS16 = NA,
                                                     Bivalent = NA, Active = NA, Divergent = NA, Heterochromatin = NA, Quies = NA)

List.CS.Simmilarity.LECIF.PhyloP.ZmOs <- list(After_5050 = CS.Simmilarity.LECIF.PhyloP.ZmOs.After,
                                              After_6040 = CS.Simmilarity.LECIF.PhyloP.ZmOs.After,
                                              After_7030 = CS.Simmilarity.LECIF.PhyloP.ZmOs.After,
                                              After_8020 = CS.Simmilarity.LECIF.PhyloP.ZmOs.After,
                                              After_9010 = CS.Simmilarity.LECIF.PhyloP.ZmOs.After)
bincolumns <- c(8:12)
for(z in 1:length(bincolumns)){ # 4 different cutoffs
  cat(paste0(colnames(ZmOs.LECIFscore.PhyloPsinglent.mean)[bincolumns[z]], " cutoff processing ... \n"))
  for(i in 1:21){# 16 CS and 5 CS functional groups
    cat(paste0(i / 21 * 100, ' % CS groups started \n'))
    for(j in 1:4){ # lecif score percentile bins
      cat(paste0(j / 4 * 100, ' % cutoffs bins started \n'))
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], CS.Os = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i])), method = 1)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+4,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], CS.Os = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i])), method = 2) 
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+8,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], CS.Os = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i])), method = 3)  
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+12,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], CS.Os = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i])), method = 4) 
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+16,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], CS.Os = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i])), method = 5) 
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+20,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], CS.Os = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i])), method = 6) 
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+24,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], CS.Os = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i])), method = 7) 
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+28,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], CS.Os = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i])), method = 8) 
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+32,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], CS.Os = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i])), method = 9) 
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+36,2+i] <- ade4::dist.binary(t(data.frame(CS.Zm = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], CS.Os = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i])), method = 10) 
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+40,2+i] <- philentropy::euclidean(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F) 
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+44,2+i] <- philentropy::sorensen(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+48,2+i] <- philentropy::canberra(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+52,2+i] <- philentropy::wave_hedges(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+56,2+i] <- philentropy::tanimoto(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+60,2+i] <- philentropy::fidelity(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+64,2+i] <- philentropy::squared_chord(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+68,2+i] <- philentropy::squared_chi_sq(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+72,2+i] <- philentropy::additive_symm_chi_sq(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+76,2+i] <- philentropy::topsoe(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F,unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+80,2+i] <- philentropy::kumar_johnson(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F, epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+84,2+i] <- philentropy::manhattan(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+88,2+i] <- philentropy::lorentzian(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F, unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+92,2+i] <- philentropy::czekanowski(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+96,2+i] <- philentropy::kumar_hassebrook(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+100,2+i] <- philentropy::squared_euclidean(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+104,2+i] <- philentropy::prob_symm_chi_sq(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+108,2+i] <- philentropy::kullback_leibler_distance(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F, unit = "log2", epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+112,2+i] <- philentropy::jensen_shannon(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F, unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+116,2+i] <- philentropy::minkowski(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F, n = 3)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+120,2+i] <- philentropy::soergel(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+124,2+i] <- philentropy::inner_product(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+128,2+i] <- philentropy::hellinger(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+132,2+i] <- philentropy::divergence_sq(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+136,2+i] <- philentropy::jeffreys(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F, unit = "log2", epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+140,2+i] <- philentropy::jensen_difference(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F, unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+144,2+i] <- philentropy::chebyshev(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+148,2+i] <- philentropy::kulczynski_d(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F, epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+152,2+i] <- philentropy::harmonic_mean_dist(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+156,2+i] <- philentropy::matusita(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+160,2+i] <- philentropy::neyman_chi_sq(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F, epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+164,2+i] <- philentropy::clark_sq(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+168,2+i] <- philentropy::k_divergence(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F, unit = "log2")
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+172,2+i] <- philentropy::taneja(P = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], Q = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i], testNA = F, unit = "log2", epsilon = 0.01)
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+176,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Zm = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], CS.Os = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i])), method = "maximum", margin = 2)[1,2]
      List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j+180,2+i] <- proxyC::dist(as.matrix(data.frame(CS.Zm = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),12+i], CS.Os = ZmOs.LECIFscore.PhyloPsinglent.mean[which(ZmOs.LECIFscore.PhyloPsinglent.mean[,bincolumns[z]] == List.CS.Simmilarity.LECIF.PhyloP.ZmOs[[z]][j,1]),33+i])), method = "hamming", margin = 2)[1,2]
      cat(paste0(j / 4 * 100, " % cutoff bins completed \n"))
    }
    cat(paste0(i / 21 * 100, " % CS groups completed \n"))
  }
  cat(paste0(colnames(ZmOs.LECIFscore.PhyloPsinglent.mean)[bincolumns[z]], " cutoff computation completed \n"))
}

write.table(List.CS.Simmilarity.LECIF.PhyloP.ZmOs$After_5050, file = "H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/ZmOs_CSsimmilarity.DF.Distances.After5050.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(List.CS.Simmilarity.LECIF.PhyloP.ZmOs$After_6040, file = "H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/ZmOs_CSsimmilarity.DF.Distances.After6040.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(List.CS.Simmilarity.LECIF.PhyloP.ZmOs$After_7030, file = "H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/ZmOs_CSsimmilarity.DF.Distances.After7030.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(List.CS.Simmilarity.LECIF.PhyloP.ZmOs$After_8020, file = "H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/ZmOs_CSsimmilarity.DF.Distances.After8020.txt", col.names = T, row.names = F, sep = "\t", quote = F)
write.table(List.CS.Simmilarity.LECIF.PhyloP.ZmOs$After_9010, file = "H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/ZmOs_CSsimmilarity.DF.Distances.After9010.txt", col.names = T, row.names = F, sep = "\t", quote = F)


### 1.4 Correlation of LECIF score and PhyloP score (include only the numbers in the figure)

### 1.4.1 Arabidopsis-Oryza

cor(x = AtOs.LECIFscore.PhyloPsinglent.mean$LECIF, y = AtOs.LECIFscore.PhyloPsinglent.mean$PhyloP, method = "pearson") # 0.004074123
cor(x = AtOs.LECIFscore.PhyloPsinglent.mean$LECIF, y = AtOs.LECIFscore.PhyloPsinglent.mean$PhyloP, method = "spearman") # 0.00472545

### 1.4.2 Arabidopsis-Zea

cor(x = AtZm.LECIFscore.PhyloPsinglent.mean$LECIF, y = AtZm.LECIFscore.PhyloPsinglent.mean$PhyloP, method = "pearson") # 0.02191947
cor(x = AtZm.LECIFscore.PhyloPsinglent.mean$LECIF, y = AtZm.LECIFscore.PhyloPsinglent.mean$PhyloP, method = "spearman") # 0.0146263

### 1.4.3 Oryza-Arabidopsis

cor(x = OsAt.LECIFscore.PhyloPsinglent.mean$LECIF, y = OsAt.LECIFscore.PhyloPsinglent.mean$PhyloP, method = "pearson") # 0.042502777
cor(x = OsAt.LECIFscore.PhyloPsinglent.mean$LECIF, y = OsAt.LECIFscore.PhyloPsinglent.mean$PhyloP, method = "spearman") # 0.03056301

### 1.4.4 Oryza-Zea

cor(x = OsZm.LECIFscore.PhyloPsinglent.mean$LECIF, y = OsZm.LECIFscore.PhyloPsinglent.mean$PhyloP, method = "pearson") # 0.1197667
cor(x = OsZm.LECIFscore.PhyloPsinglent.mean$LECIF, y = OsZm.LECIFscore.PhyloPsinglent.mean$PhyloP, method = "spearman") # 0.1182813

### 1.4.5 Zea-Arabidopsis

cor(x = ZmAt.LECIFscore.PhyloPsinglent.mean$LECIF, y = ZmAt.LECIFscore.PhyloPsinglent.mean$PhyloP, method = "pearson") # 0.006604811
cor(x = ZmAt.LECIFscore.PhyloPsinglent.mean$LECIF, y = ZmAt.LECIFscore.PhyloPsinglent.mean$PhyloP, method = "spearman") # 0.007103574

### 1.4.6 Zea-Oryza

cor(x = ZmOs.LECIFscore.PhyloPsinglent.mean$LECIF, y = ZmOs.LECIFscore.PhyloPsinglent.mean$PhyloP, method = "pearson") # 0.02857325
cor(x = ZmOs.LECIFscore.PhyloPsinglent.mean$LECIF, y = ZmOs.LECIFscore.PhyloPsinglent.mean$PhyloP, method = "spearman") # 0.02857325

###############################################
######### 2. LECIF score Conservation ######### 
##############################################

### 2. LECIF score over conservation elements from different multiple genome alignments for each species (PhastCons and CNEs)

### 2.1 Arabidopsis

AtOs_regions <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/At-Os_LECIFscore.bedgraph", header = F, sep = "\t")
AtOs_CNEs_Score <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/CNEs_At-Os_LECIF.unique.sorted.bed", header = F, sep = "\t")
AtOs_CEs_Score <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/CE_At-Os_LECIF.unique.sorted.bed", header = F, sep = "\t")
colnames(AtOs_CNEs_Score) <- c("Chr", "Start", "End", "Score")
colnames(AtOs_CEs_Score) <- c("Chr", "Start", "End", "Score")
AtOs_CNEs_Score$Class <- rep("CNEs", nrow(AtOs_CNEs_Score))
AtOs_CEs_Score$Class <- rep("PhastCons", nrow(AtOs_CEs_Score))
AtOs_Conservation_Coverage <- data.frame(Class = c("CNEs", "PhastCons"), Coverage = c((sum(AtOs_CNEs_Score$End - AtOs_CNEs_Score$Start)/sum(AtOs_regions$V3 - AtOs_regions$V2))*100,
                                                                                      (sum(AtOs_CEs_Score$End - AtOs_CEs_Score$Start)/sum(AtOs_regions$V3 - AtOs_regions$V2))*100),
                                         CoverageRounder = c(round((sum(AtOs_CNEs_Score$End - AtOs_CNEs_Score$Start)/sum(AtOs_regions$V3 - AtOs_regions$V2))*100,2),
                                                             round((sum(AtOs_CEs_Score$End - AtOs_CEs_Score$Start)/sum(AtOs_regions$V3 - AtOs_regions$V2))*100,2)))
AtOs_Conservation_Scores <- rbind(AtOs_CNEs_Score, AtOs_CEs_Score)

AtZm_regions <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/At-Zm_LECIFscore.bedgraph", header = F, sep = "\t")
AtZm_CNEs_Score <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/CNEs_At-Zm_LECIF.unique.sorted.bed", header = F, sep = "\t")
AtZm_CEs_Score  <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/CE_At-Zm_LECIF.unique.sorted.bed", header = F, sep = "\t")
colnames(AtZm_CNEs_Score) <- c("Chr", "Start", "End", "Score")
colnames(AtZm_CEs_Score) <- c("Chr", "Start", "End", "Score")
AtZm_CNEs_Score$Class <- rep("CNEs", nrow(AtZm_CNEs_Score))
AtZm_CEs_Score$Class <- rep("PhastCons", nrow(AtZm_CEs_Score))
AtZm_Conservation_Coverage <- data.frame(Class = c("CNEs", "PhastCons"), Coverage = c((sum(AtZm_CNEs_Score$End - AtZm_CNEs_Score$Start)/sum(AtZm_regions$V3 - AtZm_regions$V2))*100,
                                                                                      (sum(AtZm_CEs_Score$End - AtZm_CEs_Score$Start)/sum(AtZm_regions$V3 - AtZm_regions$V2))*100),
                                         CoverageRounder = c(round((sum(AtZm_CNEs_Score$End - AtZm_CNEs_Score$Start)/sum(AtZm_regions$V3 - AtZm_regions$V2))*100,2),
                                                             round((sum(AtZm_CEs_Score$End - AtZm_CEs_Score$Start)/sum(AtZm_regions$V3 - AtZm_regions$V2))*100,2)))
AtZm_Conservation_Scores <- rbind(AtZm_CNEs_Score, AtZm_CEs_Score)

### 2.2 Oryza

OsAt_regions <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/Os-At_LECIFscore.bedgraph", header = F, sep = "\t")
OsAt_CNEs_Score <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/CNEs_Os-At_LECIF.unique.sorted.bed", header = F, sep = "\t")
OsAt_CEs_Score <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/CE_Os-At_LECIF.unique.sorted.bed", header = F, sep = "\t")
colnames(OsAt_CNEs_Score) <- c("Chr", "Start", "End", "Score")
colnames(OsAt_CEs_Score) <- c("Chr", "Start", "End", "Score")
OsAt_CNEs_Score$Class <- rep("CNEs", nrow(OsAt_CNEs_Score))
OsAt_CEs_Score$Class <- rep("PhastCons", nrow(OsAt_CEs_Score))
OsAt_Conservation_Coverage <- data.frame(Class = c("CNEs", "PhastCons"), Coverage = c((sum(OsAt_CNEs_Score$End - OsAt_CNEs_Score$Start)/sum(OsAt_regions$V3 - OsAt_regions$V2))*100,
                                                                                      (sum(OsAt_CEs_Score$End - OsAt_CEs_Score$Start)/sum(OsAt_regions$V3 - OsAt_regions$V2))*100),
                                         CoverageRounder = c(round((sum(OsAt_CNEs_Score$End - OsAt_CNEs_Score$Start)/sum(OsAt_regions$V3 - OsAt_regions$V2))*100,2),
                                                             round((sum(OsAt_CEs_Score$End - OsAt_CEs_Score$Start)/sum(OsAt_regions$V3 - OsAt_regions$V2))*100,2)))
OsAt_Conservation_Scores <- rbind(OsAt_CNEs_Score, OsAt_CEs_Score)

OsZm_regions <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/Os-Zm_LECIFscore.bedgraph", header = F, sep = "\t")
OsZm_CNEs_Score <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/CNEs_Os-Zm_LECIF.unique.sorted.bed", header = F, sep = "\t")
OsZm_CEs_Score  <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/CE_Os-Zm_LECIF.unique.sorted.bed", header = F, sep = "\t")
colnames(OsZm_CNEs_Score) <- c("Chr", "Start", "End", "Score")
colnames(OsZm_CEs_Score) <- c("Chr", "Start", "End", "Score")
OsZm_CNEs_Score$Class <- rep("CNEs", nrow(OsZm_CNEs_Score))
OsZm_CEs_Score$Class <- rep("PhastCons", nrow(OsZm_CEs_Score))
OsZm_Conservation_Coverage <- data.frame(Class = c("CNEs", "PhastCons"), Coverage = c((sum(OsZm_CNEs_Score$End - OsZm_CNEs_Score$Start)/sum(OsZm_regions$V3 - OsZm_regions$V2))*100,
                                                                                      (sum(OsZm_CEs_Score$End - OsZm_CEs_Score$Start)/sum(OsZm_regions$V3 - OsZm_regions$V2))*100),
                                         CoverageRounder = c(round((sum(OsZm_CNEs_Score$End - OsZm_CNEs_Score$Start)/sum(OsZm_regions$V3 - OsZm_regions$V2))*100,2),
                                                             round((sum(OsZm_CEs_Score$End - OsZm_CEs_Score$Start)/sum(OsZm_regions$V3 - OsZm_regions$V2))*100,2)))
OsZm_Conservation_Scores <- rbind(OsZm_CNEs_Score, OsZm_CEs_Score)

### 2.3 Zea

ZmAt_regions <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/Zm-At_LECIFscore.bedgraph", header = F, sep = "\t")
ZmAt_CNEs_Score <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/CNEs_Zm-At_LECIF.unique.sorted.bed", header = F, sep = "\t")
ZmAt_CEs_Score <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/CE_Zm-At_LECIF.unique.sorted.bed", header = F, sep = "\t")
colnames(ZmAt_CNEs_Score) <- c("Chr", "Start", "End", "Score")
colnames(ZmAt_CEs_Score) <- c("Chr", "Start", "End", "Score")
ZmAt_CNEs_Score$Class <- rep("CNEs", nrow(ZmAt_CNEs_Score))
ZmAt_CEs_Score$Class <- rep("PhastCons", nrow(ZmAt_CEs_Score))
ZmAt_Conservation_Coverage <- data.frame(Class = c("CNEs", "PhastCons"), Coverage = c((sum(ZmAt_CNEs_Score$End - ZmAt_CNEs_Score$Start)/sum(ZmAt_regions$V3 - ZmAt_regions$V2))*100,
                                                                                      (sum(ZmAt_CEs_Score$End - ZmAt_CEs_Score$Start)/sum(ZmAt_regions$V3 - ZmAt_regions$V2))*100),
                                         CoverageRounder = c(round((sum(ZmAt_CNEs_Score$End - ZmAt_CNEs_Score$Start)/sum(ZmAt_regions$V3 - ZmAt_regions$V2))*100,2),
                                                             round((sum(ZmAt_CEs_Score$End - ZmAt_CEs_Score$Start)/sum(ZmAt_regions$V3 - ZmAt_regions$V2))*100,2)))
ZmAt_Conservation_Scores <- rbind(ZmAt_CNEs_Score, ZmAt_CEs_Score)

ZmOs_regions <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/Zm-Os_LECIFscore.bedgraph", header = F, sep = "\t")
ZmOs_CNEs_Score <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/CNEs_Zm-Os_LECIF.unique.sorted.bed", header = F, sep = "\t")
ZmOs_CEs_Score <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/CE_Zm-Os_LECIF.unique.sorted.bed", header = F, sep = "\t")
colnames(ZmOs_CNEs_Score) <- c("Chr", "Start", "End", "Score")
colnames(ZmOs_CEs_Score) <- c("Chr", "Start", "End", "Score")
ZmOs_CNEs_Score$Class <- rep("CNEs", nrow(ZmOs_CNEs_Score))
ZmOs_CEs_Score$Class <- rep("PhastCons", nrow(ZmOs_CEs_Score))
ZmOs_Conservation_Coverage <- data.frame(Class = c("CNEs", "PhastCons"), Coverage = c((sum(ZmOs_CNEs_Score$End - ZmOs_CNEs_Score$Start)/sum(ZmOs_regions$V3 - ZmOs_regions$V2))*100,
                                                                                      (sum(ZmOs_CEs_Score$End - ZmOs_CEs_Score$Start)/sum(ZmOs_regions$V3 - ZmOs_regions$V2))*100),
                                         CoverageRounder = c(round((sum(ZmOs_CNEs_Score$End - ZmOs_CNEs_Score$Start)/sum(ZmOs_regions$V3 - ZmOs_regions$V2))*100,2),
                                                             round((sum(ZmOs_CEs_Score$End - ZmOs_CEs_Score$Start)/sum(ZmOs_regions$V3 - ZmOs_regions$V2))*100,2)))
ZmOs_Conservation_Scores <- rbind(ZmOs_CNEs_Score, ZmOs_CEs_Score)

###############################################
######### 3. LECIF score Variability ######### 
##############################################

### 3. Check if GWAS SNPs are overlap enriched in different LECIF bins than common SNPs

### 3.1 Arabidopsis

AtOs.df.F3 <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/AtOs_CSsimmilarity.DF3.txt", header = T, sep = "\t", row.names = F)
for(i in c("First", "Second", "Thirdth", "Fourth", "Fifth")){
  write.table(AtOs.df.F3[which(AtOs.df.F3$Percentile.Rank.Bin.Five == i),c(1,2,3)], file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/Annotation/0.LOLA/regionDB/At/LECIF_Os/LECIF_", i, ".bed"), col.names = F, row.names = F, sep = "\t", quote = F)
}

regionDB <- loadRegionDB(dbLocation = "H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/Annotation/0.LOLA/regionDB/At/", collections = "LECIF_Os")

Universe <- read.table("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/SNPs_Common_At-Os_LECIF.unique.sorted.bed")
Universe <- Universe[,c(1,2,3)]
colnames(Universe) <- c("chr", "start", "end")
Universe <- makeGRangesFromDataFrame(Universe)

UserSet <- read.table("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/SNPs_GWAS_At-Os_LECIF.unique.sorted.bed")
UserSet <- UserSet[,c(1,2,3)]
colnames(UserSet) <- c("chr", "start", "end")
UserSet <- makeGRangesFromDataFrame(UserSet)

RegionResults_AtOs <- runLOLA(UserSet, Universe, regionDB, cores = 16)

AtZm.df.F3 <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/AtZm_CSsimmilarity.DF3.txt", header = T, sep = "\t", row.names = F)
for(i in c("First", "Second", "Thirdth", "Fourth", "Fifth")){
  write.table(AtZm.df.F3[which(AtZm.df.F3$Percentile.Rank.Bin.Five == i),c(1,2,3)], file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/Annotation/0.LOLA/regionDB/At/LECIF_Zm/regions/LECIF_", i, ".bed"), col.names = F, row.names = F, sep = "\t", quote = F)
}

regionDB <- loadRegionDB(dbLocation = "H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/Annotation/0.LOLA/regionDB/At/", collections = "LECIF_Zm")

Universe <- read.table("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/SNPs_Common_At-Zm_LECIF.unique.sorted.bed")
Universe <- Universe[,c(1,2,3)]
colnames(Universe) <- c("chr", "start", "end")
Universe <- makeGRangesFromDataFrame(Universe)

UserSet <- read.table("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/SNPs_GWAS_At-Zm_LECIF.unique.sorted.bed")
UserSet <- UserSet[,c(1,2,3)]
colnames(UserSet) <- c("chr", "start", "end")
UserSet <- makeGRangesFromDataFrame(UserSet)

RegionResults_AtZm <- runLOLA(UserSet, Universe, regionDB, cores = 16)

### 3.2 Oryza

OsAt.df.F3 <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/OsAt_CSsimmilarity.DF3.txt", header = T, sep = "\t")
for(i in c("First", "Second", "Thirdth", "Fourth", "Fifth")){
  write.table(OsAt.df.F3[which(OsAt.df.F3$Percentile.Rank.Bin.Five == i),c(1,2,3)], file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/Annotation/0.LOLA/regionDB/Os/LECIF_At/LECIF_", i, ".bed"), col.names = F, row.names = F, sep = "\t", quote = F)
}

regionDB <- loadRegionDB(dbLocation = "H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/Annotation/0.LOLA/regionDB/Os/", collections = "LECIF_At")

Universe <- read.table("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/SNPs_Common_Os-At_LECIF.unique.sorted.bed")
Universe <- Universe[,c(1,2,3)]
colnames(Universe) <- c("chr", "start", "end")
Universe <- makeGRangesFromDataFrame(Universe)

UserSet <- read.table("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/SNPs_GWAS_Os-At_LECIF.unique.sorted.bed")
UserSet <- UserSet[,c(1,2,3)]
colnames(UserSet) <- c("chr", "start", "end")
UserSet <- makeGRangesFromDataFrame(UserSet)

RegionResults_OsAt <- runLOLA(UserSet, Universe, regionDB, cores = 16)

OsZm.df.F3 <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/OsZm_CSsimmilarity.DF3.txt", header = T, sep = "\t")
for(i in c("First", "Second", "Thirdth", "Fourth", "Fifth")){
  write.table(OsZm.df.F3[which(OsZm.df.F3$Percentile.Rank.Bin.Five == i),c(1,2,3)], file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/Annotation/0.LOLA/regionDB/Os/LECIF_Zm/regions/LECIF_", i, ".bed"), col.names = F, row.names = F, sep = "\t", quote = F)
}

regionDB <- loadRegionDB(dbLocation = "H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/Annotation/0.LOLA/regionDB/Os/", collections = "LECIF_Zm")

Universe <- read.table("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/SNPs_Common_Os-Zm_LECIF.unique.sorted.bed")
Universe <- Universe[,c(1,2,3)]
colnames(Universe) <- c("chr", "start", "end")
Universe <- makeGRangesFromDataFrame(Universe)

UserSet <- read.table("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/SNPs_GWAS_Os-Zm_LECIF.unique.sorted.bed")
UserSet <- UserSet[,c(1,2,3)]
colnames(UserSet) <- c("chr", "start", "end")
UserSet <- makeGRangesFromDataFrame(UserSet)

RegionResults_OsZm <- runLOLA(UserSet, Universe, regionDB, cores = 16)

### 3.3 Zea

ZmAt.df.F3 <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/ZmAt_CSsimmilarity.DF3.txt", header = T, sep = "\t")
for(i in c("First", "Second", "Thirdth", "Fourth", "Fifth")){
  write.table(ZmAt.df.F3[which(ZmAt.df.F3$Percentile.Rank.Bin.Five == i),c(1,2,3)], file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/Annotation/0.LOLA/regionDB/Zm/LECIF_At/LECIF_", i, ".bed"), col.names = F, row.names = F, sep = "\t", quote = F)
}

regionDB <- loadRegionDB(dbLocation = "H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/Annotation/0.LOLA/regionDB/Zm/", collections = "LECIF_At")

Universe <- read.table("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/SNPs_Common_Zm-At_LECIF.unique.sorted.bed")
Universe <- Universe[,c(1,2,3)]
colnames(Universe) <- c("chr", "start", "end")
Universe <- makeGRangesFromDataFrame(Universe)

UserSet <- read.table("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/SNPs_GWAS_Zm-At_LECIF.unique.sorted.bed")
UserSet <- UserSet[,c(1,2,3)]
colnames(UserSet) <- c("chr", "start", "end")
UserSet <- makeGRangesFromDataFrame(UserSet)

RegionResults_ZmAt <- runLOLA(UserSet, Universe, regionDB, cores = 16)

ZmOs.df.F3 <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/ZmOs_CSsimmilarity.DF3.txt", header = T, sep = "\t")
for(i in c("First", "Second", "Thirdth", "Fourth", "Fifth")){
  write.table(ZmOs.df.F3[which(ZmOs.df.F3$Percentile.Rank.Bin.Five == i),c(1,2,3)], file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/Annotation/0.LOLA/regionDB/Zm/LECIF_Os/regions/LECIF_", i, ".bed"), col.names = F, row.names = F, sep = "\t", quote = F)
}

regionDB <- loadRegionDB(dbLocation = "H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/Annotation/0.LOLA/regionDB/Zm/", collections = "LECIF_Os")

Universe <- read.table("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/SNPs_Common_Zm-Os_LECIF.unique.sorted.bed")
Universe <- Universe[,c(1,2,3)]
colnames(Universe) <- c("chr", "start", "end")
Universe <- makeGRangesFromDataFrame(Universe)

UserSet <- read.table("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/SNPs_GWAS_Zm-Os_LECIF.unique.sorted.bed")
UserSet <- UserSet[,c(1,2,3)]
colnames(UserSet) <- c("chr", "start", "end")
UserSet <- makeGRangesFromDataFrame(UserSet)

RegionResults_ZmOs <- runLOLA(UserSet, Universe, regionDB, cores = 16)

###################################################################
#################### 4. Plotting (Fusion) zone #################### 
###################################################################

### 4.1.1 Arabisopsis-Oryza Chromatin States-based

## Violin-plot CS score distribution

data <- data.frame(name = LECIF_AtOs$V5, value = LECIF_AtOs$V4)
data$names <- factor(data$name , levels=c("CS16", "CS15", "CS14", "CS13",
                                          "CS12", "CS11", "CS10", "CS9",
                                          "CS8", "CS7", "CS6", "CS5",
                                          "CS4", "CS3", "CS2", "CS1"), ordered = T)
r1<-with(data, tapply(value, name, mean))
r2<-with(data, tapply(value, name, median))
g <-
  ggplot(data, aes(x = names, y = value,
                   color = value)) +
  labs(x = "Chromatin states", y = "LECIF score") +
  scale_color_brewer(palette = "Dark2", guide = "none")
A <- g + geom_violin(aes(fill = names), size = 1, alpha = .5) +
  geom_boxplot(outlier.alpha = 0, coef = 0, color = "black", width = .2) +
  scale_fill_manual(values=c(rep("gray70",3),  rep("gold",3), rep("#4DBBD5FF",1), rep("#00A087FF",4), rep("#3C5488FF",5))) +
  coord_flip() + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  ylim(0,1)

## LECIF score distribution

B <- ggplot(data, aes(x=value)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="gray70") +
  geom_density(alpha=.4, fill="black") + theme_bw() + ggtitle("LECIF score distribution") + 
  labs(y = "Proportion") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) + 
  xlim(0,1) +
  annotate("text", x = 0.80, y = 7, label=paste0("Mean: ", round(mean(data$value), digits = 2)), size = 5, color = "black") +
  annotate("text", x = 0.82, y = 6, label=paste0("Median: ", round(median(data$value), digits = 2)), size = 5, color = "black")

## Coverage of each CS in across aligning regions (first species is the target genome)

aligningRegions <- sum(LECIF_AtOs$V3 - LECIF_AtOs$V2)
Coverage <- data.frame("Chromatin States" = c("CS1","CS2", "CS3", "CS4", "CS5", "CS6", "CS7", "CS8", "CS9", "CS10", "CS11", "CS12", "CS13", "CS14", "CS15", "CS16"),
                       Coverage = NA)
for(i in 1:nrow(Coverage)){
  Coverage[i,2] <- (sum(LECIF_AtOs$V3[which(LECIF_AtOs$V5 == Coverage[i,1])] - LECIF_AtOs$V2[which(LECIF_AtOs$V5 == Coverage[i,1])])/aligningRegions)*100 
}
Coverage$Chromatin.States <- factor(Coverage$Chromatin.States , levels=c("CS16", "CS15", "CS14", "CS13",
                                                                         "CS12", "CS11", "CS10", "CS9",
                                                                         "CS8", "CS7", "CS6", "CS5",
                                                                         "CS4", "CS3", "CS2", "CS1"), ordered = TRUE) 
Coverage$CoverageRounded <- round(Coverage$Coverage, digits = 2)
C <- ggplot(Coverage, aes(x = Chromatin.States, y = Coverage))+
  geom_col(width = 0.7) + coord_flip() + 
  theme_bw() + ylab("Coverage (%)") + xlab("Chromatin states") + ylim(0,100) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(data=Coverage,aes(x=Chromatin.States,y=Coverage,label=CoverageRounded),vjust=0.5, hjust = -0.2) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
 ### ADD TEXT CORRELATION


## Lineplot Chromatin State Functional group Simmilairty in three LECIF percentile bins

Lineplot.table.AtOs <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/AtOs_CSsimmilarity.DF.Distances3.txt", header = T, sep = "\t")
Lineplot.table.AtOs <- Lineplot.table.AtOs[c(13,14,15),]
Lineplot.table.AtOs <- Lineplot.table.AtOs[,c(19,20,21,22,23)]
Lineplot.df.AtOs <- data.frame(CSGroup = c(rep("Bivalent",3), rep("Active",3), rep("Divergent",3), rep("Heterochromatin",3), rep("Quies",3)),
                                  LECIF.Bin = c(rep(c(0.5,1.5,2.5),5)),
                                  Simmilarity = c(Lineplot.table.AtOs$Bivalent, Lineplot.table.AtOs$Active, Lineplot.table.AtOs$Divergent,
                                                  Lineplot.table.AtOs$Heterochromatin, Lineplot.table.AtOs$Quies))
Lineplot.df.AtOs$Simmilarity <- 1 - Lineplot.df.AtOs$Simmilarity
Lineplot.df.AtOs$LECIF.Bin <- as.numeric(as.vector(Lineplot.df.AtOs$LECIF.Bin))
D <- ggplot(data=Lineplot.df.AtOs, aes(x=LECIF.Bin, y=Simmilarity, group=CSGroup, color=CSGroup)) +
  geom_line(size=1.0) + geom_point(size=3)+
  scale_color_manual(values=c("#00A087FF", "#3C5488FF", "#4DBBD5FF", "gold", "gray70"))+
  theme_bw() + coord_flip() + theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlim(0.5,2.5) + ylim(0,0.5) + ggtitle("Simmilarity in 3 LECIF bins")


## Grouped-barplot (two bars) Chromatin  state Simmilarity between high;low and low;high LECIF;PhyloP scores (50;50 | 60;40 cutoffs)

Barplot.table.AtOs <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/AtOs_CSsimmilarity.DF.Distances.After6040.txt", header = T, sep = "\t")
Barplot.table.AtOs <- Barplot.table.AtOs[c(19,20),]
Barplot.table.AtOs <- Barplot.table.AtOs[,c(1:18)]
Barplot.df.AtOs <- data.frame(CS = c(rep("CS1",2),rep("CS2",2),rep("CS3",2),rep("CS4",2),
                                     rep("CS5",2),rep("CS6",2),rep("CS7",2),rep("CS8",2),
                                     rep("CS9",2),rep("CS10",2),rep("CS11",2),rep("CS12",2),
                                     rep("CS13",2),rep("CS14",2),rep("CS15",2),rep("CS16",2)), Class = c(rep(c("High,Low", "Low,High"),16)), Simmilarity = c(Barplot.table.AtOs$CS1,Barplot.table.AtOs$CS2,Barplot.table.AtOs$CS3,Barplot.table.AtOs$CS4,
                                                                                                                       Barplot.table.AtOs$CS5,Barplot.table.AtOs$CS6,Barplot.table.AtOs$CS7,Barplot.table.AtOs$CS8,
                                                                                                                       Barplot.table.AtOs$CS9,Barplot.table.AtOs$CS10,Barplot.table.AtOs$CS11,Barplot.table.AtOs$CS12,
                                                                                                                       Barplot.table.AtOs$CS13,Barplot.table.AtOs$CS14,Barplot.table.AtOs$CS15,Barplot.table.AtOs$CS16))
Barplot.df.AtOs$Simmilarity <- 1 - Barplot.df.AtOs$Simmilarity
Barplot.df.AtOs$CS<-factor(Barplot.df.AtOs$CS, levels=c("CS16", "CS15", "CS14", "CS13",
                                        "CS12", "CS11", "CS10", "CS9",
                                        "CS8", "CS7", "CS6", "CS5",
                                        "CS4", "CS3", "CS2", "CS1"), ordered = T)
Barplot.df.AtOs$Simmilarity <- 1 - Barplot.df.AtOs$Simmilarity
Barplot.df.AtOs$Class <- as.factor(Barplot.df.AtOs$Class)
# Grouped
E <- ggplot(Barplot.df.AtOs, aes(fill=Class, y=Simmilarity, x=CS)) + 
  geom_bar(position="dodge", stat="identity") + coord_flip() + theme_bw() + ylim(0,0.5) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")


### Arabidopsis-Oryza Conservation Boxplot CNEs|PhastCons Score distribution and coverage

eFe <- ggplot(AtOs_Conservation_Scores, aes(x = Class, y = Score)) +
  geom_boxplot(lwd=2) + theme_bw() + coord_flip() + ylim(0,1) + 
  theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none") +
  labs( y = "LECIF score") + geom_hline(aes(yintercept=0.16), colour="darkgray") + geom_hline(aes(yintercept=0.16), colour="darkgray")


H <- ggplot(AtOs_Conservation_Coverage, aes(x = Class, y = Coverage))+
  geom_col(width = 0.7) + coord_flip() + 
  theme_bw() + ylab("Coverage (%)") + xlab("Conserved Elements") + ylim(0,100) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(data=AtOs_Conservation_Coverage,aes(x=Class,y=Coverage,label=CoverageRounder),vjust=0.5, hjust = -0.2) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())


### 4.1.2 Arabidopsis-Zea Chromatin States-based (add correlation text)

## Violin-plot CS score distribution

data2 <- data.frame(name = LECIF_AtZm$V5, value = LECIF_AtZm$V4)
data2$names <- factor(data2$name , levels=c("CS16", "CS15", "CS14", "CS13",
                                            "CS12", "CS11", "CS10", "CS9",
                                            "CS8", "CS7", "CS6", "CS5",
                                            "CS4", "CS3", "CS2", "CS1"), ordered = T)
r3<-with(data2, tapply(value, name, mean))
r4<-with(data2, tapply(value, name, median))
g <-
  ggplot(data2, aes(x = names, y = value,
                    color = value)) +
  labs(x = "Chromatin states", y = "LECIF score") +
  scale_color_brewer(palette = "Dark2", guide = "none")
A2 <- g + geom_violin(aes(fill = names), size = 1, alpha = .5) +
  geom_boxplot(outlier.alpha = 0, coef = 0, color = "black", width = .2) +
  scale_fill_manual(values=c(rep("gray70",3),  rep("gold",3), rep("#4DBBD5FF",1), rep("#00A087FF",4), rep("#3C5488FF",5))) +
  coord_flip() + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  ylim(0,1)

## LECIF score distribution

B2 <- ggplot(data2, aes(x=value)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="gray70") +
  geom_density(alpha=.4, fill="black") + theme_bw() + ggtitle("LECIF score distribution") + 
  labs(y = "Proportion") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) + 
  xlim(0,1) +
  annotate("text", x = 0.80, y = 5.2, label=paste0("Mean: ", round(mean(data2$value), digits = 2)), size = 5, color = "black") +
  annotate("text", x = 0.82, y = 4.2, label=paste0("Median: ", round(median(data2$value), digits = 2)), size = 5, color = "black")

## Coverage of each CS in across aligning regions (first species is the target genome)

aligningRegions2 <- sum(LECIF_AtZm$V3 - LECIF_AtZm$V2)
Coverage2 <- data.frame("Chromatin States" = c("CS1","CS2", "CS3", "CS4", "CS5", "CS6", "CS7", "CS8", "CS9", "CS10", "CS11", "CS12", "CS13", "CS14", "CS15", "CS16"),
                        Coverage = NA)
for(i in 1:nrow(Coverage2)){
  Coverage2[i,2] <- (sum(LECIF_AtZm$V3[which(LECIF_AtZm$V5 == Coverage2[i,1])] - LECIF_AtZm$V2[which(LECIF_AtZm$V5 == Coverage2[i,1])])/aligningRegions2)*100 
}
Coverage2$Chromatin.States <- factor(Coverage2$Chromatin.States , levels=c("CS16", "CS15", "CS14", "CS13",
                                                                           "CS12", "CS11", "CS10", "CS9",
                                                                           "CS8", "CS7", "CS6", "CS5",
                                                                           "CS4", "CS3", "CS2", "CS1"), ordered = TRUE) 
Coverage2$CoverageRounded <- round(Coverage2$Coverage, digits = 2)
C2 <- ggplot(Coverage2, aes(x = Chromatin.States, y = Coverage))+
  geom_col(width = 0.7) + coord_flip() + 
  theme_bw() + ylab("Coverage (%)") + xlab("Chromatin states") + ylim(0,100) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(data=Coverage2,aes(x=Chromatin.States,y=Coverage,label=CoverageRounded),vjust=0.5, hjust = -0.2) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
 ### ADD TEXT CORRELATION


## Lineplot Chromatin State Functional group Simmilairty in three LECIF percentile bins

Lineplot.table.AtZm <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/AtZm_CSsimmilarity.DF.Distances3.txt", header = T, sep = "\t")
Lineplot.table.AtZm <- Lineplot.table.AtZm[c(13,14,15),]
Lineplot.table.AtZm <- Lineplot.table.AtZm[,c(19,20,21,22,23)]
Lineplot.df.AtZm <- data.frame(CSGroup = c(rep("Bivalent",3), rep("Active",3), rep("Divergent",3), rep("Heterochromatin",3), rep("Quies",3)),
                               LECIF.Bin = c(rep(c(0.5,1.5,2.5),5)),
                               Simmilarity = c(Lineplot.table.AtZm$Bivalent, Lineplot.table.AtZm$Active, Lineplot.table.AtZm$Divergent,
                                               Lineplot.table.AtZm$Heterochromatin, Lineplot.table.AtZm$Quies))
Lineplot.df.AtZm$Simmilarity <- 1 - Lineplot.df.AtZm$Simmilarity
Lineplot.df.AtZm$LECIF.Bin <- as.numeric(as.vector(Lineplot.df.AtZm$LECIF.Bin))
D2 <- ggplot(data=Lineplot.df.AtZm, aes(x=LECIF.Bin, y=Simmilarity, group=CSGroup, color=CSGroup)) +
  geom_line(size=1.0) + geom_point(size=3)+
  scale_color_manual(values=c("#00A087FF", "#3C5488FF", "#4DBBD5FF", "gold", "gray70"))+
  theme_bw() + coord_flip() + theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlim(0.5,2.5) + ylim(0,0.5) + ggtitle("Simmilarity in 3 LECIF bins")


## Grouped-barplot (two bars) Chromatin  state Simmilarity between high;low and low;high LECIF;PhyloP scores (50;50 | 60;40 cutoffs)

Barplot.table.AtZm <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/AtZm_CSsimmilarity.DF.Distances.After6040.txt", header = T, sep = "\t")
Barplot.table.AtZm <- Barplot.table.AtZm[c(19,20),]
Barplot.table.AtZm <- Barplot.table.AtZm[,c(1:18)]
Barplot.df.AtZm <- data.frame(CS = c(rep("CS1",2),rep("CS2",2),rep("CS3",2),rep("CS4",2),
                                     rep("CS5",2),rep("CS6",2),rep("CS7",2),rep("CS8",2),
                                     rep("CS9",2),rep("CS10",2),rep("CS11",2),rep("CS12",2),
                                     rep("CS13",2),rep("CS14",2),rep("CS15",2),rep("CS16",2)), Class = c(rep(c("High,Low", "Low,High"),16)), Simmilarity = c(Barplot.table.AtZm$CS1,Barplot.table.AtZm$CS2,Barplot.table.AtZm$CS3,Barplot.table.AtZm$CS4,
                                                                                                                                                             Barplot.table.AtZm$CS5,Barplot.table.AtZm$CS6,Barplot.table.AtZm$CS7,Barplot.table.AtZm$CS8,
                                                                                                                                                             Barplot.table.AtZm$CS9,Barplot.table.AtZm$CS10,Barplot.table.AtZm$CS11,Barplot.table.AtZm$CS12,
                                                                                                                                                             Barplot.table.AtZm$CS13,Barplot.table.AtZm$CS14,Barplot.table.AtZm$CS15,Barplot.table.AtZm$CS16))
Barplot.df.AtZm$Simmilarity <- 1 - Barplot.df.AtZm$Simmilarity
Barplot.df.AtZm$CS<-factor(Barplot.df.AtZm$CS, levels=c("CS16", "CS15", "CS14", "CS13",
                                                        "CS12", "CS11", "CS10", "CS9",
                                                        "CS8", "CS7", "CS6", "CS5",
                                                        "CS4", "CS3", "CS2", "CS1"), ordered = T)
# Grouped
E2 <- ggplot(Barplot.df.AtZm, aes(fill=Class, y=Simmilarity, x=CS)) + 
  geom_bar(position="dodge", stat="identity") + coord_flip() + theme_bw() + ylim(0,0.5) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")

### Arabidopsis-Oryza Conservation Boxplot CNEs|PhastCons Score distribution and coverage

eFe2 <- ggplot(AtZm_Conservation_Scores, aes(x = Class, y = Score)) +
  geom_boxplot(lwd=2) + theme_bw() + coord_flip() + ylim(0,1) + 
  theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none") +
  labs( y = "LECIF score") + geom_hline(aes(yintercept=0.17), colour="darkgray") + geom_hline(aes(yintercept=0.15), colour="darkgray")


H2 <- ggplot(AtZm_Conservation_Coverage, aes(x = Class, y = Coverage))+
  geom_col(width = 0.7) + coord_flip() + 
  theme_bw() + ylab("Coverage (%)") + xlab("Conserved Elements") + ylim(0,100) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(data=AtZm_Conservation_Coverage,aes(x=Class,y=Coverage,label=CoverageRounder),vjust=0.5, hjust = -0.2) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

### 4.1 Arabidopsis block composition

try <- ggdraw() +
  draw_plot(B, x = 0.03, y = 0.8, width = 0.42, height = .20) +
  draw_plot(A, x = 0, y = 0.2, width = 0.45, height = .60) +
  draw_plot(D, x = 0.46, y = 0.8, width = 0.35, height = .20) +
  draw_plot(E, x = 0.46, y = 0.2, width = 0.35, height = .60) +
  draw_plot(C, x = 0.81, y = 0.2, width = 0.15, height = .60) + 
  draw_plot(eFe, x = 0, y = 0, width = 0.8, height = .20) +
  draw_plot(H, x = 0.81, y = 0, width = 0.15, height = .20)

try2 <- ggdraw() +
  draw_plot(B2, x = 0.03, y = 0.8, width = 0.42, height = .20) +
  draw_plot(A2, x = 0, y = 0.2, width = 0.45, height = .60) +
  draw_plot(D2, x = 0.46, y = 0.8, width = 0.35, height = .20) +
  draw_plot(E2, x = 0.46, y = 0.2, width = 0.35, height = .60) +
  draw_plot(C2, x = 0.81, y = 0.2, width = 0.15, height = .60) + 
  draw_plot(eFe2, x = 0, y = 0, width = 0.8, height = .20) +
  draw_plot(H2, x = 0.81, y = 0, width = 0.15, height = .20)

ggdraw() +
  draw_plot(try, x = 0, y = 0, width = 0.5, height = 1) +
  draw_plot(try2, x = 0.5, y = 0, width = 0.5, height = 1) # horizontal

### 4.1.4 Arabidopsis Variability

SNPs_AtOs <- ggplot(RegionResults_AtOs, aes(x=description, y=oddsRatio, fill=description)) + 
  geom_bar(stat='identity') +
  scale_fill_manual(values=c('gray60', 'gray60', 'gray60', 'gray60', 'gray60')) +
  theme_bw() +
  xlab("At LECIF score bins") + ylab("LOLA OddsRatio") + 
  ggtitle("Os LECIF score GWAS/Common Overlap Enrichments") + ylim(0,4) + theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

SNPs_AtZm <- ggplot(RegionResults_AtZm, aes(x=description, y=oddsRatio, fill=description)) + 
  geom_bar(stat='identity') +
  scale_fill_manual(values=c('gray60', 'gray60', 'gray60', 'gray60', 'gray60')) +
  theme_bw() +
  xlab("At LECIF score bins") + ylab("LOLA OddsRatio") + 
  ggtitle("Zm LECIF score GWAS/Common Overlap Enrichments") + ylim(0,4) + theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

Variability_At <- ggdraw() +
  draw_plot(SNPs_AtOs, x = 0, y = 0, width = 0.5, height = 1) +
  draw_plot(SNPs_AtZm, x = 0.5, y = 0, width = 0.5, height = 1)

### 4.2.1 Oryza-Arabidopsis Chromatin states-based

## Violin-plot CS score distribution

data <- data.frame(name = LECIF_OsAt$V5, value = LECIF_OsAt$V4)
data$names <- factor(data$name , levels=c("CS16", "CS15", "CS14", "CS13",
                                          "CS12", "CS11", "CS10", "CS9",
                                          "CS8", "CS7", "CS6", "CS5",
                                          "CS4", "CS3", "CS2", "CS1"), ordered = T)
r1<-with(data, tapply(value, name, mean))
r2<-with(data, tapply(value, name, median))
g <-
  ggplot(data, aes(x = names, y = value,
                   color = value)) +
  labs(x = "Chromatin states", y = "LECIF score") +
  scale_color_brewer(palette = "Dark2", guide = "none")
A3 <- g + geom_violin(aes(fill = names), size = 1, alpha = .5) +
  geom_boxplot(outlier.alpha = 0, coef = 0, color = "black", width = .2) +
  scale_fill_manual(values=c(rep("gray70",3),  rep("gold",3), rep("#4DBBD5FF",1), rep("#00A087FF",4), rep("#3C5488FF",5))) +
  coord_flip() + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  ylim(0,1)

## LECIF score distribution

B3 <- ggplot(data, aes(x=value)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="gray70") +
  geom_density(alpha=.4, fill="black") + theme_bw() + ggtitle("LECIF score distribution") + 
  labs(y = "Proportion") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) + 
  xlim(0,1) +
  annotate("text", x = 0.80, y = 7, label=paste0("Mean: ", round(mean(data$value), digits = 2)), size = 5, color = "black") +
  annotate("text", x = 0.82, y = 6, label=paste0("Median: ", round(median(data$value), digits = 2)), size = 5, color = "black")

## Coverage of each CS in across aligning regions (first species is the target genome)

aligningRegions <- sum(LECIF_OsAt$V3 - LECIF_OsAt$V2)
Coverage <- data.frame("Chromatin States" = c("CS1","CS2", "CS3", "CS4", "CS5", "CS6", "CS7", "CS8", "CS9", "CS10", "CS11", "CS12", "CS13", "CS14", "CS15", "CS16"),
                       Coverage = NA)
for(i in 1:nrow(Coverage)){
  Coverage[i,2] <- (sum(LECIF_OsAt$V3[which(LECIF_OsAt$V5 == Coverage[i,1])] - LECIF_OsAt$V2[which(LECIF_OsAt$V5 == Coverage[i,1])])/aligningRegions)*100 
}
Coverage$Chromatin.States <- factor(Coverage$Chromatin.States , levels=c("CS16", "CS15", "CS14", "CS13",
                                                                         "CS12", "CS11", "CS10", "CS9",
                                                                         "CS8", "CS7", "CS6", "CS5",
                                                                         "CS4", "CS3", "CS2", "CS1"), ordered = TRUE) 
Coverage$CoverageRounded <- round(Coverage$Coverage, digits = 2)
C3 <- ggplot(Coverage, aes(x = Chromatin.States, y = Coverage))+
  geom_col(width = 0.7) + coord_flip() + 
  theme_bw() + ylab("Coverage (%)") + xlab("Chromatin states") + ylim(0,100) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(data=Coverage,aes(x=Chromatin.States,y=Coverage,label=CoverageRounded),vjust=0.5, hjust = -0.2) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
### ADD TEXT CORRELATION


## Lineplot Chromatin State Functional group Simmilairty in three LECIF percentile bins

Lineplot.table.OsAt <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/OsAt_CSsimmilarity.DF.Distances3.txt", header = T, sep = "\t")
Lineplot.table.OsAt <- Lineplot.table.OsAt[c(13,14,15),]
Lineplot.table.OsAt <- Lineplot.table.OsAt[,c(19,20,21,22,23)]
Lineplot.df.OsAt <- data.frame(CSGroup = c(rep("Bivalent",3), rep("Active",3), rep("Divergent",3), rep("Heterochromatin",3), rep("Quies",3)),
                               LECIF.Bin = c(rep(c(0.5,1.5,2.5),5)),
                               Simmilarity = c(Lineplot.table.OsAt$Bivalent, Lineplot.table.OsAt$Active, Lineplot.table.OsAt$Divergent,
                                               Lineplot.table.OsAt$Heterochromatin, Lineplot.table.OsAt$Quies))
Lineplot.df.OsAt$Simmilarity <- 1 - Lineplot.df.OsAt$Simmilarity
Lineplot.df.OsAt$LECIF.Bin <- as.numeric(as.vector(Lineplot.df.OsAt$LECIF.Bin))
D3 <- ggplot(data=Lineplot.df.OsAt, aes(x=LECIF.Bin, y=Simmilarity, group=CSGroup, color=CSGroup)) +
  geom_line(size=1.0) + geom_point(size=3)+
  scale_color_manual(values=c("#00A087FF", "#3C5488FF", "#4DBBD5FF", "gold", "gray70"))+
  theme_bw() + coord_flip() + theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlim(0.5,2.5) + ylim(0,0.5) + ggtitle("Simmilarity in 3 LECIF bins")


## Grouped-barplot (two bars) Chromatin  state Simmilarity between high;low and low;high LECIF;PhyloP scores (50;50 | 60;40 cutoffs)

Barplot.table.OsAt <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/OsAt_CSsimmilarity.DF.Distances.After6040.txt", header = T, sep = "\t")
Barplot.table.OsAt <- Barplot.table.OsAt[c(19,20),]
Barplot.table.OsAt <- Barplot.table.OsAt[,c(1:18)]
Barplot.df.OsAt <- data.frame(CS = c(rep("CS1",2),rep("CS2",2),rep("CS3",2),rep("CS4",2),
                                     rep("CS5",2),rep("CS6",2),rep("CS7",2),rep("CS8",2),
                                     rep("CS9",2),rep("CS10",2),rep("CS11",2),rep("CS12",2),
                                     rep("CS13",2),rep("CS14",2),rep("CS15",2),rep("CS16",2)), Class = c(rep(c("High,Low", "Low,High"),16)), Simmilarity = c(Barplot.table.OsAt$CS1,Barplot.table.OsAt$CS2,Barplot.table.OsAt$CS3,Barplot.table.OsAt$CS4,
                                                                                                                                                             Barplot.table.OsAt$CS5,Barplot.table.OsAt$CS6,Barplot.table.OsAt$CS7,Barplot.table.OsAt$CS8,
                                                                                                                                                             Barplot.table.OsAt$CS9,Barplot.table.OsAt$CS10,Barplot.table.OsAt$CS11,Barplot.table.OsAt$CS12,
                                                                                                                                                             Barplot.table.OsAt$CS13,Barplot.table.OsAt$CS14,Barplot.table.OsAt$CS15,Barplot.table.OsAt$CS16))
Barplot.df.OsAt$Simmilarity <- 1 - Barplot.df.OsAt$Simmilarity
Barplot.df.OsAt$CS<-factor(Barplot.df.OsAt$CS, levels=c("CS16", "CS15", "CS14", "CS13",
                                                        "CS12", "CS11", "CS10", "CS9",
                                                        "CS8", "CS7", "CS6", "CS5",
                                                        "CS4", "CS3", "CS2", "CS1"), ordered = T)
Barplot.df.OsAt$Class <- as.factor(Barplot.df.OsAt$Class)
# Grouped
E3 <- ggplot(Barplot.df.OsAt, aes(fill=Class, y=Simmilarity, x=CS)) + 
  geom_bar(position="dodge", stat="identity") + coord_flip() + theme_bw() + ylim(0,0.5) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")


### Oryza-Arabidopsis Conservation Boxplot CNEs|PhastCons Score distribution and coverage

eFe3 <- ggplot(OsAt_Conservation_Scores, aes(x = Class, y = Score)) +
  geom_boxplot(lwd=2) + theme_bw() + coord_flip() + ylim(0,1) + 
  theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none") +
  labs( y = "LECIF score") + geom_hline(aes(yintercept=0.15), colour="darkgray") + geom_hline(aes(yintercept=0.16), colour="darkgray")


H3 <- ggplot(OsAt_Conservation_Coverage, aes(x = Class, y = Coverage))+
  geom_col(width = 0.7) + coord_flip() + 
  theme_bw() + ylab("Coverage (%)") + xlab("Conserved Elements") + ylim(0,100) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(data=OsAt_Conservation_Coverage,aes(x=Class,y=Coverage,label=CoverageRounder),vjust=0.5, hjust = -0.2) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())


### 4.2.2 Oryza-Zea Chromatin states-based (add correlation text)

## Violin-plot CS score distribution

data <- data.frame(name = LECIF_OsZm$V5, value = LECIF_OsZm$V4)
data$names <- factor(data$name , levels=c("CS16", "CS15", "CS14", "CS13",
                                          "CS12", "CS11", "CS10", "CS9",
                                          "CS8", "CS7", "CS6", "CS5",
                                          "CS4", "CS3", "CS2", "CS1"), ordered = T)
r1<-with(data, tapply(value, name, mean))
r2<-with(data, tapply(value, name, median))
g <-
  ggplot(data, aes(x = names, y = value,
                   color = value)) +
  labs(x = "Chromatin states", y = "LECIF score") +
  scale_color_brewer(palette = "Dark2", guide = "none")
A4 <- g + geom_violin(aes(fill = names), size = 1, alpha = .5) +
  geom_boxplot(outlier.alpha = 0, coef = 0, color = "black", width = .2) +
  scale_fill_manual(values=c(rep("gray70",3),  rep("gold",3), rep("#4DBBD5FF",1), rep("#00A087FF",4), rep("#3C5488FF",5))) +
  coord_flip() + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  ylim(0,1)

## LECIF score distribution

B4 <- ggplot(data, aes(x=value)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="gray70") +
  geom_density(alpha=.4, fill="black") + theme_bw() + ggtitle("LECIF score distribution") + 
  labs(y = "Proportion") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) + 
  xlim(0,1) +
  annotate("text", x = 0.80, y = 7, label=paste0("Mean: ", round(mean(data$value), digits = 2)), size = 5, color = "black") +
  annotate("text", x = 0.82, y = 6, label=paste0("Median: ", round(median(data$value), digits = 2)), size = 5, color = "black")

## Coverage of each CS in across aligning regions (first species is the target genome)

aligningRegions <- sum(LECIF_OsZm$V3 - LECIF_OsZm$V2)
Coverage <- data.frame("Chromatin States" = c("CS1","CS2", "CS3", "CS4", "CS5", "CS6", "CS7", "CS8", "CS9", "CS10", "CS11", "CS12", "CS13", "CS14", "CS15", "CS16"),
                       Coverage = NA)
for(i in 1:nrow(Coverage)){
  Coverage[i,2] <- (sum(LECIF_OsZm$V3[which(LECIF_OsZm$V5 == Coverage[i,1])] - LECIF_OsZm$V2[which(LECIF_OsZm$V5 == Coverage[i,1])])/aligningRegions)*100 
}
Coverage$Chromatin.States <- factor(Coverage$Chromatin.States , levels=c("CS16", "CS15", "CS14", "CS13",
                                                                         "CS12", "CS11", "CS10", "CS9",
                                                                         "CS8", "CS7", "CS6", "CS5",
                                                                         "CS4", "CS3", "CS2", "CS1"), ordered = TRUE) 
Coverage$CoverageRounded <- round(Coverage$Coverage, digits = 2)
C4 <- ggplot(Coverage, aes(x = Chromatin.States, y = Coverage))+
  geom_col(width = 0.7) + coord_flip() + 
  theme_bw() + ylab("Coverage (%)") + xlab("Chromatin states") + ylim(0,100) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(data=Coverage,aes(x=Chromatin.States,y=Coverage,label=CoverageRounded),vjust=0.5, hjust = -0.2) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
### ADD TEXT CORRELATION


## Lineplot Chromatin State Functional group Simmilairty in three LECIF percentile bins

Lineplot.table.OsZm <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/OsZm_CSsimmilarity.DF.Distances3.txt", header = T, sep = "\t")
Lineplot.table.OsZm <- Lineplot.table.OsZm[c(13,14,15),]
Lineplot.table.OsZm <- Lineplot.table.OsZm[,c(19,20,21,22,23)]
Lineplot.df.OsZm <- data.frame(CSGroup = c(rep("Bivalent",3), rep("Active",3), rep("Divergent",3), rep("Heterochromatin",3), rep("Quies",3)),
                               LECIF.Bin = c(rep(c(0.5,1.5,2.5),5)),
                               Simmilarity = c(Lineplot.table.OsZm$Bivalent, Lineplot.table.OsZm$Active, Lineplot.table.OsZm$Divergent,
                                               Lineplot.table.OsZm$Heterochromatin, Lineplot.table.OsZm$Quies))
Lineplot.df.OsZm$Simmilarity <- 1 - Lineplot.df.OsZm$Simmilarity
Lineplot.df.OsZm$LECIF.Bin <- as.numeric(as.vector(Lineplot.df.OsZm$LECIF.Bin))
D4 <- ggplot(data=Lineplot.df.OsZm, aes(x=LECIF.Bin, y=Simmilarity, group=CSGroup, color=CSGroup)) +
  geom_line(size=1.0) + geom_point(size=3)+
  scale_color_manual(values=c("#00A087FF", "#3C5488FF", "#4DBBD5FF", "gold", "gray70"))+
  theme_bw() + coord_flip() + theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlim(0.5,2.5) + ylim(0,0.57) + ggtitle("Simmilarity in 3 LECIF bins")


## Grouped-barplot (two bars) Chromatin  state Simmilarity between high;low and low;high LECIF;PhyloP scores (50;50 | 60;40 cutoffs)

Barplot.table.OsZm <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/OsZm_CSsimmilarity.DF.Distances.After6040.txt", header = T, sep = "\t")
Barplot.table.OsZm <- Barplot.table.OsZm[c(19,20),]
Barplot.table.OsZm <- Barplot.table.OsZm[,c(1:18)]
Barplot.df.OsZm <- data.frame(CS = c(rep("CS1",2),rep("CS2",2),rep("CS3",2),rep("CS4",2),
                                     rep("CS5",2),rep("CS6",2),rep("CS7",2),rep("CS8",2),
                                     rep("CS9",2),rep("CS10",2),rep("CS11",2),rep("CS12",2),
                                     rep("CS13",2),rep("CS14",2),rep("CS15",2),rep("CS16",2)), Class = c(rep(c("High,Low", "Low,High"),16)), Simmilarity = c(Barplot.table.OsZm$CS1,Barplot.table.OsZm$CS2,Barplot.table.OsZm$CS3,Barplot.table.OsZm$CS4,
                                                                                                                                                             Barplot.table.OsZm$CS5,Barplot.table.OsZm$CS6,Barplot.table.OsZm$CS7,Barplot.table.OsZm$CS8,
                                                                                                                                                             Barplot.table.OsZm$CS9,Barplot.table.OsZm$CS10,Barplot.table.OsZm$CS11,Barplot.table.OsZm$CS12,
                                                                                                                                                             Barplot.table.OsZm$CS13,Barplot.table.OsZm$CS14,Barplot.table.OsZm$CS15,Barplot.table.OsZm$CS16))
Barplot.df.OsZm$Simmilarity <- 1 - Barplot.df.OsZm$Simmilarity
Barplot.df.OsZm$CS<-factor(Barplot.df.OsZm$CS, levels=c("CS16", "CS15", "CS14", "CS13",
                                                        "CS12", "CS11", "CS10", "CS9",
                                                        "CS8", "CS7", "CS6", "CS5",
                                                        "CS4", "CS3", "CS2", "CS1"), ordered = T)
Barplot.df.OsZm$Class <- as.factor(Barplot.df.OsZm$Class)
# Grouped
E4 <- ggplot(Barplot.df.OsZm, aes(fill=Class, y=Simmilarity, x=CS)) + 
  geom_bar(position="dodge", stat="identity") + coord_flip() + theme_bw() + ylim(0,0.5) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")


### Oryza-Zea Conservation Boxplot CNEs|PhastCons Score distribution and coverage

eFe4 <- ggplot(OsZm_Conservation_Scores, aes(x = Class, y = Score)) +
  geom_boxplot(lwd=2) + theme_bw() + coord_flip() + ylim(0,1) + 
  theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none") +
  labs( y = "LECIF score") + geom_hline(aes(yintercept=0.19), colour="darkgray") + geom_hline(aes(yintercept=0.2), colour="darkgray")


H4 <- ggplot(OsZm_Conservation_Coverage, aes(x = Class, y = Coverage))+
  geom_col(width = 0.7) + coord_flip() + 
  theme_bw() + ylab("Coverage (%)") + xlab("Conserved Elements") + ylim(0,100) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(data=OsZm_Conservation_Coverage,aes(x=Class,y=Coverage,label=CoverageRounder),vjust=0.5, hjust = -0.2) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

### 4.2 Oryza block composition

try3 <- ggdraw() +
  draw_plot(B3, x = 0.03, y = 0.8, width = 0.42, height = .20) +
  draw_plot(A3, x = 0, y = 0.2, width = 0.45, height = .60) +
  draw_plot(D3, x = 0.46, y = 0.8, width = 0.35, height = .20) +
  draw_plot(E3, x = 0.46, y = 0.2, width = 0.35, height = .60) +
  draw_plot(C3, x = 0.81, y = 0.2, width = 0.15, height = .60) + 
  draw_plot(eFe3, x = 0, y = 0, width = 0.8, height = .20) +
  draw_plot(H3, x = 0.81, y = 0, width = 0.15, height = .20)

try4 <- ggdraw() +
  draw_plot(B4, x = 0.03, y = 0.8, width = 0.42, height = .20) +
  draw_plot(A4, x = 0, y = 0.2, width = 0.45, height = .60) +
  draw_plot(D4, x = 0.46, y = 0.8, width = 0.35, height = .20) +
  draw_plot(E4, x = 0.46, y = 0.2, width = 0.35, height = .60) +
  draw_plot(C4, x = 0.81, y = 0.2, width = 0.15, height = .60) + 
  draw_plot(eFe4, x = 0, y = 0, width = 0.8, height = .20) +
  draw_plot(H4, x = 0.81, y = 0, width = 0.15, height = .20)

ggdraw() +
  draw_plot(try3, x = 0, y = 0, width = 0.5, height = 1) +
  draw_plot(try4, x = 0.5, y = 0, width = 0.5, height = 1) # horizontal

### 4.2.4 Oryza Variability

SNPs_OsAt <- ggplot(RegionResults_OsAt, aes(x=description, y=oddsRatio, fill=description)) + 
  geom_bar(stat='identity') +
  scale_fill_manual(values=c('black', 'black', 'gray60', 'gray60', 'gray60')) +
  theme_bw() +
  xlab("Os LECIF score bins") + ylab("LOLA OddsRatio") + 
  ggtitle("At LECIF score GWAS/Common Overlap Enrichments") + ylim(0,4) + theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

SNPs_OsZm <- ggplot(RegionResults_OsZm, aes(x=description, y=oddsRatio, fill=description)) + 
  geom_bar(stat='identity') +
  scale_fill_manual(values=c('gray60', 'gray60', 'gray60', 'black', 'black')) +
  theme_bw() +
  xlab("Os LECIF score bins") + ylab("LOLA OddsRatio") + 
  ggtitle("Zm LECIF score GWAS/Common Overlap Enrichments") + ylim(0,4) + theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

Variability_Os <- ggdraw() +
  draw_plot(SNPs_OsAt, x = 0, y = 0, width = 0.5, height = 1) +
  draw_plot(SNPs_OsZm, x = 0.5, y = 0, width = 0.5, height = 1)

### 4.3.1 Zea-Arabidopsis Chromatin states-based

## Violin-plot CS score distribution

data <- data.frame(name = LECIF_ZmAt$V5, value = LECIF_ZmAt$V4)
data$names <- factor(data$name , levels=c("CS16", "CS15", "CS14", "CS13",
                                          "CS12", "CS11", "CS10", "CS9",
                                          "CS8", "CS7", "CS6", "CS5",
                                          "CS4", "CS3", "CS2", "CS1"), ordered = T)
r1<-with(data, tapply(value, name, mean))
r2<-with(data, tapply(value, name, median))
g <-
  ggplot(data, aes(x = names, y = value,
                   color = value)) +
  labs(x = "Chromatin states", y = "LECIF score") +
  scale_color_brewer(palette = "Dark2", guide = "none")
A5 <- g + geom_violin(aes(fill = names), size = 1, alpha = .5) +
  geom_boxplot(outlier.alpha = 0, coef = 0, color = "black", width = .2) +
  scale_fill_manual(values=c(rep("gray70",3),  rep("gold",3), rep("#4DBBD5FF",1), rep("#00A087FF",4), rep("#3C5488FF",5))) +
  coord_flip() + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  ylim(0,1)

## LECIF score distribution

B5 <- ggplot(data, aes(x=value)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="gray70") +
  geom_density(alpha=.4, fill="black") + theme_bw() + ggtitle("LECIF score distribution") + 
  labs(y = "Proportion") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) + 
  xlim(0,1) +
  annotate("text", x = 0.80, y = 7, label=paste0("Mean: ", round(mean(data$value), digits = 2)), size = 5, color = "black") +
  annotate("text", x = 0.82, y = 6, label=paste0("Median: ", round(median(data$value), digits = 2)), size = 5, color = "black")

## Coverage of each CS in across aligning regions (first species is the target genome)

aligningRegions <- sum(LECIF_ZmAt$V3 - LECIF_ZmAt$V2)
Coverage <- data.frame("Chromatin States" = c("CS1","CS2", "CS3", "CS4", "CS5", "CS6", "CS7", "CS8", "CS9", "CS10", "CS11", "CS12", "CS13", "CS14", "CS15", "CS16"),
                       Coverage = NA)
for(i in 1:nrow(Coverage)){
  Coverage[i,2] <- (sum(LECIF_ZmAt$V3[which(LECIF_ZmAt$V5 == Coverage[i,1])] - LECIF_ZmAt$V2[which(LECIF_ZmAt$V5 == Coverage[i,1])])/aligningRegions)*100 
}
Coverage$Chromatin.States <- factor(Coverage$Chromatin.States , levels=c("CS16", "CS15", "CS14", "CS13",
                                                                         "CS12", "CS11", "CS10", "CS9",
                                                                         "CS8", "CS7", "CS6", "CS5",
                                                                         "CS4", "CS3", "CS2", "CS1"), ordered = TRUE) 
Coverage$CoverageRounded <- round(Coverage$Coverage, digits = 2)
C5 <- ggplot(Coverage, aes(x = Chromatin.States, y = Coverage))+
  geom_col(width = 0.7) + coord_flip() + 
  theme_bw() + ylab("Coverage (%)") + xlab("Chromatin states") + ylim(0,100) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(data=Coverage,aes(x=Chromatin.States,y=Coverage,label=CoverageRounded),vjust=0.5, hjust = -0.2) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
### ADD TEXT CORRELATION


## Lineplot Chromatin State Functional group Simmilairty in three LECIF percentile bins

Lineplot.table.ZmAt <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/ZmAt_CSsimmilarity.DF.Distances3.txt", header = T, sep = "\t")
Lineplot.table.ZmAt <- Lineplot.table.ZmAt[c(13,14,15),]
Lineplot.table.ZmAt <- Lineplot.table.ZmAt[,c(19,20,21,22,23)]
Lineplot.df.ZmAt <- data.frame(CSGroup = c(rep("Bivalent",3), rep("Active",3), rep("Divergent",3), rep("Heterochromatin",3), rep("Quies",3)),
                               LECIF.Bin = c(rep(c(0.5,1.5,2.5),5)),
                               Simmilarity = c(Lineplot.table.ZmAt$Bivalent, Lineplot.table.ZmAt$Active, Lineplot.table.ZmAt$Divergent,
                                               Lineplot.table.ZmAt$Heterochromatin, Lineplot.table.ZmAt$Quies))
Lineplot.df.ZmAt$Simmilarity <- 1 - Lineplot.df.ZmAt$Simmilarity
Lineplot.df.ZmAt$LECIF.Bin <- as.numeric(as.vector(Lineplot.df.ZmAt$LECIF.Bin))
D5 <- ggplot(data=Lineplot.df.ZmAt, aes(x=LECIF.Bin, y=Simmilarity, group=CSGroup, color=CSGroup)) +
  geom_line(size=1.0) + geom_point(size=3)+
  scale_color_manual(values=c("#00A087FF", "#3C5488FF", "#4DBBD5FF", "gold", "gray70"))+
  theme_bw() + coord_flip() + theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlim(0.5,2.5) + ylim(0,0.5) + ggtitle("Simmilarity in 3 LECIF bins")


## Grouped-barplot (two bars) Chromatin  state Simmilarity between high;low and low;high LECIF;PhyloP scores (50;50 | 60;40 cutoffs)

Barplot.table.ZmAt <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/ZmAt_CSsimmilarity.DF.Distances.After6040.txt", header = T, sep = "\t")
Barplot.table.ZmAt <- Barplot.table.ZmAt[c(19,20),]
Barplot.table.ZmAt <- Barplot.table.ZmAt[,c(1:18)]
Barplot.df.ZmAt <- data.frame(CS = c(rep("CS1",2),rep("CS2",2),rep("CS3",2),rep("CS4",2),
                                     rep("CS5",2),rep("CS6",2),rep("CS7",2),rep("CS8",2),
                                     rep("CS9",2),rep("CS10",2),rep("CS11",2),rep("CS12",2),
                                     rep("CS13",2),rep("CS14",2),rep("CS15",2),rep("CS16",2)), Class = c(rep(c("High,Low", "Low,High"),16)), Simmilarity = c(Barplot.table.ZmAt$CS1,Barplot.table.ZmAt$CS2,Barplot.table.ZmAt$CS3,Barplot.table.ZmAt$CS4,
                                                                                                                                                             Barplot.table.ZmAt$CS5,Barplot.table.ZmAt$CS6,Barplot.table.ZmAt$CS7,Barplot.table.ZmAt$CS8,
                                                                                                                                                             Barplot.table.ZmAt$CS9,Barplot.table.ZmAt$CS10,Barplot.table.ZmAt$CS11,Barplot.table.ZmAt$CS12,
                                                                                                                                                             Barplot.table.ZmAt$CS13,Barplot.table.ZmAt$CS14,Barplot.table.ZmAt$CS15,Barplot.table.ZmAt$CS16))
Barplot.df.ZmAt$Simmilarity <- 1 - Barplot.df.ZmAt$Simmilarity
Barplot.df.ZmAt$CS<-factor(Barplot.df.ZmAt$CS, levels=c("CS16", "CS15", "CS14", "CS13",
                                                        "CS12", "CS11", "CS10", "CS9",
                                                        "CS8", "CS7", "CS6", "CS5",
                                                        "CS4", "CS3", "CS2", "CS1"), ordered = T)
Barplot.df.ZmAt$Class <- as.factor(Barplot.df.ZmAt$Class)
# Grouped
E5 <- ggplot(Barplot.df.ZmAt, aes(fill=Class, y=Simmilarity, x=CS)) + 
  geom_bar(position="dodge", stat="identity") + coord_flip() + theme_bw() + ylim(0,0.5) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")


### Zea-Arabidopsis Conservation Boxplot CNEs|PhastCons Score distribution and coverage

eFe5 <- ggplot(ZmAt_Conservation_Scores, aes(x = Class, y = Score)) +
  geom_boxplot(lwd=2) + theme_bw() + coord_flip() + ylim(0,1) + 
  theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none") +
  labs( y = "LECIF score") + geom_hline(aes(yintercept=0.15), colour="darkgray") + geom_hline(aes(yintercept=0.14), colour="darkgray")


H5 <- ggplot(ZmAt_Conservation_Coverage, aes(x = Class, y = Coverage))+
  geom_col(width = 0.7) + coord_flip() + 
  theme_bw() + ylab("Coverage (%)") + xlab("Conserved Elements") + ylim(0,100) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(data=ZmAt_Conservation_Coverage,aes(x=Class,y=Coverage,label=CoverageRounder),vjust=0.5, hjust = -0.2) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

### 4.3.2 Zea-Oryza Chromatin states-based (add correlation text)

## Violin-plot CS score distribution

data <- data.frame(name = LECIF_ZmOs$V5, value = LECIF_ZmOs$V4)
data$names <- factor(data$name , levels=c("CS16", "CS15", "CS14", "CS13",
                                          "CS12", "CS11", "CS10", "CS9",
                                          "CS8", "CS7", "CS6", "CS5",
                                          "CS4", "CS3", "CS2", "CS1"), ordered = T)
r1<-with(data, tapply(value, name, mean))
r2<-with(data, tapply(value, name, median))
g <-
  ggplot(data, aes(x = names, y = value,
                   color = value)) +
  labs(x = "Chromatin states", y = "LECIF score") +
  scale_color_brewer(palette = "Dark2", guide = "none")
A6 <- g + geom_violin(aes(fill = names), size = 1, alpha = .5) +
  geom_boxplot(outlier.alpha = 0, coef = 0, color = "black", width = .2) +
  scale_fill_manual(values=c(rep("gray70",3),  rep("gold",3), rep("#4DBBD5FF",1), rep("#00A087FF",4), rep("#3C5488FF",5))) +
  coord_flip() + theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
  ylim(0,1)

## LECIF score distribution

B6 <- ggplot(data, aes(x=value)) + 
  geom_histogram(aes(y=..density..), colour="black", fill="gray70") +
  geom_density(alpha=.4, fill="black") + theme_bw() + ggtitle("LECIF score distribution") + 
  labs(y = "Proportion") +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank()) + 
  xlim(0,1) +
  annotate("text", x = 0.80, y = 7, label=paste0("Mean: ", round(mean(data$value), digits = 2)), size = 5, color = "black") +
  annotate("text", x = 0.82, y = 6, label=paste0("Median: ", round(median(data$value), digits = 2)), size = 5, color = "black")

## Coverage of each CS in across aligning regions (first species is the target genome)

aligningRegions <- sum(LECIF_ZmOs$V3 - LECIF_ZmOs$V2)
Coverage <- data.frame("Chromatin States" = c("CS1","CS2", "CS3", "CS4", "CS5", "CS6", "CS7", "CS8", "CS9", "CS10", "CS11", "CS12", "CS13", "CS14", "CS15", "CS16"),
                       Coverage = NA)
for(i in 1:nrow(Coverage)){
  Coverage[i,2] <- (sum(LECIF_ZmOs$V3[which(LECIF_ZmOs$V5 == Coverage[i,1])] - LECIF_ZmOs$V2[which(LECIF_ZmOs$V5 == Coverage[i,1])])/aligningRegions)*100 
}
Coverage$Chromatin.States <- factor(Coverage$Chromatin.States , levels=c("CS16", "CS15", "CS14", "CS13",
                                                                         "CS12", "CS11", "CS10", "CS9",
                                                                         "CS8", "CS7", "CS6", "CS5",
                                                                         "CS4", "CS3", "CS2", "CS1"), ordered = TRUE) 
Coverage$CoverageRounded <- round(Coverage$Coverage, digits = 2)
C6 <- ggplot(Coverage, aes(x = Chromatin.States, y = Coverage))+
  geom_col(width = 0.7) + coord_flip() + 
  theme_bw() + ylab("Coverage (%)") + xlab("Chromatin states") + ylim(0,100) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(data=Coverage,aes(x=Chromatin.States,y=Coverage,label=CoverageRounded),vjust=0.5, hjust = -0.2) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
### ADD TEXT CORRELATION


## Lineplot Chromatin State Functional group Simmilairty in three LECIF percentile bins

Lineplot.table.ZmOs <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/ZmOs_CSsimmilarity.DF.Distances3.txt", header = T, sep = "\t")
Lineplot.table.ZmOs <- Lineplot.table.ZmOs[c(13,14,15),]
Lineplot.table.ZmOs <- Lineplot.table.ZmOs[,c(19,20,21,22,23)]
Lineplot.df.ZmOs <- data.frame(CSGroup = c(rep("Bivalent",3), rep("Active",3), rep("Divergent",3), rep("Heterochromatin",3), rep("Quies",3)),
                               LECIF.Bin = c(rep(c(0.5,1.5,2.5),5)),
                               Simmilarity = c(Lineplot.table.ZmOs$Bivalent, Lineplot.table.ZmOs$Active, Lineplot.table.ZmOs$Divergent,
                                               Lineplot.table.ZmOs$Heterochromatin, Lineplot.table.ZmOs$Quies))
Lineplot.df.ZmOs$Simmilarity <- 1 - Lineplot.df.ZmOs$Simmilarity
Lineplot.df.ZmOs$LECIF.Bin <- as.numeric(as.vector(Lineplot.df.ZmOs$LECIF.Bin))
D6 <- ggplot(data=Lineplot.df.ZmOs, aes(x=LECIF.Bin, y=Simmilarity, group=CSGroup, color=CSGroup)) +
  geom_line(size=1.0) + geom_point(size=3)+
  scale_color_manual(values=c("#00A087FF", "#3C5488FF", "#4DBBD5FF", "gold", "gray70"))+
  theme_bw() + coord_flip() + theme(axis.title.x = element_blank(), legend.position = "none", axis.title.y = element_blank(),plot.title = element_text(hjust = 0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank(),axis.text.x = element_blank(), axis.ticks.x = element_blank()) + xlim(0.5,2.5) + ylim(0,0.5) + ggtitle("Simmilarity in 3 LECIF bins")


## Grouped-barplot (two bars) Chromatin  state Simmilarity between high;low and low;high LECIF;PhyloP scores (50;50 | 60;40 cutoffs)

Barplot.table.ZmOs <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/ZmOs_CSsimmilarity.DF.Distances.After6040.txt", header = T, sep = "\t")
Barplot.table.ZmOs <- Barplot.table.ZmOs[c(19,20),]
Barplot.table.ZmOs <- Barplot.table.ZmOs[,c(1:18)]
Barplot.df.ZmOs <- data.frame(CS = c(rep("CS1",2),rep("CS2",2),rep("CS3",2),rep("CS4",2),
                                     rep("CS5",2),rep("CS6",2),rep("CS7",2),rep("CS8",2),
                                     rep("CS9",2),rep("CS10",2),rep("CS11",2),rep("CS12",2),
                                     rep("CS13",2),rep("CS14",2),rep("CS15",2),rep("CS16",2)), Class = c(rep(c("High,Low", "Low,High"),16)), Simmilarity = c(Barplot.table.ZmOs$CS1,Barplot.table.ZmOs$CS2,Barplot.table.ZmOs$CS3,Barplot.table.ZmOs$CS4,
                                                                                                                                                             Barplot.table.ZmOs$CS5,Barplot.table.ZmOs$CS6,Barplot.table.ZmOs$CS7,Barplot.table.ZmOs$CS8,
                                                                                                                                                             Barplot.table.ZmOs$CS9,Barplot.table.ZmOs$CS10,Barplot.table.ZmOs$CS11,Barplot.table.ZmOs$CS12,
                                                                                                                                                             Barplot.table.ZmOs$CS13,Barplot.table.ZmOs$CS14,Barplot.table.ZmOs$CS15,Barplot.table.ZmOs$CS16))
Barplot.df.ZmOs$Simmilarity <- 1 - Barplot.df.ZmOs$Simmilarity
Barplot.df.ZmOs$CS<-factor(Barplot.df.ZmOs$CS, levels=c("CS16", "CS15", "CS14", "CS13",
                                                        "CS12", "CS11", "CS10", "CS9",
                                                        "CS8", "CS7", "CS6", "CS5",
                                                        "CS4", "CS3", "CS2", "CS1"), ordered = T)
Barplot.df.ZmOs$Class <- as.factor(Barplot.df.ZmOs$Class)
# Grouped
E6 <- ggplot(Barplot.df.ZmOs, aes(fill=Class, y=Simmilarity, x=CS)) + 
  geom_bar(position="dodge", stat="identity") + coord_flip() + theme_bw() + ylim(0,0.5) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")


### Zea-Arabidopsis Conservation Boxplot CNEs|PhastCons Score distribution and coverage

eFe6 <- ggplot(ZmOs_Conservation_Scores, aes(x = Class, y = Score)) +
  geom_boxplot(lwd=2) + theme_bw() + coord_flip() + ylim(0,1) + 
  theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none") +
  labs( y = "LECIF score") + geom_hline(aes(yintercept=0.2), colour="darkgray") + geom_hline(aes(yintercept=0.14), colour="darkgray")


H6 <- ggplot(ZmOs_Conservation_Coverage, aes(x = Class, y = Coverage))+
  geom_col(width = 0.7) + coord_flip() + 
  theme_bw() + ylab("Coverage (%)") + xlab("Conserved Elements") + ylim(0,100) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(data=ZmOs_Conservation_Coverage,aes(x=Class,y=Coverage,label=CoverageRounder),vjust=0.5, hjust = -0.2) + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

### 4.3 Zea block composition

try5 <- ggdraw() +
  draw_plot(B5, x = 0.03, y = 0.8, width = 0.42, height = .20) +
  draw_plot(A5, x = 0, y = 0.2, width = 0.45, height = .60) +
  draw_plot(D5, x = 0.46, y = 0.8, width = 0.35, height = .20) +
  draw_plot(E5, x = 0.46, y = 0.2, width = 0.35, height = .60) +
  draw_plot(C5, x = 0.81, y = 0.2, width = 0.15, height = .60) + 
  draw_plot(eFe5, x = 0, y = 0, width = 0.8, height = .20) +
  draw_plot(H5, x = 0.81, y = 0, width = 0.15, height = .20)

try6 <- ggdraw() +
  draw_plot(B6, x = 0.03, y = 0.8, width = 0.42, height = .20) +
  draw_plot(A6, x = 0, y = 0.2, width = 0.45, height = .60) +
  draw_plot(D6, x = 0.46, y = 0.8, width = 0.35, height = .20) +
  draw_plot(E6, x = 0.46, y = 0.2, width = 0.35, height = .60) +
  draw_plot(C6, x = 0.81, y = 0.2, width = 0.15, height = .60) + 
  draw_plot(eFe6, x = 0, y = 0, width = 0.8, height = .20) +
  draw_plot(H6, x = 0.81, y = 0, width = 0.15, height = .20)

ggdraw() +
  draw_plot(try5, x = 0, y = 0, width = 0.5, height = 1) +
  draw_plot(try6, x = 0.5, y = 0, width = 0.5, height = 1) # horizontal

### 4.3.4 Zea Variability

SNPs_ZmAt <- ggplot(RegionResults_ZmAt, aes(x=description, y=oddsRatio, fill=description)) + 
  geom_bar(stat='identity') +
  scale_fill_manual(values=c('gray60', 'gray60', 'gray60', 'black', 'black')) +
  theme_bw() +
  xlab("Zm LECIF score bins") + ylab("LOLA OddsRatio") + 
  ggtitle("At LECIF score GWAS/Common Overlap Enrichments") + ylim(0,4) + theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

SNPs_ZmOs <- ggplot(RegionResults_ZmOs, aes(x=description, y=oddsRatio, fill=description)) + 
  geom_bar(stat='identity') +
  scale_fill_manual(values=c('gray60', 'gray60', 'gray60', 'black', 'black')) +
  theme_bw() +
  xlab("Zm LECIF score bins") + ylab("LOLA OddsRatio") + 
  ggtitle("Os LECIF score GWAS/Common Overlap Enrichments") + ylim(0,4) + theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

Variability_Zm <- ggdraw() +
  draw_plot(SNPs_ZmAt, x = 0, y = 0, width = 0.5, height = 1) +
  draw_plot(SNPs_ZmOs, x = 0.5, y = 0, width = 0.5, height = 1)

###############################################
####### 5. CIRCOS FINAL REPRESENTATION ######## 
##############################################

### 5. Circos representation of LECIF score,PhyloP score + synteny interactions + states + gene density
## At Blue/Green || Os Gold || Zm Purple

### Option 1) Concatenating genomes At*10 and Os*7 | Lines or nested with points of different color
### TRY make wide windows and map LECIF to that ones || TRACK GAP

## 1 Initialize

At_ChromInfo <- read.table("H:/Taiotactical/Taiotactical_ChIPseq/Utils/Genomes/0.At/At_genome.sizes")
Os_ChromInfo <- read.table("H:/Taiotactical/Taiotactical_ChIPseq/Utils/Genomes/0.Os/Os_genome.sizes")
Zm_ChromInfo <- read.table("H:/Taiotactical/Taiotactical_ChIPseq/Utils/Genomes/0.Zm/Zm_genome_circos.sizes")
At_ChromInfo$V1 <- paste0("At_", At_ChromInfo$V1)
Os_ChromInfo$V1 <- paste0("Os_", Os_ChromInfo$V1)
Zm_ChromInfo$V1 <- paste0("Zm_", Zm_ChromInfo$V1)
ChromInfo <- rbind(At_ChromInfo, Os_ChromInfo, Zm_ChromInfo)
df <- data.frame(name = ChromInfo$V1, start = rep(0, nrow(ChromInfo)), end = ChromInfo$V2)
df$end[1:5] <- df$end[1:5]*10
df$end[6:19] <- df$end[6:19]*7
circos.par(gap.after = c(rep(2, 4), 30, rep(2, 13), 30, rep(2, 9),30),start.degree= 10)
circos.genomicInitialize(df)

colors <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "yellow", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2", "gray71")

#2) Gene density

At <- makeTxDbFromGFF(paste0("H:/Taiotactical/Taiotactical_ChIPseq/Utils/Annotation/0.At/Annotation/0.At.gene_exons.gff3"))
Os <- makeTxDbFromGFF(paste0("H:/Taiotactical/Taiotactical_ChIPseq/Utils/Annotation/0.Os/Annotation/0.Os.gene_exons.gff3"))
Zm <- makeTxDbFromGFF(paste0("H:/Taiotactical/Taiotactical_ChIPseq/Utils/Annotation/0.Zm/Annotation/0.Zm.gene_exons.gff3"))

Os_Genes <- genes(Os)
Zm_Genes <- genes(Zm)
At_Genes <- data.frame(chr = as.character(At_Genes@seqnames), start = At_Genes@ranges@start*10, end = (At_Genes@ranges@start + At_Genes@ranges@width - 1)*10)
At_Genes$chr <- paste0("At_", At_Genes$chr)
Os_Genes <- data.frame(chr = as.character(Os_Genes@seqnames), start = Os_Genes@ranges@start*7, end = (Os_Genes@ranges@start + Os_Genes@ranges@width - 1)*7)
Os_Genes$chr <- paste0("Os_", Os_Genes$chr)
Zm_Genes <- data.frame(chr = as.character(Zm_Genes@seqnames), start = Zm_Genes@ranges@start, end = (Zm_Genes@ranges@start + Zm_Genes@ranges@width - 1))
Zm_Genes$chr <- paste0("Zm_", Zm_Genes$chr)
Genes <- rbind(At_Genes, Os_Genes, Zm_Genes)
colnames(Genes) <- c("chr", "start", "end")
Genes$value <- rep(1, nrow(Genes))
circos.genomicTrack(Genes, ylim=c(0,1),cell.padding=c(0,0,0,0),track.height = uh(1.5, "mm"),bg.border = colors[1],track.margin=c(0.001,0.001),
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, type = "h",col=colors[1],border=NA)
                    })
set_track_gap(gap = 0.02)

#2) LECIF At

LECIF_OsAt_Circos <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/Os-At_LECIFscore.merged.bedgraph", header = F)
LECIF_ZmAt_Circos <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/Zm-At_LECIFscore.merged.bedgraph", header = F)
LECIF_OsAt_Circos$V1 <- paste0("Os_",LECIF_OsAt_Circos$V1)
LECIF_ZmAt_Circos$V1 <- paste0("Zm_",LECIF_ZmAt_Circos$V1)
LECIF_OsAt_Circos$V2 <- LECIF_OsAt_Circos$V2*7
LECIF_OsAt_Circos$V3 <- LECIF_OsAt_Circos$V3*7
LECIF_At <- rbind(LECIF_OsAt_Circos, LECIF_ZmAt_Circos)
colnames(LECIF_At) <- c("chr", "start", "end", "value1")

circos.genomicTrack(LECIF_At, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, type = "l", area = TRUE, col = colors[7], border = colors[7], ...)
})

circos.genomicTrack(LECIF_At,panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, col = colors[7], pch = 16, cex = 0.5, ...)
})

circos.genomicDensity(LECIF_At[,1:3], col = colors[7], track.height = 0.1)
set_track_gap(gap = 0.02)

#3) LECIF Os

LECIF_AtOs_Circos <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/At-Os_LECIFscore.merged.bedgraph", header = F)
LECIF_ZmOs_Circos <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/Zm-Os_LECIFscore.merged.bedgraph", header = F)
LECIF_AtOs_Circos$V1 <- paste0("At_",LECIF_AtOs_Circos$V1)
LECIF_ZmOs_Circos$V1 <- paste0("Zm_",LECIF_ZmOs_Circos$V1)
LECIF_AtOs_Circos$V2 <- LECIF_AtOs_Circos$V2*10
LECIF_AtOs_Circos$V3 <- LECIF_AtOs_Circos$V3*10
LECIF_Os <- rbind(LECIF_AtOs_Circos, LECIF_ZmOs_Circos)
colnames(LECIF_Os) <- c("chr", "start", "end", "value1")

circos.genomicTrack(LECIF_At, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, type = "l", area = TRUE, col = colors[7], border = colors[7], ...)
})

circos.genomicTrack(LECIF_At,panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, col = colors[7], pch = 16, cex = 0.5, ...)
})

circos.genomicDensity(LECIF_Os[,1:3], col = colors[8], track.height = 0.1)
set_track_gap(gap = 0.02)

#3) LECIF Zm

LECIF_AtZm_Circos <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/At-Zm_LECIFscore.merged.bedgraph", header = F)
LECIF_OsZm_Circos <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/Os-Zm_LECIFscore.merged.bedgraph", header = F)
LECIF_AtZm_Circos$V1 <- paste0("At_",LECIF_AtZm_Circos$V1)
LECIF_OsZm_Circos$V1 <- paste0("Os_",LECIF_OsZm_Circos$V1)
LECIF_AtZm_Circos$V2 <- LECIF_AtZm_Circos$V2*10
LECIF_AtZm_Circos$V3 <- LECIF_AtZm_Circos$V3*10
LECIF_OsZm_Circos$V2 <- LECIF_OsZm_Circos$V2*7
LECIF_OsZm_Circos$V3 <- LECIF_OsZm_Circos$V3*7
LECIF_Zm <- rbind(LECIF_AtZm_Circos, LECIF_OsZm_Circos)
colnames(LECIF_Os) <- c("chr", "start", "end", "value1")

circos.genomicTrack(LECIF_At, panel.fun = function(region, value, ...) {
  circos.genomicLines(region, value, type = "l", area = TRUE, col = colors[7], border = colors[7], ...)
})

circos.genomicTrack(LECIF_At,panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, col = colors[7], pch = 16, cex = 0.5, ...)
})

circos.genomicDensity(LECIF_Zm[,1:3], col = colors[9], track.height = 0.1)
set_track_gap(gap = 0.08)

#4) hiHMM segments

At_hiHMM <- read.table("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/hihmm.model1.K15.At.ReMapped.ReNamed.Reduced.Reordered.bed", header = F)
Os_hiHMM <- read.table("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/hihmm.model1.K15.Os.ReMapped.ReNamed.Reduced.Reordered.bed", header = F)
Zm_hiHMM <- read.table("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/hihmm.model1.K15.Zm.ReMapped.ReNamed.Reduced.Reordered.bed", header = F)
At_hiHMM$V1 <- paste0("At_",At_hiHMM$V1)
Os_hiHMM$V1 <- paste0("Os_",Os_hiHMM$V1)
Zm_hiHMM$V1 <- paste0("Zm_",Zm_hiHMM$V1)
At_hiHMM$V2 <- At_hiHMM$V2*10
At_hiHMM$V3 <- At_hiHMM$V3*10
Os_hiHMM$V2 <- Os_hiHMM$V2*7
Os_hiHMM$V3 <- Os_hiHMM$V3*7
hiHMM <- rbind(At_hiHMM, Os_hiHMM, Zm_hiHMM)
colors <- hiHMM$V4
hiHMM$V4 <- gsub("CS", "", hiHMM$V4)
c("#00A087FF", "#3C5488FF", "#4DBBD5FF", "gold", "gray70")
colors <- gsub("CS2", "#3C5488FF", colors, fixed = T)
colors <- gsub("CS3", "#3C5488FF", colors, fixed = T)
colors <- gsub("CS4", "#3C5488FF", colors, fixed = T)
colors <- gsub("CS5", "#3C5488FF", colors, fixed = T)
colors <- gsub("CS6", "#00A087FF", colors, fixed = T)
colors <- gsub("CS7", "#00A087FF", colors, fixed = T)
colors <- gsub("CS8", "#00A087FF", colors, fixed = T)
colors <- gsub("CS9", "#00A087FF", colors, fixed = T)
colors <- gsub("CS10", "#4DBBD5FF", colors, fixed = T)
colors <- gsub("CS11", "#FFD700", colors, fixed = T)
colors <- gsub("CS12", "#FFD700", colors, fixed = T)
colors <- gsub("CS13", "#FFD700", colors, fixed = T)
colors <- gsub("CS14", "#666B6C", colors, fixed = T)
colors <- gsub("CS15", "#666B6C", colors, fixed = T)
colors <- gsub("CS16", "#666B6C", colors, fixed = T)
colors <- gsub("CS1", "#3C5488FF", colors, fixed = T)
colnames(hiHMM) <- c("chr", "start", "end", "value1")
hiHMM$value1 <- as.numeric(hiHMM$value1)
hiHMM_AtChr1 <- hiHMM[which(hiHMM$chr == "At_Chr1"),]
colors_AtChr1 <- colors[1:8558]
hiHMM$value1 <- 1


colorsList <- list(At_Chr1 = colors[which(hiHMM$chr == "At_Chr1")], At_Chr2 = colors[which(hiHMM$chr == "At_Chr2")],
                   At_Chr3 = colors[which(hiHMM$chr == "At_Chr3")], At_Chr4 = colors[which(hiHMM$chr == "At_Chr4")],
                   At_Chr5 = colors[which(hiHMM$chr == "At_Chr5")],
                   Os_Chr1 = colors[which(hiHMM$chr == "Os_Chr1")], Os_Chr2 = colors[which(hiHMM$chr == "Os_Chr2")],
                   Os_Chr3 = colors[which(hiHMM$chr == "Os_Chr3")], Os_Chr4 = colors[which(hiHMM$chr == "Os_Chr4")],
                   Os_Chr5 = colors[which(hiHMM$chr == "Os_Chr5")], Os_Chr6 = colors[which(hiHMM$chr == "Os_Chr6")],
                   Os_Chr7 = colors[which(hiHMM$chr == "Os_Chr7")], Os_Chr8 = colors[which(hiHMM$chr == "Os_Chr8")],
                   Os_Chr9 = colors[which(hiHMM$chr == "Os_Chr9")], Os_Chr10 = colors[which(hiHMM$chr == "Os_Chr10")],
                   Os_Chr11 = colors[which(hiHMM$chr == "Os_Chr11")], Os_Chr12 = colors[which(hiHMM$chr == "Os_Chr12")],
                   Os_ChrSy = colors[which(hiHMM$chr == "Os_ChrSy")], Os_ChrUn = colors[which(hiHMM$chr == "Os_ChrUn")],
                   Zm_1 = colors[which(hiHMM$chr == "Zm_1")], Zm_2 = colors[which(hiHMM$chr == "Zm_2")],
                   Zm_3 = colors[which(hiHMM$chr == "Zm_3")], Zm_4 = colors[which(hiHMM$chr == "Zm_4")],
                   Zm_5 = colors[which(hiHMM$chr == "Zm_5")], Zm_6 = colors[which(hiHMM$chr == "Zm_6")],
                   Zm_7 = colors[which(hiHMM$chr == "Zm_7")], Zm_8 = colors[which(hiHMM$chr == "Zm_8")],
                   Zm_9 = colors[which(hiHMM$chr == "Zm_9")], Zm_10 = colors[which(hiHMM$chr == "Zm_10")]
                   )
hiHMMList <- list(At_Chr1 = hiHMM[which(hiHMM$chr == "At_Chr1"),], At_Chr2 = hiHMM[which(hiHMM$chr == "At_Chr2"),],
                   At_Chr3 = hiHMM[which(hiHMM$chr == "At_Chr3"),], At_Chr4 = hiHMM[which(hiHMM$chr == "At_Chr4"),],
                   At_Chr5 = hiHMM[which(hiHMM$chr == "At_Chr5"),],
                   Os_Chr1 = hiHMM[which(hiHMM$chr == "Os_Chr1"),], Os_Chr2 = hiHMM[which(hiHMM$chr == "Os_Chr2"),],
                   Os_Chr3 = hiHMM[which(hiHMM$chr == "Os_Chr3"),], Os_Chr4 = hiHMM[which(hiHMM$chr == "Os_Chr4"),],
                   Os_Chr5 = hiHMM[which(hiHMM$chr == "Os_Chr5"),], Os_Chr6 = hiHMM[which(hiHMM$chr == "Os_Chr6"),],
                   Os_Chr7 = hiHMM[which(hiHMM$chr == "Os_Chr7"),], Os_Chr8 = hiHMM[which(hiHMM$chr == "Os_Chr8"),],
                   Os_Chr9 = hiHMM[which(hiHMM$chr == "Os_Chr9"),], Os_Chr10 = hiHMM[which(hiHMM$chr == "Os_Chr10"),],
                   Os_Chr11 = hiHMM[which(hiHMM$chr == "Os_Chr11"),], Os_Chr12 = hiHMM[which(hiHMM$chr == "Os_Chr12"),],
                   Os_ChrSy = hiHMM[which(hiHMM$chr == "Os_ChrSy"),], Os_ChrUn = hiHMM[which(hiHMM$chr == "Os_ChrUn"),],
                   Zm_1 = hiHMM[which(hiHMM$chr == "Zm_1"),], Zm_2 = hiHMM[which(hiHMM$chr == "Zm_2"),],
                   Zm_3 = hiHMM[which(hiHMM$chr == "Zm_3"),], Zm_4 = hiHMM[which(hiHMM$chr == "Zm_4"),],
                   Zm_5 = hiHMM[which(hiHMM$chr == "Zm_5"),], Zm_6 = hiHMM[which(hiHMM$chr == "Zm_6"),],
                   Zm_7 = hiHMM[which(hiHMM$chr == "Zm_7"),], Zm_8 = hiHMM[which(hiHMM$chr == "Zm_8"),],
                   Zm_9 = hiHMM[which(hiHMM$chr == "Zm_9"),], Zm_10 = hiHMM[which(hiHMM$chr == "Zm_10"),]
)

hiHMM$pain <- colors
#bg.border = colorsList[[v]]
circos.genomicTrack(hiHMMList, ylim=c(0,1),cell.padding=c(0,0,0,0),track.height = uh(7.5, "mm"),track.margin=c(0.001,0.001),
                    panel.fun = function(region, value, ...) {
                      i <- getI(...)
                      print(i)
                      circos.genomicLines(region, value, type = "h",col=colorsList[[i]],border=NA)
                    })
set_track_gap(gap = 0.02)

##5) Links (At,Os,Zm)
## better with list and may be merge with overlap with at leat 10000 bp | too complex
Align_AtOs <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Oryza/output/AtOs_CSsimmilarity.DF3.txt")
bed1_AtOs <- Align_AtOs[,c(1,2,3)]
bed2_AtOs <- Align_AtOs[,c(4,5,6)]
bed1_AtOs[,1] <- paste0("At_", bed1_AtOs[,1])
bed2_AtOs[,1] <- paste0("Os_", bed2_AtOs[,1])
bed1_AtOs[,2] <- bed1_AtOs[,2]*10
bed1_AtOs[,3] <- bed1_AtOs[,3]*10
bed2_AtOs[,2] <- bed2_AtOs[,2]*7
bed2_AtOs[,3] <- bed2_AtOs[,3]*7
colnames(bed1_AtOs) <- c("chr", "start", "end")
colnames(bed2_AtOs) <- c("chr", "start", "end")

Align_AtZm <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Arabidopsis_Zea/output/AtZm_CSsimmilarity.DF3.txt")
bed1_AtZm <- Align_AtZm[,c(1,2,3)]
bed2_AtZm <- Align_AtZm[,c(4,5,6)]
bed1_AtZm[,1] <- paste0("At_", bed1_AtZm[,1])
bed2_AtZm[,1] <- paste0("Zm_", bed2_AtZm[,1])
bed1_AtZm[,2] <- bed1_AtZm[,2]*10
bed1_AtZm[,3] <- bed1_AtZm[,3]*10
colnames(bed1_AtZm) <- c("chr", "start", "end")
colnames(bed2_AtZm) <- c("chr", "start", "end")

Align_OsAt <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Arabidopsis/output/OsAt_CSsimmilarity.DF3.txt")
bed1_OsAt <- Align_OsAt[,c(1,2,3)]
bed2_OsAt <- Align_OsAt[,c(4,5,6)]
bed1_OsAt[,1] <- paste0("Os_", bed1_OsAt[,1])
bed2_OsAt[,1] <- paste0("At_", bed2_OsAt[,1])
bed1_OsAt[,2] <- bed1_OsAt[,2]*7
bed1_OsAt[,3] <- bed1_OsAt[,3]*7
bed2_OsAt[,2] <- bed2_OsAt[,2]*10
bed2_OsAt[,3] <- bed2_OsAt[,3]*10
colnames(bed1_OsAt) <- c("chr", "start", "end")
colnames(bed2_OsAt) <- c("chr", "start", "end")

Align_OsZm <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Oryza_Zea/output/OsZm_CSsimmilarity.DF3.txt")
bed1_OsZm <- Align_OsZm[,c(1,2,3)]
bed2_OsZm <- Align_OsZm[,c(4,5,6)]
bed1_OsZm[,1] <- paste0("Os_", bed1_OsZm[,1])
bed2_OsZm[,1] <- paste0("Zm_", bed2_OsZm[,1])
bed1_OsZm[,2] <- bed1_OsZm[,2]*7
bed1_OsZm[,3] <- bed1_OsZm[,3]*7
colnames(bed1_OsZm) <- c("chr", "start", "end")
colnames(bed2_OsZm) <- c("chr", "start", "end")

Align_ZmAt <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis/output/ZmAt_CSsimmilarity.DF3.txt")
bed1_ZmAt <- Align_ZmAt[,c(1,2,3)]
bed2_ZmAt <- Align_ZmAt[,c(4,5,6)]
bed1_ZmAt[,1] <- paste0("Zm_", bed1_ZmAt[,1])
bed2_ZmAt[,1] <- paste0("At_", bed2_ZmAt[,1])
bed2_ZmAt[,2] <- bed2_ZmAt[,2]*10
bed2_ZmAt[,3] <- bed2_ZmAt[,3]*10
colnames(bed1_ZmAt) <- c("chr", "start", "end")
colnames(bed2_ZmAt) <- c("chr", "start", "end")

Align_ZmOs <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Oryza/output/ZmOs_CSsimmilarity.DF3.txt")
bed1_ZmOs <- Align_ZmOs[,c(1,2,3)]
bed2_ZmOs <- Align_ZmOs[,c(4,5,6)]
bed1_ZmOs[,1] <- paste0("Zm_", bed1_ZmOs[,1])
bed2_ZmOs[,1] <- paste0("Os_", bed2_ZmOs[,1])
bed2_ZmOs[,2] <- bed2_ZmOs[,2]*7
bed2_ZmOs[,3] <- bed2_ZmOs[,3]*7
colnames(bed1_ZmOs) <- c("chr", "start", "end")
colnames(bed2_ZmOs) <- c("chr", "start", "end")

colLinks <- c(rep(adjustcolor("#66C2A5", alpha.f = 1), nrow(bed1_AtOs) + nrow(bed1_AtZm)),
              rep(adjustcolor("#3288BD", alpha.f = 1), nrow(bed1_OsAt) + nrow(bed1_OsZm)),
              rep(adjustcolor("#5E4FA2", alpha.f = 1), nrow(bed1_ZmAt) + nrow(bed1_ZmOs)))

bed1 <- rbind(bed1_AtOs, bed1_AtZm, bed1_OsAt, bed1_OsZm, bed1_ZmAt, bed1_ZmOs)
bed2 <- rbind(bed2_AtOs, bed2_AtZm, bed2_OsAt, bed2_OsZm, bed2_ZmAt, bed2_ZmOs)

circos.genomicLink(bed1, bed2, col = colLinks, border = NA)


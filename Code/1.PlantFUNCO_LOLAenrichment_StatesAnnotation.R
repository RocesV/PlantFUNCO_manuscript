###################################
# TAIOTACTICAL Other Annotations  #
###################################

### Load libs

library(c("simpleCache", "LOLA", "GenomicRanges", "dplyr", "data.table", "ggplot2", "reshape2", "pheatmap", "RColorBrewer", "scales"))
library(simpleCache)
library(GenomicRanges)
library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(LOLA)
library(ggsci)

### Load regions dbs

regionDB <- loadRegionDB(dbLocation = "H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/Annotation/0.LOLA/regionDB/Zm/", collections = "TFMOTIFs")
regionDB <- loadRegionDB(dbLocation = "H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/Annotation/0.LOLA/regionDB/Zm/", collections = "Conservation")
regionDB <- loadRegionDB(dbLocation = "H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/Annotation/0.LOLA/regionDB/Zm/", collections = "SNPs")
regionDB <- loadRegionDB(dbLocation = "H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/Annotation/0.LOLA/regionDB/Zm/", collections = "ChromatinProteins")
regionDB <- loadRegionDB(dbLocation = "H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/Annotation/0.LOLA/regionDB/Zm/", collections = "HistoneMarks")

### Run enrichment

Universe <- read.table("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/bed/hihmm.model1.K15.Zm.Remapped.ReNamed.Reduced.Reordered.ALL.sorted.bed")
Universe <- Universe[,c(1,2,3)]
colnames(Universe) <- c("chr", "start", "end")
Universe <- makeGRangesFromDataFrame(Universe)

inputs <- list()
for(i in 1:length(list.files("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/bed/Zm/"))){
  input <- read.table(paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/bed/Zm/",list.files("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/bed/Zm/")[i]))
  input <- input[,c(1,2,3)]
  colnames(input) <- c("chr", "start", "end")
  inputs[[i]] <- makeGRangesFromDataFrame(input)
}
UserSets <- GRangesList(inputs)
names(UserSets) <- c("CS1", "CS10", "CS11", "CS12", "CS13", "CS14", "CS15", "CS16", "CS2", "CS3", "CS4", "CS5", "CS6", "CS7", "CS8", "CS9")

regionResults <- lapply(UserSets, FUN = function(x){
  Results <- runLOLA(x, Universe, regionDB, cores = 16)
})


### Plotting results and export

## TFMOTIFs

for(i in 1:length(regionResults)){
  name <- names(regionResults[i])
  regionResults[[i]]$file <- rep(name, nrow(regionResults[[i]]))
  regionResults[[i]]$oddsRatio[which(regionResults[[i]]$pValueLog < 1.30103)] <- 0
  regionResults[[i]]$oddsRatio[which(regionResults[[i]]$oddsRatio > 10)] <- 10
}

regionResults <- do.call(rbind, regionResults)
df <- data.frame(id = regionResults$file, variable = regionResults$description, value = regionResults$oddsRatio, condition = regionResults$treatment)
df <- reshape(df, idvar = "id", v.names = "value",  timevar = "variable", direction = "wide")
colnames(df) <- gsub("value.", "", colnames(df), fixed = T)
rownames(df) <- df$id
df <- df[,-1]
col <- colorRampPalette(brewer.pal(9, "Blues"))(250)
df <- df[,-1]
df <- t(as.matrix(df))
df <- df[,c(1,9,10,11,12,13,14,15,16,2,3,4,5,6,7,8)]
title <- "Zm TFMOTIFs"
pdf(file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/figures_model1_15/LOLA_Zm_TFMOTIFs.pdf"), paper = "a4r", height = 21, width = 28, onefile = T)
pheatmap(df, scale = "none", color = col,  cluster_rows = F, cluster_cols = F,  clustering_method = "ward.D2", cellwidth = 45, cellheight = 45,main = title, border_color = "white")
dev.off()


## Conservation 

for(i in 1:length(regionResults)){
  name <- names(regionResults[i])
  regionResults[[i]]$file <- rep(name, nrow(regionResults[[i]]))
  regionResults[[i]]$oddsRatio[which(regionResults[[i]]$pValueLog < 1.30103)] <- 0
  regionResults[[i]]$oddsRatio[which(regionResults[[i]]$oddsRatio > 10)] <- 10
}

regionResults <- do.call(rbind, regionResults)
df <- data.frame(id = regionResults$file, variable = regionResults$description, value = regionResults$oddsRatio, condition = regionResults$treatment)
df <- reshape(df, idvar = "id", v.names = "value",  timevar = "variable", direction = "wide")
colnames(df) <- gsub("value.", "", colnames(df), fixed = T)
rownames(df) <- df$id
df <- df[,-1]
col <- colorRampPalette(brewer.pal(7, "Greens"))(250)
df <- df[,-1]
df <- t(as.matrix(df))
df <- df[,c(1,9,10,11,12,13,14,15,16,2,3,4,5,6,7,8)]
title <- "Zm Conservation"
pdf(file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/figures_model1_15/LOLA_Zm_Conservation.pdf"), paper = "a4r", height = 21, width = 28, onefile = T)
pheatmap(df, scale = "none", color = col,  cluster_rows = F, cluster_cols = F,  clustering_method = "ward.D2", cellwidth = 45, cellheight = 45,main = title, border_color = "white")
dev.off()

## SNPs

for(i in 1:length(regionResults)){
  name <- names(regionResults[i])
  regionResults[[i]]$file <- rep(name, nrow(regionResults[[i]]))
  regionResults[[i]]$oddsRatio[which(regionResults[[i]]$pValueLog < 1.30103)] <- 0
  regionResults[[i]]$oddsRatio[which(regionResults[[i]]$oddsRatio > 10)] <- 10
}

regionResults <- do.call(rbind, regionResults)
df <- data.frame(id = regionResults$file, variable = regionResults$description, value = regionResults$oddsRatio, condition = regionResults$treatment)
df <- reshape(df, idvar = "id", v.names = "value",  timevar = "variable", direction = "wide")
colnames(df) <- gsub("value.", "", colnames(df), fixed = T)
rownames(df) <- df$id
df <- df[,-1]
col <- colorRampPalette(brewer.pal(7, "Greys"))(250)
names <- rownames(df)
df <- df[,-1]
df <- t(df)
colnames(df) <- names
rownames(df) <- "GWAS_SNPs"
df <- df[,c(1,9,10,11,12,13,14,15,16,2,3,4,5,6,7,8)]
title <- "Zm GWAS SNPs"
pdf(file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/figures_model1_15/LOLA_Zm_GWAS_SNPs.pdf"), paper = "a4r", height = 21, width = 28, onefile = T)
pheatmap(t(df), scale = "none", color = col,  cluster_rows = F, cluster_cols = F,  clustering_method = "ward.D2", cellwidth = 45, cellheight = 45,main = title, border_color = "white")
dev.off()

## Chromatin Proteins


for(i in 1:length(regionResults)){
  name <- names(regionResults[i])
  regionResults[[i]]$file <- rep(name, nrow(regionResults[[i]]))
  regionResults[[i]]$oddsRatio[which(regionResults[[i]]$pValueLog < 1.30103)] <- 0
  regionResults[[i]]$oddsRatio[which(regionResults[[i]]$oddsRatio > 10)] <- 10
}

regionResults <- do.call(rbind, regionResults)
df <- data.frame(id = regionResults$file, variable = regionResults$description, value = regionResults$oddsRatio, condition = regionResults$treatment)
df <- reshape(df, idvar = "id", v.names = "value",  timevar = "variable", direction = "wide")
colnames(df) <- gsub("value.", "", colnames(df), fixed = T)
rownames(df) <- df$id
df <- df[,-1]
col <- colorRampPalette(brewer.pal(9, "Purples"))(250)
names <- rownames(df)
df <- df[,-1]
df <- t(df)
colnames(df) <- names
df <- df[,c(1,9,10,11,12,13,14,15,16,2,3,4,5,6,7,8)]
title <- "Zm Chromatin Proteins"
pdf(file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/figures_model1_15/LOLA_Zm_OtherChromatinProteins.pdf"), paper = "a4r", height = 28, width = 28, onefile = T)
pheatmap(df, scale = "none", color = col,  cluster_rows = T, cluster_cols = F,  clustering_method = "ward.D2", cellwidth = 18, cellheight = 12,main = title, border_color = "white" )
dev.off()


## Other Histone marks

for(i in 1:length(regionResults)){
  name <- names(regionResults[i])
  regionResults[[i]]$file <- rep(name, nrow(regionResults[[i]]))
  regionResults[[i]]$oddsRatio[which(regionResults[[i]]$pValueLog < 1.30103)] <- 0
  regionResults[[i]]$oddsRatio[which(regionResults[[i]]$oddsRatio > 10)] <- 10
}
regionResults <- do.call(rbind, regionResults)
df <- data.frame(id = regionResults$file, condition = regionResults$treatment, variable =  regionResults$description,value = regionResults$oddsRatio)
Conditions <- data.frame(Condition = df$condition[which(df$id == "CS1")])
rownames(Conditions) <- df$variable[which(df$id == "CS1")]
df <- reshape(df[,-2], idvar = "id", v.names = "value",  timevar = "variable", direction = "wide")
df[is.na(df)] <- 1
rownames(df) <- df$id
df <- df[,-1]
col <- colorRampPalette(brewer.pal(9, "Reds"))(250)
colnames(df) <- gsub("value.", "", colnames(df), fixed = T)
df <- t(as.matrix(df))
title <- "Zm Histone Marks"
Annotation_colors <- list(Conditions = c(active = "darkgreen", repressive = "darksalmon", Nucleosomes = "lightsteelblue1", notclear = "grey"))
df <- df[,c(1,9,10,11,12,13,14,15,16,2,3,4,5,6,7,8)]
pdf(file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/figures_model1_15/LOLA_Os_HistoneMarks.pdf"), paper = "a4r", height = 21, width = 28, onefile = T)
pheatmap(df, scale = "none", color = col,  cluster_rows = F, cluster_cols = F,  clustering_method = "ward.D2", cellwidth = 35, cellheight = 35,main = title, border_color = "white", annotation_colors = Annotation_colors, annotation_row = Conditions)
dev.off()

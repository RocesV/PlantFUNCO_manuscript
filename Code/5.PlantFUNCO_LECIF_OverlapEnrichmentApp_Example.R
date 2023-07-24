
regionDB <- list("Arabidopsis" = list(), "Oryza" = list(), "Zea" = list())

regionDB$Zea$PhyloP <- loadRegionDB(dbLocation = "./regionDB/Zea/", collections = "PhyloP")

At.PhyloP <- read.delim("H:/Taiotactical/Taiotactical_LECIF/0.Zea_Arabidopsis//output/Zma_PhyloP_merged.bedGraph", sep = "\t", header = F)
At.PhyloP1 <- At.PhyloP[which(percent_rank(x = At.PhyloP$V4) >= 0 & percent_rank(x = At.PhyloP$V4) < 0.2),]
At.PhyloP2 <- At.PhyloP[which(percent_rank(x = At.PhyloP$V4) >= 0.2 & percent_rank(x = At.PhyloP$V4) < 0.4),]
At.PhyloP3 <- At.PhyloP[which(percent_rank(x = At.PhyloP$V4) >= 0.4 & percent_rank(x = At.PhyloP$V4) < 0.6),]
At.PhyloP4 <- At.PhyloP[which(percent_rank(x = At.PhyloP$V4) >= 0.6 & percent_rank(x = At.PhyloP$V4) < 0.8),]
At.PhyloP5 <- At.PhyloP[which(percent_rank(x = At.PhyloP$V4) >= 0.8 & percent_rank(x = At.PhyloP$V4) <= 1),]

write.table(At.PhyloP1, "H:/Taiotactical/Taiotactical_LECIF/3.Applications/OverlapEnrichment/regionDB/Zea/PhyloP/regions/PhyloP_Bin1.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(At.PhyloP2, "H:/Taiotactical/Taiotactical_LECIF/3.Applications/OverlapEnrichment/regionDB/Zea/PhyloP/regions/PhyloP_Bin2.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(At.PhyloP3, "H:/Taiotactical/Taiotactical_LECIF/3.Applications/OverlapEnrichment/regionDB/Zea/PhyloP/regions/PhyloP_Bin3.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(At.PhyloP4, "H:/Taiotactical/Taiotactical_LECIF/3.Applications/OverlapEnrichment/regionDB/Zea/PhyloP/regions/PhyloP_Bin4.bed", sep = "\t", col.names = F, row.names = F, quote = F)
write.table(At.PhyloP5, "H:/Taiotactical/Taiotactical_LECIF/3.Applications/OverlapEnrichment/regionDB/Zea/PhyloP/regions/PhyloP_Bin5.bed", sep = "\t", col.names = F, row.names = F, quote = F)

At <- makeTxDbFromGFF(paste0("H:/Taiotactical/Taiotactical_ChIPseq/Utils/Annotation/0.At/Annotation/0.At.gene_exons.gff3"))
Genes.At <- genes(At)
Genes.At.df <- data.frame(Chr = as.character(seqnames(Genes.At)), Start = Genes.At@ranges@start, End = Genes.At@ranges@start + Genes.At@ranges@width - 1, GeneID = Genes.At$gene_id)
write.table(Genes.At.df[,c(1,2,3)], "H:/Taiotactical/Taiotactical_LECIF/3.Applications/OverlapEnrichment/example/ArabidopsisGenome.background.bed", quote = F, sep = "\t", row.names = F, col.names = F)

DiffGenes.Abiotic <- read.delim("clipboard", header = T)
length(which((Genes.At.df$GeneID %in% DiffGenes.Abiotic$GENE_ID) == TRUE))
DiffGenes.Abiotic <- Genes.At.df[Genes.At.df$GeneID %in% DiffGenes.Abiotic$GENE_ID,]
write.table(DiffGenes.Abiotic[,c(1,2,3)], "H:/Taiotactical/Taiotactical_LECIF/3.Applications/OverlapEnrichment/example/ArabidopsisDiffGenes.abiotic.bed", quote = F, sep = "\t", row.names = F, col.names = F)

DiffGenes.Biotic <- read.delim("clipboard", header = T)
length(which((Genes.At.df$GeneID %in% DiffGenes.Biotic$GENE_ID) == TRUE))
DiffGenes.Biotic <- Genes.At.df[Genes.At.df$GeneID %in% DiffGenes.Biotic$GENE_ID,]
write.table(DiffGenes.Biotic[,c(1,2,3)], "H:/Taiotactical/Taiotactical_LECIF/3.Applications/OverlapEnrichment/example/ArabidopsisDiffGenes.biotic.bed", quote = F, sep = "\t", row.names = F, col.names = F)

DiffGenes.Tissues <- read.delim("clipboard", header = T)
length(which((Genes.At.df$GeneID %in% DiffGenes.Tissues$GENE_ID) == TRUE))
DiffGenes.Tissues <- Genes.At.df[Genes.At.df$GeneID %in% DiffGenes.Tissues$GENE_ID,]
write.table(DiffGenes.Tissues[,c(1,2,3)], "H:/Taiotactical/Taiotactical_LECIF/3.Applications/OverlapEnrichment/example/ArabidopsisDiffGenes.tissues.bed", quote = F, sep = "\t", row.names = F, col.names = F)

DiffSplicing.Abiotic <- read.delim("clipboard", header = T)
DiffSplicing.Abiotic$Chr <- do.call(rbind,strsplit(DiffSplicing.Abiotic$COORD, split = ":", fixed = T))[,1]
DiffSplicing.Abiotic$Chr <- gsub("chr", "Chr", fixed = T, x = DiffSplicing.Abiotic$Chr)
DiffSplicing.Abiotic$COORD <- do.call(rbind,strsplit(DiffSplicing.Abiotic$COORD, split = ":", fixed = T))[,2]
DiffSplicing.Abiotic$Start <- do.call(rbind,strsplit(DiffSplicing.Abiotic$COORD, split = "-", fixed = T))[,1]
DiffSplicing.Abiotic$End <- do.call(rbind,strsplit(DiffSplicing.Abiotic$COORD, split = "-", fixed = T))[,2]
DiffSplicing.Abiotic$Start <- as.numeric(DiffSplicing.Abiotic$Start)
DiffSplicing.Abiotic$End <- as.numeric(DiffSplicing.Abiotic$End)
DiffSplicing.Abiotic$Start[is.na(DiffSplicing.Abiotic$Start)] <- DiffSplicing.Abiotic$End[is.na(DiffSplicing.Abiotic$Start)] - 1
DiffSplicing.Abiotic$End[is.na(DiffSplicing.Abiotic$End)] <- DiffSplicing.Abiotic$End[is.na(DiffSplicing.Abiotic$End)] + 1
write.table(DiffSplicing.Abiotic[,c(2,3,4)], "H:/Taiotactical/Taiotactical_LECIF/3.Applications/OverlapEnrichment/example/ArabidopsisDiffSplicing.abiotic.bed", quote = F, sep = "\t", row.names = F, col.names = F)

DiffSplicing.Biotic <- read.delim("clipboard", header = T)
DiffSplicing.Biotic$Chr <- do.call(rbind,strsplit(DiffSplicing.Biotic$COORD, split = ":", fixed = T))[,1]
DiffSplicing.Biotic$Chr <- gsub("chr", "Chr", fixed = T, x = DiffSplicing.Biotic$Chr)
DiffSplicing.Biotic$COORD <- do.call(rbind,strsplit(DiffSplicing.Biotic$COORD, split = ":", fixed = T))[,2]
DiffSplicing.Biotic$Start <- do.call(rbind,strsplit(DiffSplicing.Biotic$COORD, split = "-", fixed = T))[,1]
DiffSplicing.Biotic$End <- do.call(rbind,strsplit(DiffSplicing.Biotic$COORD, split = "-", fixed = T))[,2]
DiffSplicing.Biotic$Start <- as.numeric(DiffSplicing.Biotic$Start)
DiffSplicing.Biotic$End <- as.numeric(DiffSplicing.Biotic$End)
DiffSplicing.Biotic$Start[is.na(DiffSplicing.Biotic$Start)] <- DiffSplicing.Biotic$End[is.na(DiffSplicing.Biotic$Start)] - 1
DiffSplicing.Biotic$End[is.na(DiffSplicing.Biotic$End)] <- DiffSplicing.Biotic$End[is.na(DiffSplicing.Biotic$End)] + 1
write.table(DiffSplicing.Biotic[,c(2,3,4)], "H:/Taiotactical/Taiotactical_LECIF/3.Applications/OverlapEnrichment/example/ArabidopsisDiffSplicing.biotic.bed", quote = F, sep = "\t", row.names = F, col.names = F)

DiffSplicing.Tissues <- read.delim("clipboard", header = T)
DiffSplicing.Tissues$Chr <- do.call(rbind,strsplit(DiffSplicing.Tissues$COORD, split = ":", fixed = T))[,1]
DiffSplicing.Tissues$Chr <- gsub("chr", "Chr", fixed = T, x = DiffSplicing.Tissues$Chr)
DiffSplicing.Tissues$COORD <- do.call(rbind,strsplit(DiffSplicing.Tissues$COORD, split = ":", fixed = T))[,2]
DiffSplicing.Tissues$Start <- do.call(rbind,strsplit(DiffSplicing.Tissues$COORD, split = "-", fixed = T))[,1]
DiffSplicing.Tissues$End <- do.call(rbind,strsplit(DiffSplicing.Tissues$COORD, split = "-", fixed = T))[,2]
DiffSplicing.Tissues$Start <- as.numeric(DiffSplicing.Tissues$Start)
DiffSplicing.Tissues$End <- as.numeric(DiffSplicing.Tissues$End)
DiffSplicing.Tissues$Start[is.na(DiffSplicing.Tissues$Start)] <- DiffSplicing.Tissues$End[is.na(DiffSplicing.Tissues$Start)] - 1
DiffSplicing.Tissues$End[is.na(DiffSplicing.Tissues$End)] <- DiffSplicing.Tissues$End[is.na(DiffSplicing.Tissues$End)] + 1
write.table(DiffSplicing.Tissues[,c(2,3,4)], "H:/Taiotactical/Taiotactical_LECIF/3.Applications/OverlapEnrichment/example/ArabidopsisDiffSplicing.tissues.bed", quote = F, sep = "\t", row.names = F, col.names = F)



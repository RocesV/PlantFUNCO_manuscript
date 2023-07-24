###################################
# TAIOTACTICAL-Annotation seeker #
###################################

## 0. Load libs and pkgs
library(GenomicFeatures)


## 2 COORDS
species <- c("0.At","0.Os", "0.Zm")
OutputDir <- c("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/Annotation/")
for(i in species){
  
  Nuclear <- makeTxDbFromGFF(paste0("H:/Taiotactical/Taiotactical_ChIPseq/Utils/Annotation/", i, "/Annotation/", i, ".gene_exons.gff3"))
  
  Exons.bed <- list()
  Introns.bed <- list()
  genes.bed <- list()
  Intergenic.bed <- list()
  III.UTR.bed <- list()
  V.UTR.bed <- list()
  TSS.bed <- list()
  TSS.2kb.bed <- list()
  TES.bed <- list()
  
  
  # Exons
  Exons <- exonicParts(Nuclear, linked.to.single.gene.only = F)
  Exons.bed[["Nuclear"]] <- data.frame(chr = as.character(Exons@seqnames), start = Exons@ranges@start, end = (Exons@ranges@start + Exons@ranges@width - 1))
  write.table(Exons.bed[["Nuclear"]], file = paste0(OutputDir,"/PlantinaExon.",i,".bed"), col.names = F, row.names = F, sep = "\t", 
              quote = F)
  # Introns
  Introns <- intronicParts(Nuclear)
  Introns.bed[['Nuclear']] <- data.frame(chr = as.character(Introns@seqnames), start = Introns@ranges@start, end = (Introns@ranges@start + Introns@ranges@width - 1))
  write.table(Introns.bed[["Nuclear"]], file = paste0(OutputDir,"/PlantinaIntron.",i,".bed"), col.names = F, row.names = F, sep = "\t", 
              quote = F)
  
  # Genes and Intergenic
  Genes <- genes(Nuclear)
  Intergenic <- gaps(reduce(Genes, ignore.strand=T))
  Intergenic <- Intergenic[strand(Intergenic) == "*"]
  genes.bed[['Nuclear']] <- data.frame(chr = as.character(Genes@seqnames), start = Genes@ranges@start, end = (Genes@ranges@start + Genes@ranges@width - 1))
  genes.bed$Nuclear$chr <- factor(genes.bed$Nuclear$chr, levels(as.factor(genes.bed$Nuclear$chr)), ordered = TRUE)
  genes.bed$Nuclear <- genes.bed$Nuclear[order(genes.bed$Nuclear$chr, genes.bed$Nuclear$start),]
  Intergenic.bed[['Nuclear']] <- data.frame(chr = as.character(Intergenic@seqnames), start = Intergenic@ranges@start, end = (Intergenic@ranges@start + Intergenic@ranges@width - 1))
  Intergenic.bed$Nuclear$chr <- factor(Intergenic.bed$Nuclear$chr, levels(as.factor(Intergenic.bed$Nuclear$chr)), ordered = TRUE)
  Intergenic.bed$Nuclear <- Intergenic.bed$Nuclear[order(Intergenic.bed$Nuclear$chr, Intergenic.bed$Nuclear$start),]
  write.table(genes.bed$Nuclear, file = paste0(OutputDir,"/PlantinaGene.",i,".bed"),col.names = F, row.names = F, sep = "\t", 
              quote = F)
  write.table(Intergenic.bed$Nuclear, file = paste0(OutputDir,"/PlantinaIntergenic.",i,".bed"),col.names = F, row.names = F, sep = "\t", 
              quote = F)
  
  # UTRs/TSS/TSS2kb/TES/TSS Anchor/TES Anchoer
  
  III.UTR <- unlist(threeUTRsByTranscript(Nuclear), use.names = F)
  III.UTR.bed[['Nuclear']] <- data.frame(chr = as.character(III.UTR@seqnames), start = III.UTR@ranges@start, end = (III.UTR@ranges@start + III.UTR@ranges@width - 1))
  III.UTR.bed$Nuclear <- III.UTR.bed$Nuclear[!duplicated(III.UTR.bed$Nuclear),]
  III.UTR.bed$Nuclear$chr <- factor(III.UTR.bed$Nuclear$chr, levels(as.factor(III.UTR.bed$Nuclear$chr)), ordered = TRUE)
  III.UTR.bed$Nuclear <- III.UTR.bed$Nuclear[order(III.UTR.bed$Nuclear$chr, III.UTR.bed$Nuclear$start),]
  write.table(III.UTR.bed$Nuclear, file = paste0(OutputDir,"/Plantina3UTR.",i,".bed"), col.names = F, row.names = F, sep = "\t", 
              quote = F)
  
  V.UTR <- unlist(fiveUTRsByTranscript(Nuclear), use.names = F)
  V.UTR.bed[['Nuclear']] <- data.frame(chr = as.character(V.UTR@seqnames), start = V.UTR@ranges@start, end = (V.UTR@ranges@start + V.UTR@ranges@width - 1))
  V.UTR.bed$Nuclear <- V.UTR.bed$Nuclear[!duplicated(V.UTR.bed$Nuclear),]
  V.UTR.bed$Nuclear$chr <- factor(V.UTR.bed$Nuclear$chr, levels(as.factor(V.UTR.bed$Nuclear$chr)), ordered = TRUE)
  V.UTR.bed$Nuclear <- V.UTR.bed$Nuclear[order(V.UTR.bed$Nuclear$chr, V.UTR.bed$Nuclear$start),]
  write.table(V.UTR.bed$Nuclear, file = paste0(OutputDir,"/Plantina5UTR.",i,".bed"), col.names = F, row.names = F, sep = "\t", 
              quote = F)
  
  # TES
  
  V.UTR <- fiveUTRsByTranscript(Nuclear)
  III.UTR <- threeUTRsByTranscript(Nuclear)
  
  TES.chr <- list()
  TES.start <- list()
  TES.end <- list()
  TES.strand <- list()
  for(j in 1:length(III.UTR)){
    TES <- III.UTR[j]
    TES.chr[[j]] <- as.character(TES@unlistData@seqnames@values)
    TES.strand[[j]] <- TES@unlistData@strand@values
    if(TES@unlistData@strand@values != "-"){
      TES.end[[j]] <- max(TES@unlistData@ranges@start + TES@unlistData@ranges@width - 1)
      TES.start[[j]] <- TES.end[[j]] - 1
    }else if(TES@unlistData@strand@values == "-"){
      TES.start[[j]] <- min(TES@unlistData@ranges@start + TES@unlistData@ranges@width - 1)
      TES.end[[j]] <- TES.start[[j]] + 1
    }
  }
  
  TES.COORD <- data.frame(chr = unlist(TES.chr), start = unlist(TES.start), end = unlist(TES.end), strad = unlist(TES.strand))
  TES.COORD <-  TES.COORD[!duplicated(TES.COORD),]
  TES.COORD$chr <- factor(TES.COORD$chr, levels(as.factor(TES.COORD$chr)), ordered = TRUE)
  TES.COORD <- TES.COORD[order(TES.COORD$chr, TES.COORD$start),]
  write.table(TES.COORD[,1:3], file = paste0(OutputDir,"/PlantinaTES.",i,".bed"), col.names = F, row.names = F, sep = "\t", 
              quote = F)
  
  TES.ANCHOR.start <- rep(NA, nrow(TES.COORD)) 
  TES.ANCHOR.start[which(TES.COORD$strad != "-")] <- TES.COORD$end[which(TES.COORD$strad != "-")]
  TES.ANCHOR.start[which(TES.COORD$strad == "-")] <- TES.COORD$start[which(TES.COORD$strad == "-")]
  TES.ANCHOR <- data.frame(chr = TES.COORD$chr, start = TES.ANCHOR.start, strand = TES.COORD$strad)
  TES.ANCHOR <- TES.ANCHOR[order(TES.ANCHOR$chr, TES.ANCHOR$start),]
  write.table(TES.ANCHOR, file = paste0(OutputDir,"/ANCHOR/PlantinaTES.",i,".bed"), col.names = F, row.names = F, sep = "\t", 
              quote = F)
  
  # TSS
  
  TSS.chr <- list()
  TSS.start <- list()
  TSS.end <- list()
  TSS.strand <- list()
  for(j in 1:length(V.UTR)){
    TSS <- V.UTR[j]
    TSS.chr[[j]] <- as.character(TSS@unlistData@seqnames@values)
    TSS.strand[[j]] <- TSS@unlistData@strand@values
    if(TSS@unlistData@strand@values != "-"){
      TSS.start[[j]] <- min(TSS@unlistData@ranges@start + TSS@unlistData@ranges@width - 1)
      TSS.end[[j]] <- TSS.start[[j]] + 1
    }else if(TSS@unlistData@strand@values == "-"){
      TSS.end[[j]] <- max(TSS@unlistData@ranges@start + TSS@unlistData@ranges@width - 1)
      TSS.start[[j]] <- TSS.end[[j]] - 1
    }
  }
  
  TSS.COORD <- data.frame(chr = unlist(TSS.chr), start = unlist(TSS.start), end = unlist(TSS.end), strand = unlist(TSS.strand))
  TSS.COORD <-  TSS.COORD[!duplicated(TSS.COORD),]
  TSS.COORD$chr <- factor(TSS.COORD$chr, levels(as.factor(TSS.COORD$chr)), ordered = TRUE)
  TSS.COORD <- TSS.COORD[order(TSS.COORD$chr, TSS.COORD$start),]
  write.table(TSS.COORD[,1:3], file = paste0(OutputDir,"/PlantinaTSS.",i,".bed"), col.names = F, row.names = F, sep = "\t", 
              quote = F)
  
  TSS.ANCHOR.start <- rep(NA, nrow(TSS.COORD)) 
  TSS.ANCHOR.start[which(TSS.COORD$strand != "-")] <- TSS.COORD$start[which(TSS.COORD$strand != "-")]
  TSS.ANCHOR.start[which(TSS.COORD$strand == "-")] <- TSS.COORD$end[which(TSS.COORD$strand == "-")]
  TSS.ANCHOR <- data.frame(chr = TSS.COORD$chr, start = TSS.ANCHOR.start, strand = TSS.COORD$strand)
  TSS.ANCHOR <- TSS.ANCHOR[order(TSS.ANCHOR$chr, TSS.ANCHOR$start),]
  write.table(TSS.ANCHOR, file = paste0(OutputDir,"/ANCHOR/PlantinaTSS.",i,".bed"), col.names = F, row.names = F, sep = "\t", 
              quote = F)
  
  TSS.2kb.COORD <- TSS.COORD
  TSS.2kb.COORD$start <- TSS.2kb.COORD$start - 1000
  TSS.2kb.COORD$end <- TSS.2kb.COORD$end + 1000
  write.table(TSS.2kb.COORD[,1:3], file = paste0(OutputDir,"/PlantinaTSS2Kb.",i,".bed"), col.names = F, row.names = F, sep = "\t", 
              quote = F)
  
  # CDS
  CDS.bed <- list()
  CDS <- unlist(cdsBy(Nuclear, by = "gene"),use.names = F)
  CDS.bed[["Nuclear"]] <- data.frame(chr = as.character(CDS@seqnames), start = CDS@ranges@start, end = (CDS@ranges@start + CDS@ranges@width - 1))
  CDS.bed$Nuclear <- CDS.bed$Nuclear[!duplicated(CDS.bed$Nuclear),]
  CDS.bed$Nuclear$chr <- factor(CDS.bed$Nuclear$chr, levels(as.factor(CDS.bed$Nuclear$chr)), ordered = TRUE)
  CDS.bed$Nuclear <- CDS.bed$Nuclear[order(CDS.bed$Nuclear$chr, CDS.bed$Nuclear$start),]
  write.table(CDS.bed[["Nuclear"]], file = paste0(OutputDir,"/PlantinaCDS.",i,".bed"), col.names = F, row.names = F, sep = "\t", 
              quote = F)
  
  # Promoter
  Promoters.bed <- list()
  Promoters <- promoters(x = Nuclear, upstream = 1000, downstream = 1000, use.names = F)
  Promoters.bed[["Nuclear"]] <- data.frame(chr = as.character(Promoters@seqnames), start = Promoters@ranges@start, end = (Promoters@ranges@start + Promoters@ranges@width - 1))
  Promoters.bed$Nuclear <- Promoters.bed$Nuclear[!duplicated(Promoters.bed$Nuclear),]
  Promoters.bed$Nuclear$chr <- factor(Promoters.bed$Nuclear$chr, levels(as.factor(Promoters.bed$Nuclear$chr)), ordered = TRUE)
  Promoters.bed$Nuclear <- Promoters.bed$Nuclear[order(Promoters.bed$Nuclear$chr, Promoters.bed$Nuclear$start),]
  write.table(Promoters.bed[["Nuclear"]], file = paste0(OutputDir,"/PlantinaPromoter.",i,".bed"), col.names = F, row.names = F, sep = "\t", 
              quote = F)
  
  
}

## annotation by coords and anchor using chromHMM functionality

hiHMM_At <- read.table("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/hihmm.model1.K15.Zm.ReMapped.bed", header = F)
hiHMM_At$V4[which(hiHMM_At$V4 == "0")] <- "22"
hiHMM_At$V4 <- paste0("E", hiHMM_At$V4)
write.table(hiHMM_At, "H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/hihmm.model1.K15.Zm.ReMapped.ReNamed.bed", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(hiHMM_At[,1:4], "H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/hihmm.model1.K15.Zm.ReMapped.ReNamed.Reduced.bed", quote = F, row.names = F, col.names = F, sep = "\t")



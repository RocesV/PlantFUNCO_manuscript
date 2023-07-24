# CS to gff3

CStoGFF3 <- function(CS){
  seqlevels <- levels(as.factor(CS$V1))
  seqlevels.list <- list()
  for(i in seqlevels){
    intermediate <- CS[which(CS$V1 == i),]
    seqlevels.df <- data.frame(Chr = rep(i, nrow(intermediate)+1),
                               DB = rep("PlantFUNCO", nrow(intermediate)+1),
                               Element = c("gene", intermediate$V4),
                               Start = c(min(intermediate$V2), intermediate$V2),
                               End = c(max(intermediate$V3), intermediate$V3),
                               score = rep(".", nrow(intermediate)+1),
                               strand = rep("+", nrow(intermediate)+1),
                               codon = rep(".", nrow(intermediate)+1),
                               Tag = c(paste0("ID=CS.At.", i), rep(paste0("ID=",intermediate$V4,".",rownames(intermediate),";Parent=CS.At.", i))))
    seqlevels.list[[i]] <- seqlevels.df
  }
  gff3 <- do.call(rbind, seqlevels.list)
 return(gff3)
}

## commands 

At.CS <- read.delim("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/hihmm.model1.K15.Os.ReMapped.ReNamed.Reduced.Reordered.bed", header = F, sep = "\t")
At.GFF3 <- CStoGFF3(At.CS)
write.table(At.GFF3, "H:/Taiotactical/Taiotactical_LECIF/3.Applications/Browser/jbrowse2/PlantFUNCO_data/Oryza/1.CShiHMM/CShiHMM.gff3", quote = F, sep = "\t", col.names = F, row.names = F)

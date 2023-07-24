setwd("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/")
Species <- c("At", "Os", "Zm")

for(i in Species){
  bed <- read.table(paste0("hihmm.model1.K15.",i,".ReMapped.ReNamed.Reduced.bed"), header = F)
  States <- c("E22", "E15",  "E14",  "E13",  "E12",  "E11",  "E10",  "E9",  "E8",  "E7",  "E6",  "E5",  "E4", "E3", "E2", "E1")
  CS <- c("CS16", "CS7",  "CS13",  "CS9",  "CS8",  "CS5",  "CS10",  "CS2",  "CS6",  "CS15",  "CS14",  "CS12",  "CS11", "CS1", "CS3", "CS4")
  Conversion <- data.frame(States = States, CS = CS)
  bed.state <- list()
  bed.state.fasta <- list() 
  for(j in States){
    bed.state[[j]] <- bed[which(bed$V4 == j),]
    bed.state[[j]]$V4 <- Conversion$CS[which(Conversion$States == j)]
    bed.state.fasta[[j]] <- bed.state[[j]]
    bed.state.fasta[[j]]$V4 <- paste0(bed.state.fasta[[j]]$V4, "_", 1:nrow(bed.state.fasta[[j]]))
  }
  for(z in 1:length(bed.state.fasta)){
    write.table(x = bed.state.fasta[[z]],paste0("./splitted/bed/hihmm.model1.K15.", i,".ReMapped.ReNamed.Reduced.Reordered.",names(bed.state.fasta[z]),".bed"),sep = "\t", quote = F, row.names = F, col.names = F)
  }
  ALL <- do.call(rbind, bed.state)
  ALL <- ALL[order(ALL$V1, ALL$V2),]
  write.table(x = ALL, paste0("hihmm.model1.K15.",i,".ReMapped.ReNamed.Reduced.Reordered.bed"), quote = F, sep = "\t", row.names = F, col.names = F)
}

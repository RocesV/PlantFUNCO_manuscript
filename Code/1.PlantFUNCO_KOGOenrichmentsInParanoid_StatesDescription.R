###################################
# TAIOTACTICAL BED description   #
###################################


### Load libs

library(ggalluvial)
library(GenomicFeatures)
library(GenomicRanges)
library(stringr)
library(clusterProfiler)

test <- as.data.frame(UCBAdmissions)

### Extract alluvial data for each Species: Araport type of genes

At.df <- data.frame(Type = NA, Orthology = NA, Relation = NA, Freq = NA)
CS.list <- list()
CS.genes <- list()
Gene_info <- read.table("H:/Taiotactical/Taiotactical_ChIPseq/Utils/Annotation/0.At/Annotation/Araport11_gene_type.txt", header = T)
inputDir <- "H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/genes/At/"

for(i in 1:length(list.files(inputDir))){
  CS.list[[i]] <- read.table(paste0(inputDir, list.files(inputDir)[i]), header = F, sep = "\t")
  CS.list[[i]]$V13 <- gsub(pattern = ";", replacement = "", do.call(rbind, strsplit(CS.list[[i]]$V13, split = " ", fixed = T))[,4])
  CS.genes[[i]] <- CS.list[[i]]$V13
  CS.genes[[i]] <- data.frame(Genes = CS.genes[[i]][!duplicated(CS.genes[[i]])], Type = NA)
  CS.genes[[i]] <- CS.genes[[i]][-which((CS.genes[[i]]$Genes %in% Gene_info$name) == FALSE),]
  Iter.rows <- which((CS.genes[[i]]$Genes %in% Gene_info$name) == TRUE)
  for(j in Iter.rows){
    CS.genes[[i]]$Type[j] <- Gene_info$gene_model_type[which(Gene_info$name == CS.genes[[i]]$Genes[j])]
  }
}

# now it is orthology time: Phytozome ones

CS.genes <- lapply(CS.genes, function(x){
  x$Orthology <- NA
  x$Species <- NA
  x
})

# Simplier is the best - Orthologs Os, Orthologs Zm, Orthologs both, non-othologs
At.orthology.Zm <- read.table("H:/Taiotactical/Taiotactical_ChIPseq/Utils/Annotation/0.At/Orthology/inParanoid_Zmays_493_RefGen_V4_Athaliana_447_Araport11", header = F, sep = "\t")
At.orthology.Os <- read.table("H:/Taiotactical/Taiotactical_ChIPseq/Utils/Annotation/0.At/Orthology/inParanoid_Osativa_323_v7.0_Athaliana_447_Araport11", header = F, sep = "\t")
# Format this
At.orthology.Zm$V2 <- gsub(pattern = " ", replacement = ".",At.orthology.Zm$V2)
At.orthology.Zm$V2 <- gsub(pattern = ".", replacement = "_",At.orthology.Zm$V2, fixed = T)
for(i in 1:nrow(At.orthology.Zm)){
  V1 <- strsplit(At.orthology.Zm$V1[i], split = "_", fixed = T)[[1]][grep(x =  strsplit(At.orthology.Zm$V1[i], split = "_", fixed = T)[[1]], "Zmays:", fixed = T)]
  V2 <- strsplit(At.orthology.Zm$V2[i], split = "_", fixed = T)[[1]][grep(x =  strsplit(At.orthology.Zm$V2[i], split = "_", fixed = T)[[1]], "Athalianacolumbia:", fixed = T)]
  V1 <- strsplit(x = gsub(pattern = "Zmays:", replacement = "", x = V1), split = "_", fixed = T)[[1]]
  V2 <- gsub(pattern = "Athalianacolumbia:", replacement = "", x = V2)
  
}



# check Os and Zm genes and if both check one2one both one2one Os one2one Zm ...

for(i in 1:length(CS.genes)){
  for(j in 1:nrow(CS.genes[[i]])){
    # hit in Os
    if(length(grep(CS.genes[[i]]$Genes[j], x = At.orthology.Os$V2)) > 0){
      Orthology.Os <- "orthologous"
    }else if(length(grep(CS.genes[[i]]$Genes[j], x = At.orthology.Os$V2)) == 0){
      Orthology.Os <- "non_orthologous"
    }
    # hit in Zm
    if(length(grep(CS.genes[[i]]$Genes[j], x = At.orthology.Zm$V2)) > 0){
      Orthology.Zm <- "orthologous"
    }else if(length(grep(CS.genes[[i]]$Genes[j], x = At.orthology.Zm$V2)) == 0){
      Orthology.Zm <- "non_orthologous"
    }    
    # coherency in both 
    if(Orthology.Os == Orthology.Zm & Orthology.Os == "orthologous"){
      Species <- "both"
      Orthology <- "orthologous"
    }else if(Orthology.Os == Orthology.Zm & Orthology.Os == "non_orthologous"){
      Species <- "both"
      Orthology <- "non_orthologous"
    }else if(Orthology.Os != Orthology.Zm & Orthology.Os == "orthologous"){
      Species <- "Os"
      Orthology <- "orthologous"
    }else if(Orthology.Os != Orthology.Zm & Orthology.Zm == "orthologous"){
      Species <- "Zm"
      Orthology <- "orthologous"
    }
    # fill data.frame for each state
    CS.genes[[i]]$Orthology[j] <- Orthology
    CS.genes[[i]]$Species[j] <- Species
  }
}


## Plot loop for each state accompanied in the title by total number of genes and prot coding genes: orthologous non-orthologous
## export tables just in case

CS.names <- c("CS1","CS10", "CS11", "CS12", "CS13", "CS14", "CS15", "CS16", "CS2", 
              "CS3", "CS4", "CS5", "CS6", "CS7", "CS8", "CS9")

for(i in 1:length(CS.genes)){
  # export tables
  #write.table(as.data.frame(table(CS.genes[[i]][,2:4])), paste0(inputDir, CS.names[i], "_At_genestypeTable.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
  #write.table(CS.genes[[i]], paste0(inputDir, CS.names[i], "_At_allgenes.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
  # plot alluvial
  allplot <- ggplot(as.data.frame(table(CS.genes[[i]][,c(2:4)])),
         aes(y = Freq,
             axis1 = Type, axis2 = Orthology)) +
    geom_alluvium(aes(fill = Species),
                  width = 1/8)  +
    geom_stratum(width = 1/8, fill = "white") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("Type", "Orthology"), expand = c(.05, .05)) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    ggtitle(paste0(CS.names[i], "Total Genes: ", nrow(CS.genes[[i]]), " Prot. coding genes: ", length(which(CS.genes[[i]]$Type == "protein_coding")), " Orthologous: ", length(which(CS.genes[[i]]$Orthology == "orthologous")))) + theme_classic()
  pdf(file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/figures_model1_15/", CS.names[i],"_alluvial_At.pdf"), paper = "a4r", height = 28, width = 28, onefile = T)
  allplot
  dev.off()
  }

################################################################## Oryza sativa
### Extract alluvial data for each Species: Oryza not classified genes jut protein_coding and TEs

Os.df <- data.frame(Type = NA, Orthology = NA, Relation = NA, Freq = NA)
CS.list <- list()
CS.genes <- list()
Gene_info <- read.delim("clipboard", header = T)
inputDir <- "H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/genes/Os/"

for(i in 1:length(list.files(inputDir))){
  CS.list[[i]] <- read.table(paste0(inputDir, list.files(inputDir)[i]), header = F, sep = "\t")
  for(j in 1:nrow(CS.list[[i]])){
    if(CS.list[[i]]$V1[j] != "ChrSy" & CS.list[[i]]$V1[j] != "ChrUn"){
      CS.list[[i]]$V13[j] <- gsub(pattern = "ID=", replacement = "", strsplit(strsplit(CS.list[[i]]$V13[j], split = ";", fixed = T)[[1]][1], split = ".", fixed = T)[[1]][1])
    }else if(CS.list[[i]]$V1[j] == "ChrSy" | CS.list[[i]]$V1[j] == "ChrUn"){
      CS.list[[i]]$V13[j] <- gsub(pattern = "ID=", replacement = "", strsplit(CS.list[[i]]$V13[j], split = ";", fixed = T)[[1]][1])
    }
      }
  CS.genes[[i]] <- CS.list[[i]]$V13
  CS.genes[[i]] <- data.frame(Genes = CS.genes[[i]][!duplicated(CS.genes[[i]])], Type = NA)
  CS.genes[[i]] <- CS.genes[[i]][-which((CS.genes[[i]]$Genes %in% Gene_info$locus) == FALSE),]
  Iter.rows <- which((CS.genes[[i]]$Genes %in% Gene_info$locus) == TRUE)
  if(length(Iter.rows) != nrow(CS.genes[[i]])){ 
    stop("Something is going bad bro...")
  }
  for(j in Iter.rows){
    if(Gene_info$is_TE[which(Gene_info$locus == CS.genes[[i]]$Genes[j])] == "Y"){
      CS.genes[[i]]$Type[j] <- "transposable_element_gene"  
    }else if(Gene_info$is_TE[which(Gene_info$locus == CS.genes[[i]]$Genes[j])] != "Y"){
      CS.genes[[i]]$Type[j] <- "protein_coding"
    }
  }
}

# now it is orthology time: Phytozome ones

CS.genes <- lapply(CS.genes, function(x){
  x$Orthology <- NA
  x$Species <- NA
  x
})

# Simplier is the best - Orthologs Os, Orthologs Zm, Orthologs both, non-othologs
Os.orthology.Zm <- read.table("H:/Taiotactical/Taiotactical_ChIPseq/Utils/Annotation/0.Os/orthology/inParanoid_Zmays_493_RefGen_V4_Osativa_323_v7.0", header = F, sep = "\t")
Os.orthology.At <- read.table("H:/Taiotactical/Taiotactical_ChIPseq/Utils/Annotation/0.At/Orthology/inParanoid_Osativa_323_v7.0_Athaliana_447_Araport11", header = F, sep = "\t")
# Format this
Os.orthology.Zm$V2 <- gsub(pattern = " ", replacement = ".",Os.orthology.Zm$V2)
Os.orthology.Zm$V2 <- gsub(pattern = ".", replacement = "_",Os.orthology.Zm$V2, fixed = T)
Os.orthology.At$V2 <- gsub(pattern = " ", replacement = ".",Os.orthology.At$V2)
Os.orthology.At$V2 <- gsub(pattern = ".", replacement = "_",Os.orthology.At$V2, fixed = T)
Os.orthology.At$V1 <- gsub(pattern = " ", replacement = ".",Os.orthology.At$V1)
Os.orthology.At$V1 <- gsub(pattern = ".", replacement = "_",Os.orthology.At$V1, fixed = T)
for(i in 1:nrow(At.orthology.Zm)){
  V1 <- strsplit(At.orthology.Zm$V1[i], split = "_", fixed = T)[[1]][grep(x =  strsplit(At.orthology.Zm$V1[i], split = "_", fixed = T)[[1]], "Zmays:", fixed = T)]
  V2 <- strsplit(At.orthology.Zm$V2[i], split = "_", fixed = T)[[1]][grep(x =  strsplit(At.orthology.Zm$V2[i], split = "_", fixed = T)[[1]], "Athalianacolumbia:", fixed = T)]
  V1 <- strsplit(x = gsub(pattern = "Zmays:", replacement = "", x = V1), split = "_", fixed = T)[[1]]
  V2 <- gsub(pattern = "Athalianacolumbia:", replacement = "", x = V2)
  
}



# check Os and Zm genes and if both check one2one both one2one Os one2one Zm ...

for(i in 10:length(CS.genes)){
  for(j in 1:nrow(CS.genes[[i]])){
    # hit in At
    if(length(grep(CS.genes[[i]]$Genes[j], x = Os.orthology.At$V1)) > 0){
      Orthology.At <- "orthologous"
    }else if(length(grep(CS.genes[[i]]$Genes[j], x = Os.orthology.At$V1)) == 0){
      Orthology.At <- "non_orthologous"
    }
    # hit in Zm
    if(length(grep(CS.genes[[i]]$Genes[j], x = Os.orthology.Zm$V2)) > 0){
      Orthology.Zm <- "orthologous"
    }else if(length(grep(CS.genes[[i]]$Genes[j], x = Os.orthology.Zm$V2)) == 0){
      Orthology.Zm <- "non_orthologous"
    }    
    # coherency in both 
    if(Orthology.At == Orthology.Zm & Orthology.At == "orthologous"){
      Species <- "both"
      Orthology <- "orthologous"
    }else if(Orthology.At == Orthology.Zm & Orthology.At == "non_orthologous"){
      Species <- "both"
      Orthology <- "non_orthologous"
    }else if(Orthology.At != Orthology.Zm & Orthology.At == "orthologous"){
      Species <- "At"
      Orthology <- "orthologous"
    }else if(Orthology.At != Orthology.Zm & Orthology.Zm == "orthologous"){
      Species <- "Zm"
      Orthology <- "orthologous"
    }
    # fill data.frame for each state
    CS.genes[[i]]$Orthology[j] <- Orthology
    CS.genes[[i]]$Species[j] <- Species
  }
}


## Plot loop for each state accompanied in the title by total number of genes and prot coding genes: orthologous non-orthologous
## export tables just in case

CS.names <- c("CS1","CS10", "CS11", "CS12", "CS13", "CS14", "CS15", "CS16", "CS2", 
              "CS3", "CS4", "CS5", "CS6", "CS7", "CS8", "CS9")

for(i in 1:length(CS.genes)){
  # export tables
  #write.table(as.data.frame(table(CS.genes[[i]][,2:4])), paste0(inputDir, CS.names[i], "_Os_genestypeTable.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
  #write.table(CS.genes[[i]], paste0(inputDir, CS.names[i], "_Os_allgenes.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
  # plot alluvial
  CS.genes[[i]]$Species <- gsub("At", "cAt", CS.genes[[i]]$Species)
  allplot <- ggplot(as.data.frame(table(CS.genes[[i]][,c(2:4)])),
                    aes(y = Freq,
                        axis1 = Type, axis2 = Orthology)) +
    geom_alluvium(aes(fill = Species),
                  width = 1/8)  +
    geom_stratum(width = 1/8, fill = "white") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("Type", "Orthology"), expand = c(.05, .05)) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    ggtitle(paste0(CS.names[i], "Total Genes: ", nrow(CS.genes[[i]]), " Prot. coding genes: ", length(which(CS.genes[[i]]$Type == "protein_coding")), " Orthologous: ", length(which(CS.genes[[i]]$Orthology == "orthologous")))) + theme_classic()
  pdf(file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/figures_model1_15/", CS.names[i],"_alluvial_Os.pdf"), paper = "a4r", height = 28, width = 28, onefile = T)
  allplot
  dev.off()
}


################################################################## Zea mays
### Extract alluvial data for each Species: Zea mays mixed with TEs genes and other type of genes

Zm.df <- data.frame(Type = NA, Orthology = NA, Relation = NA, Freq = NA)
CS.list <- list()
CS.list.TE <- list()
CS.genes <- list()
CS.genes.TE <- list()
inputDir <- "H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/genes/Zm/"
CS.names <- c("CS1","CS10", "CS11", "CS12", "CS13", "CS14", "CS15", "CS16", "CS2", 
              "CS3", "CS4", "CS5", "CS6", "CS7", "CS8", "CS9")
# Create index data,frame in order to assign each to gene to different biotypes
Gene_info <- read.table(file = "H:/Taiotactical/Taiotactical_ChIPseq/Utils/Annotation/0.Zm/Annotation/Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gff3", header = T)
N.Gene_info <- Gene_info[grep("biotype", Gene_info$ID.chromosome.1),]
Gene_info <- N.Gene_info[grep("gene", N.Gene_info$chromosome),]
Genes <- list()
for(i in 1:nrow(Gene_info)){
  split <- strsplit(x = Gene_info$ID.chromosome.1[i], split = ";", )[[1]]
  Genes[[i]] <- c(gsub("ID=gene:","", split[grep(x = split, pattern = "ID=gene:")]),gsub("biotype=","", split[grep(x = split, pattern = "biotype")]))
}
Genes <- do.call(rbind, Genes)


for(i in 1:length(CS.names)){
  CS.list[[i]] <- read.table(paste0(inputDir, CS.names[i], "_genes.gff"), header = F, sep = "\t")
  CS.list.TE[[i]] <- read.table(paste0(inputDir, CS.names[i], "_TEs.gff"), header = F, sep = "\t")
  # format TE: by ID= and remove duplicates
  CS.list.TE[[i]]$V13 <- gsub("ID=", "",do.call(rbind,strsplit(CS.list.TE[[i]]$V13, split = ";", fixed = T))[,1])
  # format other type of genes
  CS.specific.genes <- list()
  for(j in 1:nrow(CS.list[[i]])){
      hit <- grep("Parent=", strsplit(CS.list[[i]]$V13[j], split = ";")[[1]])
      if(length(hit) > 0){
        CS.specific.genes[[j]] <- strsplit(strsplit(CS.list[[i]]$V13[j], split = ";")[[1]][hit], split = ":")[[1]][2]  
      }else if(length(hit) == 0){
        next
      }
  }
  CS.specific.genes <- do.call(rbind, CS.specific.genes)
  CS.specific.genes <- do.call(rbind,strsplit(CS.specific.genes, split = "_", fixed = T))[,1]
  CS.genes[[i]] <- data.frame(Genes = CS.specific.genes[!duplicated(CS.specific.genes)], Type = NA)
  CS.genes[[i]]$Genes <- gsub(pattern = "mir", replacement = "MIR", x = CS.genes[[i]]$Genes)
  Iter.rows <- which((CS.genes[[i]]$Genes %in% Genes[,1]) == TRUE)
  if(length(Iter.rows) != nrow(CS.genes[[i]])){ 
    stop("Something is going bad bro...")
  }
  for(j in Iter.rows){
      CS.genes[[i]]$Type[j] <- Genes[which(Genes[,1] == CS.genes[[i]]$Genes[j]),2] 
  }
  CS.genes.TE[[i]] <- data.frame(Genes = CS.list.TE[[i]]$V13[!duplicated(CS.list.TE[[i]]$V13)], Type = "transposable_element")
  CS.genes[[i]] <- rbind(CS.genes[[i]], CS.genes.TE[[i]])
}

# now it is orthology time: Phytozome ones

CS.genes <- lapply(CS.genes, function(x){
  x$Orthology <- NA
  x$Species <- NA
  x
})

# Simplier is the best - Orthologs Os, Orthologs Zm, Orthologs both, non-othologs
Zm.orthology.Os <- read.table("H:/Taiotactical/Taiotactical_ChIPseq/Utils/Annotation/0.Zm/orthology/inParanoid_Zmays_493_RefGen_V4_Osativa_323_v7.0", header = F, sep = "\t")
Zm.orthology.At <- read.table("H:/Taiotactical/Taiotactical_ChIPseq/Utils/Annotation/0.Zm/orthology/inParanoid_Zmays_493_RefGen_V4_Athaliana_447_Araport11", header = F, sep = "\t")
# Format this
Os.orthology.Zm$V2 <- gsub(pattern = " ", replacement = ".",Os.orthology.Zm$V2)
Os.orthology.Zm$V2 <- gsub(pattern = ".", replacement = "_",Os.orthology.Zm$V2, fixed = T)
Os.orthology.At$V2 <- gsub(pattern = " ", replacement = ".",Os.orthology.At$V2)
Os.orthology.At$V2 <- gsub(pattern = ".", replacement = "_",Os.orthology.At$V2, fixed = T)
Os.orthology.At$V1 <- gsub(pattern = " ", replacement = ".",Os.orthology.At$V1)
Os.orthology.At$V1 <- gsub(pattern = ".", replacement = "_",Os.orthology.At$V1, fixed = T)
for(i in 1:nrow(At.orthology.Zm)){
  V1 <- strsplit(At.orthology.Zm$V1[i], split = "_", fixed = T)[[1]][grep(x =  strsplit(At.orthology.Zm$V1[i], split = "_", fixed = T)[[1]], "Zmays:", fixed = T)]
  V2 <- strsplit(At.orthology.Zm$V2[i], split = "_", fixed = T)[[1]][grep(x =  strsplit(At.orthology.Zm$V2[i], split = "_", fixed = T)[[1]], "Athalianacolumbia:", fixed = T)]
  V1 <- strsplit(x = gsub(pattern = "Zmays:", replacement = "", x = V1), split = "_", fixed = T)[[1]]
  V2 <- gsub(pattern = "Athalianacolumbia:", replacement = "", x = V2)
  
}



# check Os and Zm genes and if both check one2one both one2one Os one2one Zm ...

for(i in 1:length(CS.genes)){
  for(j in 1:nrow(CS.genes[[i]])){
    # hit in At
    if(length(grep(CS.genes[[i]]$Genes[j], x = Zm.orthology.At$V1)) > 0){
      Orthology.At <- "orthologous"
    }else if(length(grep(CS.genes[[i]]$Genes[j], x = Zm.orthology.At$V1)) == 0){
      Orthology.At <- "non_orthologous"
    }
    # hit in Zm
    if(length(grep(CS.genes[[i]]$Genes[j], x = Zm.orthology.Os$V1)) > 0){
      Orthology.Os <- "orthologous"
    }else if(length(grep(CS.genes[[i]]$Genes[j], x = Zm.orthology.Os$V1)) == 0){
      Orthology.Os <- "non_orthologous"
    }    
    # coherency in both 
    if(Orthology.At == Orthology.Os & Orthology.At == "orthologous"){
      Species <- "both"
      Orthology <- "orthologous"
    }else if(Orthology.At == Orthology.Os & Orthology.At == "non_orthologous"){
      Species <- "both"
      Orthology <- "non_orthologous"
    }else if(Orthology.At != Orthology.Os & Orthology.At == "orthologous"){
      Species <- "At"
      Orthology <- "orthologous"
    }else if(Orthology.At != Orthology.Os & Orthology.Os == "orthologous"){
      Species <- "Os"
      Orthology <- "orthologous"
    }
    # fill data.frame for each state
    CS.genes[[i]]$Orthology[j] <- Orthology
    CS.genes[[i]]$Species[j] <- Species
  }
}


## Plot loop for each state accompanied in the title by total number of genes and prot coding genes: orthologous non-orthologous
## export tables just in case

CS.names <- c("CS1","CS10", "CS11", "CS12", "CS13", "CS14", "CS15", "CS16", "CS2", 
              "CS3", "CS4", "CS5", "CS6", "CS7", "CS8", "CS9")

for(i in 1:length(CS.genes)){
  # export tables
  write.table(as.data.frame(table(CS.genes[[i]][,2:4])), paste0(inputDir, CS.names[i], "_Zm_genestypeTable.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
  write.table(CS.genes[[i]], paste0(inputDir, CS.names[i], "_Zm_allgenes.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
  # plot alluvial
  CS.genes[[i]]$Species <- gsub("At", "cAt", CS.genes[[i]]$Species)
  allplot <- ggplot(as.data.frame(table(CS.genes[[i]][,c(2:4)])),
                    aes(y = Freq,
                        axis1 = Type, axis2 = Orthology)) +
    geom_alluvium(aes(fill = Species),
                  width = 1/8)  +
    geom_stratum(width = 1/8, fill = "white") +
    geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(limits = c("Type", "Orthology"), expand = c(.05, .05)) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    ggtitle(paste0(CS.names[i], "Total Genes: ", nrow(CS.genes[[i]]), " Prot. coding genes: ", length(which(CS.genes[[i]]$Type == "protein_coding")), " Orthologous: ", length(which(CS.genes[[i]]$Orthology == "orthologous")))) + theme_classic()
  pdf(file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/figures_model1_15/", CS.names[i],"_alluvial_Zm.pdf"), paper = "a4r", height = 28, width = 28, onefile = T)
  print(allplot)
  dev.off()
}


########################
# Ontology Enrichments #
########################

### GO biological Process: Only for protein coding genes

At <- read.table("H:/Taiotactical/Taiotactical_ChIPseq/Utils/Annotation/0.At/Annotation/Athaliana_447_Araport11.annotation_info.txt", sep = "\t", header = F)
Os <- read.delim("clipboard", header = T)
Zm <- read.delim("clipboard", header = T)

## 0. Generate background / universe and classify terms for biological process / term2gene

s <- strsplit(At$V10, split = ",")
At <- data.frame(Genes = rep(At$V2, sapply(s, length)), GOs = unlist(s))
At <- At[,c(2,1)]
At.GOmap <- buildGOmap(gomap = At)
At.Ont <- go2ont(At.GOmap$GO)
At.background.genes <- At.GOmap[At.GOmap$GO %in% At.Ont$go_id[which(At.Ont$Ontology == "CC")],] 

s <- strsplit(Os$GO, split = ",")
Os <- data.frame(Genes = rep(Os$locusName, sapply(s, length)), GOs = unlist(s))
Os <- Os[,c(2,1)]
Os.GOmap <- buildGOmap(gomap = Os)
Os.Ont <- go2ont(Os.GOmap$GO)
Os.background.genes <- Os.GOmap[Os.GOmap$GO %in% Os.Ont$go_id[which(Os.Ont$Ontology == "CC")],] 

s <- strsplit(Zm$GO, split = ",")
Zm <- data.frame(Genes = rep(Zm$locusName, sapply(s, length)), GOs = unlist(s))
Zm <- Zm[,c(2,1)] 
Zm.GOmap <- buildGOmap(gomap = Zm)
Zm.Ont <- go2ont(Zm.GOmap$GO)
Zm.background.genes <- Zm.GOmap[Zm.GOmap$GO %in% Zm.Ont$go_id[which(Zm.Ont$Ontology == "CC")],] 

## 2. Term2name

At.Term2name <- go2term(goid = At.background.genes$GO)
At.terms <- list()
Remove <- list()
for(i in 1:nrow(At.background.genes)){
  hit <- which(At.Term2name$go_id == At.background.genes$GO[i])
  if(length(hit) == 0){
    Remove[[i]] <- i
    next
  }
  At.terms[[i]] <- At.Term2name$Term[hit]
}

At.background.genes <- At.background.genes[-do.call(rbind, Remove)[,1],]
At.background.genes$Terms <- do.call(rbind, At.terms)[,1]
write.table(At.background.genes, file = "H:/Taiotactical/Taiotactical_ChIPseq/Utils/Annotation/0.At/At_GO_CC_Universe.txt", quote = F, sep = "\t", col.names = T, row.names = F)
At.background.genes <- read.table("H:/Taiotactical/Taiotactical_ChIPseq/Utils/Annotation/0.At/At_GO_Universe.txt", sep = "\t", header = T)


Os.Term2name <- go2term(goid = Os.background.genes$GO)
Os.terms <- list()
Remove <- list()
for(i in 1:nrow(Os.background.genes)){
  hit <- which(Os.Term2name$go_id == Os.background.genes$GO[i])
  if(length(hit) == 0){
    Remove[[i]] <- i
    next
  }
  Os.terms[[i]] <- Os.Term2name$Term[hit]
}

Os.background.genes <- Os.background.genes[-do.call(rbind, Remove)[,1],]
Os.background.genes$Terms <- do.call(rbind, Os.terms)[,1]
write.table(Os.background.genes, file = "H:/Taiotactical/Taiotactical_ChIPseq/Utils/Annotation/0.Os/Os_GO_CC_Universe.txt", quote = F, sep = "\t", col.names = T, row.names = F)
Os.background.genes <- read.table("H:/Taiotactical/Taiotactical_ChIPseq/Utils/Annotation/0.Os/Os_GO_Universe.txt", sep = "\t", header = T)

Zm.Term2name <- go2term(goid = Zm.background.genes$GO)
Zm.terms <- list()
Remove <- list()
for(i in 1:nrow(Zm.background.genes)){
  hit <- which(Zm.Term2name$go_id == Zm.background.genes$GO[i])
  if(length(hit) == 0){
    Remove[[i]] <- i
    next
  }
  Zm.terms[[i]] <- Zm.Term2name$Term[hit]
}

Zm.background.genes <- Zm.background.genes[-do.call(rbind, Remove)[,1],]
Zm.background.genes$Terms <- do.call(rbind, Zm.terms)[,1]
write.table(Zm.background.genes, file = "H:/Taiotactical/Taiotactical_ChIPseq/Utils/Annotation/0.Zm/Zm_GO_CC_Universe.txt", quote = F, sep = "\t", col.names = T, row.names = F)
Zm.background.genes <- read.table("H:/Taiotactical/Taiotactical_ChIPseq/Utils/Annotation/0.Zm/Zm_GO_Universe.txt", sep = "\t", header = T)


## 4. For each set generate enrichments, export results, plotting with and without filter by GO level

Species <- c("At", "Os", "Zm")
Background <- list(At = At.background.genes, Os = Os.background.genes, Zm = Zm.background.genes)
CS.names <- c("CS1","CS10", "CS11", "CS12", "CS13", "CS14", "CS15", "CS16", "CS2", 
              "CS3", "CS4", "CS5", "CS6", "CS7", "CS8", "CS9")

for(j in 1:length(CS.names)){
    CS.genes.At <- read.table(paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/genes/", names(Background[1]),"/", CS.names[j], "_", names(Background[1]), "_allgenes.txt"), header = T)
    CS.genes.Os <- read.table(paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/genes/", names(Background[2]),"/", CS.names[j], "_", names(Background[2]), "_allgenes.txt"), header = T)
    CS.genes.Zm <- read.table(paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/genes/", names(Background[3]),"/", CS.names[j], "_", names(Background[3]), "_allgenes.txt"), header = T)
    myGenes.At <- CS.genes.At$Genes[CS.genes.At$Genes %in% Background[[1]]$Gene]
    myGenes.Os <- CS.genes.Os$Genes[CS.genes.Os$Genes %in% Background[[2]]$Gene]
    myGenes.Zm <- CS.genes.Zm$Genes[CS.genes.Zm$Genes %in% Background[[3]]$Gene]
    myResults.At <- enricher(gene = myGenes.At, universe = levels(as.factor(Background[[1]]$Gene)), TERM2GENE = Background[[1]][,c(1,2)], TERM2NAME = Background[[1]][,c(1,3)])
    myResults.Os <- enricher(gene = myGenes.Os, universe = levels(as.factor(Background[[2]]$Gene)), TERM2GENE = Background[[2]][,c(1,2)], TERM2NAME = Background[[2]][,c(1,3)])
    myResults.Zm <- enricher(gene = myGenes.Zm, universe = levels(as.factor(Background[[3]]$Gene)), TERM2GENE = Background[[3]][,c(1,2)], TERM2NAME = Background[[3]][,c(1,3)])
    myResults.At@result$Species <- rep("At", nrow(myResults.At@result))
    myResults.Os@result$Species <- rep("Os", nrow(myResults.Os@result))
    myResults.Zm@result$Species <- rep("Zm", nrow(myResults.Zm@result))
    myResults.At@ontology <- "CC"
    myResults.Os@ontology <- "CC"
    myResults.Zm@ontology <- "CC"
    write.table(myResults.At@result, file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/enrichments/", names(Background[1]), "_",CS.names[j], "_GO_CC_enrichments.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
    write.table(myResults.Os@result, file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/enrichments/", names(Background[2]), "_",CS.names[j], "_GO_CC_enrichments.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
    write.table(myResults.Zm@result, file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/enrichments/", names(Background[3]), "_",CS.names[j], "_GO_CC_enrichments.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
    dotplot.At <- dotplot(myResults.At, showCategory = 250, size = "GeneRatio", x = "Species")
    dotplot.At.filter <- dotplot(gofilter(myResults.At, level = 5), showCategory = 200, size = "GeneRatio", x = "Species")
    dotplot.Os <- dotplot(myResults.Os, showCategory = 250, size = "GeneRatio", x = "Species")
    dotplot.Os.filter <- dotplot(gofilter(myResults.Os, level = 5), showCategory = 200, size = "GeneRatio", x = "Species")
    dotplot.Zm <- dotplot(myResults.Zm, showCategory = 250, size = "GeneRatio", x = "Species")
    dotplot.Zm.filter <- dotplot(gofilter(myResults.Zm, level = 5), showCategory = 200, size = "GeneRatio", x = "Species")
    pdf(file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/figures_model1_15/", CS.names[j],"_GO_CC_enrichments.pdf"), paper = "a4r", height = 28, width = 28, onefile = T)
    print(dotplot.At)
    print(dotplot.Os)
    print(dotplot.Zm)
    print(dotplot(gofilter(myResults.At, level = 2), showCategory = 200, size = "GeneRatio", x = "Species"))
    print(dotplot(gofilter(myResults.Os, level = 2), showCategory = 200, size = "GeneRatio", x = "Species"))
    print(dotplot(gofilter(myResults.Zm, level = 2), showCategory = 200, size = "GeneRatio", x = "Species"))
    print(dotplot(gofilter(myResults.At, level = 3), showCategory = 200, size = "GeneRatio", x = "Species"))
    print(dotplot(gofilter(myResults.Os, level = 3), showCategory = 200, size = "GeneRatio", x = "Species"))
    print(dotplot(gofilter(myResults.Zm, level = 3), showCategory = 200, size = "GeneRatio", x = "Species"))
    print(dotplot(gofilter(myResults.At, level = 4), showCategory = 200, size = "GeneRatio", x = "Species"))
    print(dotplot(gofilter(myResults.Os, level = 4), showCategory = 200, size = "GeneRatio", x = "Species"))
    print(dotplot(gofilter(myResults.Zm, level = 4), showCategory = 200, size = "GeneRatio", x = "Species"))
    print(dotplot.At.filter)
    print(dotplot.Os.filter)
    print(dotplot.Zm.filter)
    print(dotplot(gofilter(myResults.At, level = 6), showCategory = 200, size = "GeneRatio", x = "Species"))
    print(dotplot(gofilter(myResults.Os, level = 6), showCategory = 200, size = "GeneRatio", x = "Species"))
    print(dotplot(gofilter(myResults.Zm, level = 6), showCategory = 200, size = "GeneRatio", x = "Species"))
    print(dotplot(gofilter(myResults.At, level = 7), showCategory = 200, size = "GeneRatio", x = "Species"))
    print(dotplot(gofilter(myResults.Os, level = 7), showCategory = 200, size = "GeneRatio", x = "Species"))
    print(dotplot(gofilter(myResults.Zm, level = 7), showCategory = 200, size = "GeneRatio", x = "Species"))
    print(dotplot(gofilter(myResults.At, level = 8), showCategory = 200, size = "GeneRatio", x = "Species"))
    print(dotplot(gofilter(myResults.Os, level = 8), showCategory = 200, size = "GeneRatio", x = "Species"))
    print(dotplot(gofilter(myResults.Zm, level = 8), showCategory = 200, size = "GeneRatio", x = "Species"))
    dev.off()
    }

# Repeat for MF and CC above code

## 5.REVIGO Treemaps could be a better idea since some CS do not have enrichments at all (with white borders)

# CS.genes %in% GO.map and write table to REVIGO web ... 

for(j in 1:length(CS.names)){
  CS.genes.At <- read.table(paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/genes/At", "/", CS.names[j], "_", names(Background[1]), "_allgenes.txt"), header = T)
  CS.genes.Os <- read.table(paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/genes/Os", "/", CS.names[j], "_", names(Background[2]), "_allgenes.txt"), header = T)
  CS.genes.Zm <- read.table(paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/genes/Zm", "/", CS.names[j], "_", names(Background[3]), "_allgenes.txt"), header = T)
  GO.At <- At.GOmap[At.GOmap$Gene %in% CS.genes.At$Genes,1]
  GO.Os <- Os.GOmap[Os.GOmap$Gene %in% CS.genes.Os$Genes,1]
  GO.Zm <- Zm.GOmap[Zm.GOmap$Gene %in% CS.genes.Zm$Genes,1]
  write.table(GO.At, file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/enrichments/At", "_", CS.names[j], "_REVIGO.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
  write.table(GO.Os, file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/enrichments/Os", "_", CS.names[j], "_REVIGO.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
  write.table(GO.Zm, file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/enrichments/Zm", "_", CS.names[j], "_REVIGO.txt"), col.names = F, row.names = F, quote = F, sep = "\t")
  
}

# CS.genes filter enricher function results to print REVIGO associated to pvalues

CS.names <- c("CS1","CS10", "CS11", "CS12", "CS13", "CS14", "CS15", "CS16", "CS2", 
              "CS3", "CS4", "CS5", "CS6", "CS7", "CS8", "CS9")
Species <- c("At", "Os", "Zm")
for(j in 1:length(CS.names)){
  for(i in 1:length(Species)){
    Enricher.res.BP <- read.delim(paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/enrichments/",Species[i], "_", CS.names[j], "_GO_BP_enrichments.txt"), header = T, sep = "\t")
    Enricher.res.CC <- read.delim(paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/enrichments/",Species[i], "_", CS.names[j], "_GO_CC_enrichments.txt"), header = T, sep = "\t")
    Enricher.res.MF <- read.delim(paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/enrichments/",Species[i], "_", CS.names[j], "_GO_MF_enrichments.txt"), header = T, sep = "\t")
    Enricher.res <- rbind(Enricher.res.BP, Enricher.res.CC, Enricher.res.MF)
    Enricher.res <- Enricher.res[,c(1,6)]
    write.table(Enricher.res, sep = "\t", col.names = F, row.names = F, file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/enrichments/",Species[i], "_", CS.names[j], "_GO_Allfitered_enrichments.txt"), quote = F)
  }
}

### KO KEGG orthology

At <- read.delim("clipboard", header = T)
s <- strsplit(At$KO, split = ",")
At <- data.frame(Genes = rep(At$locusName, sapply(s, length)), KO = unlist(s))

Os <- read.delim("clipboard", header = T)
s <- strsplit(Os$KO, split = ",")
Os <- data.frame(Genes = rep(Os$locusName, sapply(s, length)), KO = unlist(s))

Zm <- read.delim("clipboard", header = T)
s <- strsplit(Zm$KO, split = ",")
Zm <- data.frame(Genes = rep(Zm$locusName, sapply(s, length)), KO = unlist(s))

Species <- c("At", "Os", "Zm")
Background <- list(At = At, Os = Os, Zm = Zm)
CS.names <- c("CS1","CS10", "CS11", "CS12", "CS13", "CS14", "CS15", "CS16", "CS2", 
              "CS3", "CS4", "CS5", "CS6", "CS7", "CS8", "CS9")

for(j in 1:length(CS.names)){
  CS.genes.At <- read.table(paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/genes/", names(Background[1]),"/", CS.names[j], "_", names(Background[1]), "_allgenes.txt"), header = T)
  CS.genes.Os <- read.table(paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/genes/", names(Background[2]),"/", CS.names[j], "_", names(Background[2]), "_allgenes.txt"), header = T)
  CS.genes.Zm <- read.table(paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/genes/", names(Background[3]),"/", CS.names[j], "_", names(Background[3]), "_allgenes.txt"), header = T)
  myKOs.At <- Background[[1]]$KO[Background[[1]]$Genes %in% CS.genes.At$Genes]
  myKOs.Os <- Background[[2]]$KO[Background[[2]]$Genes %in% CS.genes.Os$Genes]
  myKOs.Zm <- Background[[3]]$KO[Background[[3]]$Genes %in% CS.genes.Zm$Genes]
  myResults.At <- enrichKEGG(gene = myKOs.At,  pvalueCutoff = 0.05, keyType = "kegg", organism = "ko")
  myResults.Os <- enrichKEGG(gene = myKOs.Os,  pvalueCutoff = 0.05, keyType = "kegg", organism = "ko")
  myResults.Zm <- enrichKEGG(gene = myKOs.Zm,  pvalueCutoff = 0.05, keyType = "kegg", organism = "ko")
  myResults.M.At <- enrichMKEGG(gene = myKOs.At,  pvalueCutoff = 0.05, keyType = "kegg", organism = "ko")
  myResults.M.Os <- enrichMKEGG(gene = myKOs.Os,  pvalueCutoff = 0.05, keyType = "kegg", organism = "ko")
  myResults.M.Zm <- enrichMKEGG(gene = myKOs.Zm,  pvalueCutoff = 0.05, keyType = "kegg", organism = "ko")
  myResults.At.Universe <- enrichKEGG(gene = myKOs.At,  pvalueCutoff = 0.05, keyType = "kegg", organism = "ko", universe = levels(as.factor(Background[[1]]$KO)))
  myResults.Os.Universe <- enrichKEGG(gene = myKOs.Os,  pvalueCutoff = 0.05, keyType = "kegg", organism = "ko", universe = levels(as.factor(Background[[2]]$KO)))
  myResults.Zm.Universe <- enrichKEGG(gene = myKOs.Zm,  pvalueCutoff = 0.05, keyType = "kegg", organism = "ko", universe = levels(as.factor(Background[[3]]$KO)))
  myResults.M.At.Universe <- enrichMKEGG(gene = myKOs.At,  pvalueCutoff = 0.05, keyType = "kegg", organism = "ko", universe = levels(as.factor(Background[[1]]$KO)))
  myResults.M.Os.Universe <- enrichMKEGG(gene = myKOs.Os,  pvalueCutoff = 0.05, keyType = "kegg", organism = "ko", universe = levels(as.factor(Background[[2]]$KO)))
  myResults.M.Zm.Universe <- enrichMKEGG(gene = myKOs.Zm,  pvalueCutoff = 0.05, keyType = "kegg", organism = "ko", universe = levels(as.factor(Background[[3]]$KO)))
  write.table(myResults.At@result, file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/enrichments/", names(Background[1]), "_",CS.names[j], "_KO_enrichments.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
  write.table(myResults.Os@result, file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/enrichments/", names(Background[2]), "_",CS.names[j], "_KO_enrichments.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
  write.table(myResults.Zm@result, file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/enrichments/", names(Background[3]), "_",CS.names[j], "_KO_enrichments.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
  write.table(myResults.M.At@result, file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/enrichments/", names(Background[1]), "_",CS.names[j], "_KO_M_enrichments.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
  write.table(myResults.M.Os@result, file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/enrichments/", names(Background[2]), "_",CS.names[j], "_KO_M_enrichments.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
  write.table(myResults.M.Zm@result, file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/enrichments/", names(Background[3]), "_",CS.names[j], "_KO_M_enrichments.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
  write.table(myResults.At.Universe@result, file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/enrichments/", names(Background[1]), "_",CS.names[j], "_KO_Universe_enrichments.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
  write.table(myResults.Os.Universe@result, file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/enrichments/", names(Background[2]), "_",CS.names[j], "_KO_Universe_enrichments.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
  write.table(myResults.Zm.Universe@result, file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/enrichments/", names(Background[3]), "_",CS.names[j], "_KO_Universe_enrichments.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
  write.table(myResults.M.At.Universe@result, file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/enrichments/", names(Background[1]), "_",CS.names[j], "_KO_M_Universe_enrichments.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
  write.table(myResults.M.Os.Universe@result, file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/enrichments/", names(Background[2]), "_",CS.names[j], "_KO_M_Universe_enrichments.txt"), col.names = T, row.names = F, quote = F, sep = "\t")
  write.table(myResults.M.Zm.Universe@result, file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/model1_15/Remapped/splitted/enrichments/", names(Background[3]), "_",CS.names[j], "_KO_M_Universe_enrichments.txt"), col.names = T, row.names = F, quote = F, sep = "\t")  
  dotplot.At <- dotplot(myResults.At, showCategory = 250, size = "GeneRatio", x = "GeneRatio")
  dotplot.Os <- dotplot(myResults.Os, showCategory = 250, size = "GeneRatio", x = "GeneRatio")
  dotplot.Zm <- dotplot(myResults.Zm, showCategory = 250, size = "GeneRatio", x = "GeneRatio")
  dotplot.At.M <- dotplot(myResults.M.At, showCategory = 250, size = "GeneRatio", x = "GeneRatio")
  dotplot.Os.M <- dotplot(myResults.M.Os, showCategory = 250, size = "GeneRatio", x = "GeneRatio")
  dotplot.Zm.M <- dotplot(myResults.M.Zm, showCategory = 250, size = "GeneRatio", x = "GeneRatio")
  dotplot.At.Universe <- dotplot(myResults.At.Universe, showCategory = 250, size = "GeneRatio", x = "GeneRatio")
  dotplot.Os.Universe <- dotplot(myResults.Os.Universe, showCategory = 250, size = "GeneRatio", x = "GeneRatio")
  dotplot.Zm.Universe <- dotplot(myResults.Zm.Universe, showCategory = 250, size = "GeneRatio", x = "GeneRatio")
  dotplot.At.Universe.M <- dotplot(myResults.M.At.Universe, showCategory = 250, size = "GeneRatio", x = "GeneRatio")
  dotplot.Os.Universe.M <- dotplot(myResults.M.Os.Universe, showCategory = 250, size = "GeneRatio", x = "GeneRatio")
  dotplot.Zm.Universe.M <- dotplot(myResults.M.Zm.Universe, showCategory = 250, size = "GeneRatio", x = "GeneRatio")
  pdf(file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/1.hiHMM_states/outputs/figures_model1_15/", CS.names[j],"_KO_enrichments.pdf"), paper = "a4r", height = 28, width = 28, onefile = T)

  print(dotplot.At)
  print(dotplot.Os)
  print(dotplot.Zm)
  print(dotplot.At.M)
  print(dotplot.Os.M)
  print(dotplot.Zm.M)
  print(dotplot.At.Universe)
  print(dotplot.Os.Universe)
  print(dotplot.Zm.Universe)
  print(dotplot.At.Universe.M)
  print(dotplot.Os.Universe.M)
  print(dotplot.Zm.Universe.M)

  dev.off()
}

##### BONUS: Interspecies Corrplot

library(pheatmap)
library(RColorBrewer)

At.pearson <- read.delim("clipboard", header = T, row.names = 1)
At.spearman <- read.delim("clipboard", header = T, row.names = 1)
Os.pearson <- read.delim("clipboard", header = T, row.names = 1)
Os.spearman <- read.delim("clipboard", header = T, row.names = 1)
Zm.pearson <- read.delim("clipboard", header = T, row.names = 1)
Zm.spearman <- read.delim("clipboard", header = T, row.names = 1)

pearson <- list(At = At.pearson, Os = Os.pearson, Zm = Zm.pearson)
spearman <- list(At = At.spearman, Os = Os.spearman, Zm = Zm.spearman)

pearson <- lapply(pearson, function(x){
  colnames(x) <- rownames(x)
  x
})
spearman <- lapply(spearman, function(x){
  colnames(x) <- rownames(x)
  x
})

variance <- as.data.frame(matrix(NA, 11,11))
colnames(variance) <- colnames(pearson$At)
row.names(variance) <-  colnames(pearson$At)

averaged <- as.data.frame(matrix(NA, 11,11))
colnames(averaged) <- colnames(pearson$At)
row.names(averaged) <-  colnames(pearson$At)

for(i in rownames(pearson$At)){
  for(j in colnames(pearson$At)){
    Corr.At <- spearman$At[ which(row.names(spearman$At) == i),which(colnames(spearman$At) == j)]
    Corr.Os <- spearman$Os[ which(row.names(spearman$Os) == i),which(colnames(spearman$Os) == j)]
    Corr.Zm <- spearman$Zm[ which(row.names(spearman$Zm) == i),which(colnames(spearman$Zm) == j)]
    
    averaged[ which(row.names(averaged) == i),which(colnames(averaged) == j)] <- mean(c(Corr.At, Corr.Os, Corr.Zm))
    variance[ which(row.names(variance) == i),which(colnames(variance) == j)] <- var(c(Corr.At, Corr.Os, Corr.Zm))
    
  }
}

col <- colorRampPalette(brewer.pal(11, "RdBu"))(250)
pdf(file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/0.Taiotactical/Correlation/Correlation/Output/AveragedInterSpeciesCorrelation_Spearman.pdf"), paper = "a4r", height = 21, width = 28, onefile = T)
pheatmap(averaged, scale = "none", color = col,  cluster_rows = T, cluster_cols = T,  clustering_method = "ward.D2", cellwidth = 35, cellheight = 35,main = "Averaged Inter-species correlation spearman", border_color = "white")
dev.off()
write.table(averaged, file = "H:/Taiotactical/Taiotactical_ChIPseq/0.Taiotactical/Correlation/Correlation/Output/AveragedInterSpeciesCorrelation_Spearman.txt", sep = "\t", row.names = T, col.names = T, quote = F)

col <- colorRampPalette(brewer.pal(7, "Greys"))(250)
pdf(file = paste0("H:/Taiotactical/Taiotactical_ChIPseq/0.Taiotactical/Correlation/Correlation/Output/VarianceInterSpeciesCorrelation_Spearman.pdf"), paper = "a4r", height = 21, width = 28, onefile = T)
pheatmap(variance, scale = "none", color = col,  cluster_rows = T, cluster_cols = T,  clustering_method = "ward.D2", cellwidth = 35, cellheight = 35,main = "Variance Inter-species correlation spearman", border_color = "white")
dev.off()
write.table(variance, file = "H:/Taiotactical/Taiotactical_ChIPseq/0.Taiotactical/Correlation/Correlation/Output/VarianceInterSpeciesCorrelation_Spearman.txt", sep = "\t", row.names = T, col.names = T, quote = F)



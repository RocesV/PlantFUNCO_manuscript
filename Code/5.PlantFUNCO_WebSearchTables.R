###################################
# TAIOTACTICAL-Web Search Tables  #
###################################

#### 0. Load libs and pkgs ####

library(GenomicFeatures)
library(htmlwidgets)
library(reactable)
library(sparkline)
library(dplyr)
library(RColorBrewer)
library(htmltools)
library(fontawesome)

#### 1. Gene tables ####

### 1. Extract gene (3 species) coordinates in order to bin each with bedtools makewindows

At <- makeTxDbFromGFF(paste0("H:/Taiotactical/Taiotactical_ChIPseq/Utils/Annotation/0.At/Annotation/0.At.gene_exons.gff3"))
Os <- makeTxDbFromGFF(paste0("H:/Taiotactical/Taiotactical_ChIPseq/Utils/Annotation/0.Os/Annotation/0.Os.gene_exons.gff3"))
Zm <- makeTxDbFromGFF(paste0("H:/Taiotactical/Taiotactical_ChIPseq/Utils/Annotation/0.Zm/Annotation/0.Zm.gene_exons.gff3"))

Genes.At <- genes(At)
Genes.Os <- genes(Os)
Genes.Zm <- genes(Zm)

Genes.At.df <- data.frame(Chr = as.character(seqnames(Genes.At)), Start = Genes.At@ranges@start, End = Genes.At@ranges@start + Genes.At@ranges@width - 1, GeneID = Genes.At$gene_id)
Genes.Os.df <- data.frame(Chr = as.character(seqnames(Genes.Os)), Start = Genes.Os@ranges@start, End = Genes.Os@ranges@start + Genes.Os@ranges@width - 1, GeneID = Genes.Os$gene_id)
Genes.Zm.df <- data.frame(Chr = as.character(seqnames(Genes.Zm)), Start = Genes.Zm@ranges@start, End = Genes.Zm@ranges@start + Genes.Zm@ranges@width - 1, GeneID = Genes.Zm$gene_id)

write.table(Genes.At.df, "H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Genes.At.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(Genes.Os.df, "H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Genes.Os.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(Genes.Zm.df, "H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Genes.Zm.txt", quote = F, sep = "\t", row.names = F, col.names = F)

### 2.Bedtools makewindows 10 bins per gene

# bash

### 3.Bedtools intersect LECIF score, PhyloP, CNEs, PhastCons, hiHMM CS

# bash

#### 4.Create reactable data ####

#### 4.1 Arabidopsis-Oryza ####

### Format ###

## LECIF

At_Os.LECIF <- read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Arabidopsis-Oryza/Genes.At-Os.sorted.binned.LECIFscore.bed", header = F, sep = "\t")
At_Os.LECIF$V8[which(At_Os.LECIF$V8 == ".")] <- gsub(".", 0, fixed = T, At_Os.LECIF$V8[which(At_Os.LECIF$V8 == ".")])
At_Os.LECIF$V8 <- as.numeric(At_Os.LECIF$V8)
At_Os.LECIF <- At_Os.LECIF[,c(1,2,3,4,8)]

## PhyloP

At_Os.PhyloP <- read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Arabidopsis-Oryza/Genes.At-Os.sorted.binned.PhyloPscore.sorted.bed", header = F, sep = "\t")
At_Os.PhyloP$V8 <- as.numeric(At_Os.PhyloP$V8)
At_Os.PhyloP <- At_Os.PhyloP[,c(1,2,3,4,8)]
At_Os.PhyloP$V8[which(is.na(At_Os.PhyloP$V8) == TRUE)] <- 0

## hiHMM Chromatin States (CS)

At_Os.CS <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Arabidopsis-Oryza/Genes.At-Os.sorted.binned.CS.bed", header = F, sep = "\t")
At_Os.CS <- At_Os.CS[,c(1,2,3,4,8)]

## PhastCons

At_Os.PhastCons <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Arabidopsis-Oryza/Genes.At-Os.sorted.binned.PhastCons.bed", header = F, sep = "\t")
At_Os.PhastCons <- At_Os.PhastCons[,c(1,2,3,4,7)]

## CNEs 

At_Os.CNEs <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Arabidopsis-Oryza/Genes.At-Os.sorted.CNEs.bed", header = F, sep = "\t")
At_Os.CNEs$CNEs_100kb <- "No"
At_Os.CNEs$CNEs_100kb[which(At_Os.CNEs$V14 <= 100000)] <- "Yes"
At_Os.CNEs$CNEsNearest_Distance <- At_Os.CNEs$V14
At_Os.CNEs$CNEsID <- NA
At_Os.CNEs$CNEsID <- At_Os.CNEs$V8
At_Os.CNEs$Gene <- At_Os.CNEs$V4
At_Os.CNEs <- At_Os.CNEs[,c(18,17,15,16)]

## SNPs-GWAS

At_Os.GWAS <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Arabidopsis-Oryza/Genes.At-Os.sorted.binned.SNPsGWAS.bed", header = F, sep = "\t")
At_Os.GWAS$GWAS.SNPs_Presence <- At_Os.GWAS$V8
At_Os.GWAS$GWAS.SNPs_Location <- paste0(At_Os.GWAS$V5,":",At_Os.GWAS$V6,":",At_Os.GWAS$V7)
At_Os.GWAS <- At_Os.GWAS[,c(1,2,3,4,9,10)]

### First group by interval*gene and mean n LECIF hits ###

## LECIF

At_Os.LECIF$PASTED <- paste0(At_Os.LECIF$V1, ":", At_Os.LECIF$V2, ":", At_Os.LECIF$V3) 
data <- At_Os.LECIF %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4),LECIFscore = list(V8))
data$LECIFscoreM <- unlist(lapply(data$LECIFscore, FUN = mean), use.names = F)

## PhyloP

At_Os.PhyloP$PASTED <- paste0(At_Os.PhyloP$V1, ":", At_Os.PhyloP$V2, ":", At_Os.PhyloP$V3)
data <- At_Os.PhyloP %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4),PhyloPscore = list(V8))
data$PhyloPscoreM <- unlist(lapply(data$PhyloPscore, FUN = mean), use.names = F)

## hiHMM Chromatin States (CS)

At_Os.CS$PASTED <- paste0(At_Os.CS$V1, ":", At_Os.CS$V2, ":", At_Os.CS$V3)
data <- At_Os.CS %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4),CS_Composition = paste(V8, collapse = ","))

## PhastCons 

At_Os.PhastCons$PASTED <- paste0(At_Os.PhastCons$V1, ":", At_Os.PhastCons$V2, ":", At_Os.PhastCons$V3)
At_Os.PhastCons$V7 <- gsub("phastCons_predicted", "Yes", fixed = T, x = At_Os.PhastCons$V7)
At_Os.PhastCons$V7 <- gsub(".", "No", fixed = T, x = At_Os.PhastCons$V7)
data <- At_Os.PhastCons %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4), PhastCons_Element = paste(V7, collapse = ","))

## CNEs
## Already grouped by gene intervals

## GWAS

At_Os.GWAS$PASTED <- paste0(At_Os.GWAS$V1, ":", At_Os.GWAS$V2, ":", At_Os.GWAS$V3)
At_Os.GWAS$GWAS.SNPs_Presence <- gsub("0", "No", fixed = T, x = At_Os.GWAS$GWAS.SNPs_Presence)
At_Os.GWAS$GWAS.SNPs_Presence <- gsub("1", "Yes", fixed = T, x = At_Os.GWAS$GWAS.SNPs_Presence)
data <- At_Os.GWAS %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4), GWAS.SNPs_Presence = paste(GWAS.SNPs_Presence, collapse = ","), GWAS.SNPs_Location = paste(GWAS.SNPs_Location, collapse = ","))

### Second group by gene only and mean five bins ###

## LECIF

data2 <- data %>%
  group_by(Gene) %>%
  summarise(LECIF_Trend = list(LECIFscoreM))

data2$LECIF_Mean <- unlist(lapply(data2$LECIF_Trend, mean), use.names = F)
data2$LECIF_Mean <- round(data2$LECIF_Mean, digits = 3)
data2$LECIF_Quartile <- NA
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0 & percent_rank(x = data2$LECIF_Mean) < 0.25)] <- "Q4"
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0.25 & percent_rank(x = data2$LECIF_Mean) < 0.5)] <- "Q3"
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0.5 & percent_rank(x = data2$LECIF_Mean) < 0.75)] <- "Q2"
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0.75 & percent_rank(x = data2$LECIF_Mean) <= 1)] <- "Q1"

data2$Species <- rep("Arabidopsis thaliana", nrow(data2))
data2$Comparison <- rep("Arabidopsis-Oryza", nrow(data2))
data2 <- data2[,c(5,6,1,3,4,2)]

## PhyloP

data2 <- data %>%
  group_by(Gene) %>%
  summarise(PhyloP_Trend = list(PhyloPscoreM))

data2$PhyloP_Mean <- unlist(lapply(data2$PhyloP_Trend, mean), use.names = F)
data2$PhyloP_Mean <- round(data2$PhyloP_Mean, digits = 3)
data2$PhyloP_Quartile <- NA
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0 & percent_rank(x = data2$PhyloP_Mean) < 0.25)] <- "Q4"
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0.25 & percent_rank(x = data2$PhyloP_Mean) < 0.5)] <- "Q3"
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0.5 & percent_rank(x = data2$PhyloP_Mean) < 0.75)] <- "Q2"
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0.75 & percent_rank(x = data2$PhyloP_Mean) <= 1)] <- "Q1"

data2$Species <- rep("Arabidopsis thaliana", nrow(data2))
data2$Comparison <- rep("Arabidopsis-Oryza", nrow(data2))
data2 <- data2[,c(5,6,1,3,4,2)]

## hiHMM Chromatin States (CS)

data2 <- data %>%
  group_by(Gene) %>%
  summarise(CS_Composition = paste(CS_Composition, collapse = ","))

data2$Species <- rep("Arabidopsis thaliana", nrow(data2))
data2$Comparison <- rep("Arabidopsis-Oryza", nrow(data2))
data2 <- data2[,c(3,4,1,2)]

## PhastCons

data2 <- data %>%
  group_by(Gene) %>%
  summarise(PhastCons_Element = paste(PhastCons_Element, collapse = ","))

data2$Species <- rep("Arabidopsis thaliana", nrow(data2))
data2$Comparison <- rep("Arabidopsis-Oryza", nrow(data2))
data2 <- data2[,c(3,4,1,2)]

## CNEs

data2 <- At_Os.CNEs %>%
  group_by(Gene) %>%
  summarise(CNEs_100kb = paste(CNEs_100kb, collapse = ","),
            CNEsNearest_Distance = unique(CNEsNearest_Distance),
            CNEsID = paste(CNEsID, collapse = ","))

data2$Species <- rep("Arabidopsis thaliana", nrow(data2))
data2$Comparison <- rep("Arabidopsis-Oryza", nrow(data2))
data2 <- data2[,c(5,6,1,2,3,4)]

## GWAS

data2 <- data %>%
  group_by(Gene) %>%
  summarise(GWAS.SNPs_Presence = paste(GWAS.SNPs_Presence, collapse = ","), GWAS.SNPs_Location = paste(GWAS.SNPs_Location, collapse = ","))

data2$GWAS.SNPs_Location <- gsub(".:-1:-1", "", data2$GWAS.SNPs_Location)
data2$GWAS.SNPs_Location <- gsub(",", "", data2$GWAS.SNPs_Location)
data2$Species <- rep("Arabidopsis thaliana", nrow(data2))
data2$Comparison <- rep("Arabidopsis-Oryza", nrow(data2))
data2 <- data2[,c(4,5,1,2,3)]

### Export & joint section ###

AtOsLECIF.table <- data2
AtOsLECIF.table.jic <- data2

AtOsPhyloP.table <- data2
AtOsPhyloP.table.jic <- data2

AtOsCS.table <- data2
AtOsCS.table.jic <- data2

AtOsPhastCons.table <- data2
AtOsPhastCons.table.jic <- data2

AtOsCNEs.table <- data2
AtOsCNEs.table.jic <- data2

AtOsGWAS.table <- data2
AtOsGWAS.table.jic <- data2

if(length(which((AtOsLECIF.table$Gene == AtOsPhyloP.table$Gene) == FALSE)) == 0){ print("Ready to join")}

AtOs.LECIF.PhyloP.table <- as_tibble(cbind(AtOsLECIF.table, AtOsPhyloP.table[, 4:6]))

if(length(which((AtOsCS.table$Gene == AtOs.LECIF.PhyloP.table$Gene) == FALSE)) == 0){ print("Ready to join")}

AtOs.LECIF.PhyloP.CS.table <- as_tibble(cbind(AtOs.LECIF.PhyloP.table, AtOsCS.table[,4]))

if(length(which((AtOsPhastCons.table$Gene == AtOs.LECIF.PhyloP.CS.table$Gene) == FALSE)) == 0){ print("Ready to join")}

AtOs.LECIF.PhyloP.CS.PhastCons.table <- as_tibble(cbind(AtOs.LECIF.PhyloP.CS.table, AtOsPhastCons.table[,4]))

AtOsCNEs.table <- AtOsCNEs.table[-5445,]
if(length(which((AtOsCNEs.table$Gene == AtOs.LECIF.PhyloP.CS.PhastCons.table$Gene) == FALSE)) == 0){ print("Ready to join")}

AtOs.LECIF.PhyloP.CS.PhastCons.CNEs.table <- as_tibble(cbind(AtOs.LECIF.PhyloP.CS.PhastCons.table, AtOsCNEs.table[,4:6]))

if(length(which((AtOsGWAS.table$Gene == AtOs.LECIF.PhyloP.CS.PhastCons.CNEs.table$Gene) == FALSE)) == 0){ print("Ready to join")}

AtOs.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table <- as_tibble(cbind(AtOs.LECIF.PhyloP.CS.PhastCons.CNEs.table, AtOsGWAS.table[,4:5]))

## Custom table

data2 <- AtOs.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table
sticky_style <- list(backgroundColor = "#e8e8e8")
html.object <- htmltools::browsable(
  tagList(
    tags$button(
      tagList(fontawesome::fa("download"), "Download Table"),
      onclick = "Reactable.downloadDataCSV('GenesTableFUNCO', 'GenesTableFUNCO.csv')"
    ),
    
    reactable(data2, searchable = T, filterable = T, pagination = TRUE, paginationType = "jump", defaultPageSize = 100, defaultColDef = colDef(vAlign = "center", headerVAlign = "center", align = "center", minWidth = 80, footer = function(values, name) { htmltools::div(name, style = list(fontWeight = 600))}), 
              bordered = T, striped = TRUE, highlight = TRUE, resizable = TRUE, wrap = FALSE, sortable = TRUE, elementId = "GenesTableFUNCO", defaultSorted = list(LECIF_Mean = "desc"), height = 800,
              columns = list(
                Species = colDef(sticky = "left",
                              style = sticky_style,
                              headerStyle = sticky_style),
                Comparison = colDef(sticky = "left",
                              style = sticky_style,
                              headerStyle = sticky_style),
                Gene = colDef(sticky = "left",
                              style = sticky_style,
                              headerStyle = sticky_style),
                LECIF_Quartile = colDef(
                  style = function(value){
                    if(value == "Q1"){
                      color <- "#309143"
                    } else if(value == "Q2"){
                      color <- "#8ace7e"
                    } else if(value == "Q3"){
                      color <- "#ff684c"
                    } else if(value == "Q4"){
                      color <- "#b60a1c"
                    }
                    list(color = color, fontWeight = "bold")
                  }
                ),
                LECIF_Trend = colDef(cell = function(value, index) {
                  sparkline(data2$LECIF_Trend[[index]])
                }),
                PhyloP_Quartile = colDef(
                  style = function(value){
                    if(value == "Q1"){
                      color <- "#309143"
                    } else if(value == "Q2"){
                      color <- "#8ace7e"
                    } else if(value == "Q3"){
                      color <- "#ff684c"
                    } else if(value == "Q4"){
                      color <- "#b60a1c"
                    }
                    list(color = color, fontWeight = "bold")
                  }
                ),
                PhyloP_Trend = colDef(cell = function(value, index) {
                  sparkline(data2$PhyloP_Trend[[index]])
                }),
                PhastCons_Element = colDef(cell = function(value) {
                  if (length(grep("Yes", value)) > 0) "\u2714\ufe0f Yes" else "\u274c No"
                }),
                CNEs_100kb = colDef(cell = function(value) {
                  if (length(grep("Yes", value)) > 0) "\u2714\ufe0f Yes" else "\u274c No"
                }),
                GWAS.SNPs_Presence = colDef(cell = function(value) {
                  if (length(grep("Yes", value)) > 0) "\u2714\ufe0f Yes" else "\u274c No"
                })
              ))
  )
)

save_html(html.object, "H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Arabidopsis-Oryza/SearchGenes_At-Os.html")

#### 4.2 Arabidopsis-Zea ####

### Format ###

## LECIF

At_Zm.LECIF <- read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Arabidopsis-Zea/Genes.At-Zm.sorted.binned.LECIFscore.bed", header = F, sep = "\t")
At_Zm.LECIF$V8[which(At_Zm.LECIF$V8 == ".")] <- gsub(".", 0, fixed = T, At_Zm.LECIF$V8[which(At_Zm.LECIF$V8 == ".")])
At_Zm.LECIF$V8 <- as.numeric(At_Zm.LECIF$V8)
At_Zm.LECIF <- At_Zm.LECIF[,c(1,2,3,4,8)]

## PhyloP

At_Zm.PhyloP <- read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Arabidopsis-Zea/Genes.At-Zm.sorted.binned.PhyloPscore.bed", header = F, sep = "\t")
At_Zm.PhyloP$V8 <- as.numeric(At_Zm.PhyloP$V8)
At_Zm.PhyloP <- At_Zm.PhyloP[,c(1,2,3,4,8)]
At_Zm.PhyloP$V8[which(is.na(At_Zm.PhyloP$V8) == TRUE)] <- 0

## hiHMM Chromatin States (CS)

At_Zm.CS <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Arabidopsis-Zea/Genes.At-Zm.sorted.binned.CS.bed", header = F, sep = "\t")
At_Zm.CS <- At_Zm.CS[,c(1,2,3,4,8)]

## PhastCons

At_Zm.PhastCons <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Arabidopsis-Zea/Genes.At-Zm.sorted.binned.PhastCons.bed", header = F, sep = "\t")
At_Zm.PhastCons <- At_Zm.PhastCons[,c(1,2,3,4,7)]

## CNEs 

At_Zm.CNEs <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Arabidopsis-Zea/Genes.At-Zm.sorted.CNEs.bed", header = F, sep = "\t")
At_Zm.CNEs$CNEs_100kb <- "No"
At_Zm.CNEs$CNEs_100kb[which(At_Zm.CNEs$V14 <= 100000)] <- "Yes"
At_Zm.CNEs$CNEsNearest_Distance <- At_Zm.CNEs$V14
At_Zm.CNEs$CNEsID <- NA
At_Zm.CNEs$CNEsID <- At_Zm.CNEs$V8
At_Zm.CNEs$Gene <- At_Zm.CNEs$V4
At_Zm.CNEs <- At_Zm.CNEs[,c(18,17,15,16)]

## GWAS

At_Zm.GWAS <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Arabidopsis-Zea/Genes.At-Zm.sorted.binned.SNPsGWAS.bed", header = F, sep = "\t")
At_Zm.GWAS$GWAS.SNPs_Presence <- At_Zm.GWAS$V8
At_Zm.GWAS$GWAS.SNPs_Location <- paste0(At_Zm.GWAS$V5,":",At_Zm.GWAS$V6,":",At_Zm.GWAS$V7)
At_Zm.GWAS <- At_Zm.GWAS[,c(1,2,3,4,9,10)]

### First group by interval*gene and mean n LECIF hits ###

## LECIF

At_Zm.LECIF$PASTED <- paste0(At_Zm.LECIF$V1, ":", At_Zm.LECIF$V2, ":", At_Zm.LECIF$V3) 
data <- At_Zm.LECIF %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4),LECIFscore = list(V8))
data$LECIFscoreM <- unlist(lapply(data$LECIFscore, FUN = mean), use.names = F)

## PhyloP

At_Zm.PhyloP$PASTED <- paste0(At_Zm.PhyloP$V1, ":", At_Zm.PhyloP$V2, ":", At_Zm.PhyloP$V3)
data <- At_Zm.PhyloP %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4),PhyloPscore = list(V8))
data$PhyloPscoreM <- unlist(lapply(data$PhyloPscore, FUN = mean), use.names = F)

## hiHMM Chromatin States (CS)

At_Zm.CS$PASTED <- paste0(At_Zm.CS$V1, ":", At_Zm.CS$V2, ":", At_Zm.CS$V3)
data <- At_Zm.CS %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4),CS_Composition = paste(V8, collapse = ","))

## PhastCons 

At_Zm.PhastCons$PASTED <- paste0(At_Zm.PhastCons$V1, ":", At_Zm.PhastCons$V2, ":", At_Zm.PhastCons$V3)
At_Zm.PhastCons$V7 <- gsub("phastCons_predicted", "Yes", fixed = T, x = At_Zm.PhastCons$V7)
At_Zm.PhastCons$V7 <- gsub(".", "No", fixed = T, x = At_Zm.PhastCons$V7)
data <- At_Zm.PhastCons %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4), PhastCons_Element = paste(V7, collapse = ","))

## CNEs
## Already grouped by gene intervals

## GWAS

At_Zm.GWAS$PASTED <- paste0(At_Zm.GWAS$V1, ":", At_Zm.GWAS$V2, ":", At_Zm.GWAS$V3)
At_Zm.GWAS$GWAS.SNPs_Presence <- gsub("0", "No", fixed = T, x = At_Zm.GWAS$GWAS.SNPs_Presence)
At_Zm.GWAS$GWAS.SNPs_Presence <- gsub("1", "Yes", fixed = T, x = At_Zm.GWAS$GWAS.SNPs_Presence)
data <- At_Zm.GWAS %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4), GWAS.SNPs_Presence = paste(GWAS.SNPs_Presence, collapse = ","), GWAS.SNPs_Location = paste(GWAS.SNPs_Location, collapse = ","))

### Second group by gene only and mean five bins ###

## LECIF

data2 <- data %>%
  group_by(Gene) %>%
  summarise(LECIF_Trend = list(LECIFscoreM))

data2$LECIF_Mean <- unlist(lapply(data2$LECIF_Trend, mean), use.names = F)
data2$LECIF_Mean <- round(data2$LECIF_Mean, digits = 3)
data2$LECIF_Quartile <- NA
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0 & percent_rank(x = data2$LECIF_Mean) < 0.25)] <- "Q4"
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0.25 & percent_rank(x = data2$LECIF_Mean) < 0.5)] <- "Q3"
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0.5 & percent_rank(x = data2$LECIF_Mean) < 0.75)] <- "Q2"
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0.75 & percent_rank(x = data2$LECIF_Mean) <= 1)] <- "Q1"

data2$Species <- rep("Arabidopsis thaliana", nrow(data2))
data2$Comparison <- rep("Arabidopsis-Zea", nrow(data2))
data2 <- data2[,c(5,6,1,3,4,2)]

## PhyloP

data2 <- data %>%
  group_by(Gene) %>%
  summarise(PhyloP_Trend = list(PhyloPscoreM))

data2$PhyloP_Mean <- unlist(lapply(data2$PhyloP_Trend, mean), use.names = F)
data2$PhyloP_Mean <- round(data2$PhyloP_Mean, digits = 3)
data2$PhyloP_Quartile <- NA
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0 & percent_rank(x = data2$PhyloP_Mean) < 0.25)] <- "Q4"
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0.25 & percent_rank(x = data2$PhyloP_Mean) < 0.5)] <- "Q3"
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0.5 & percent_rank(x = data2$PhyloP_Mean) < 0.75)] <- "Q2"
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0.75 & percent_rank(x = data2$PhyloP_Mean) <= 1)] <- "Q1"

data2$Species <- rep("Arabidopsis thaliana", nrow(data2))
data2$Comparison <- rep("Arabidopsis-Zea", nrow(data2))
data2 <- data2[,c(5,6,1,3,4,2)]

## hiHMM Chromatin States (CS)

data2 <- data %>%
  group_by(Gene) %>%
  summarise(CS_Composition = paste(CS_Composition, collapse = ","))

data2$Species <- rep("Arabidopsis thaliana", nrow(data2))
data2$Comparison <- rep("Arabidopsis-Zea", nrow(data2))
data2 <- data2[,c(3,4,1,2)]

## PhastCons

data2 <- data %>%
  group_by(Gene) %>%
  summarise(PhastCons_Element = paste(PhastCons_Element, collapse = ","))

data2$Species <- rep("Arabidopsis thaliana", nrow(data2))
data2$Comparison <- rep("Arabidopsis-Zea", nrow(data2))
data2 <- data2[,c(3,4,1,2)]

## CNEs

data2 <- At_Zm.CNEs %>%
  group_by(Gene) %>%
  summarise(CNEs_100kb = paste(CNEs_100kb, collapse = ","),
            CNEsNearest_Distance = unique(CNEsNearest_Distance),
            CNEsID = paste(CNEsID, collapse = ","))

data2$Species <- rep("Arabidopsis thaliana", nrow(data2))
data2$Comparison <- rep("Arabidopsis-Zea", nrow(data2))
data2 <- data2[,c(5,6,1,2,3,4)]

## GWAS

data2 <- data %>%
  group_by(Gene) %>%
  summarise(GWAS.SNPs_Presence = paste(GWAS.SNPs_Presence, collapse = ","), GWAS.SNPs_Location = paste(GWAS.SNPs_Location, collapse = ","))

data2$GWAS.SNPs_Location <- gsub(".:-1:-1", "", data2$GWAS.SNPs_Location)
data2$GWAS.SNPs_Location <- gsub(",", "", data2$GWAS.SNPs_Location)
data2$Species <- rep("Arabidopsis thaliana", nrow(data2))
data2$Comparison <- rep("Arabidopsis-Zea", nrow(data2))
data2 <- data2[,c(4,5,1,2,3)]

### Export & joint section ###

AtZmLECIF.table <- data2
AtZmLECIF.table.jic <- data2

AtZmPhyloP.table <- data2
AtZmPhyloP.table.jic <- data2

AtZmCS.table <- data2
AtZmCS.table.jic <- data2

AtZmPhastCons.table <- data2
AtZmPhastCons.table.jic <- data2

AtZmCNEs.table <- data2
AtZmCNEs.table.jic <- data2

AtZmGWAS.table <- data2
AtZmGWAS.table.jic <- data2

if(length(which((AtZmLECIF.table$Gene == AtZmPhyloP.table$Gene) == FALSE)) == 0){ print("Ready to join")}

AtZm.LECIF.PhyloP.table <- as_tibble(cbind(AtZmLECIF.table, AtZmPhyloP.table[, 4:6]))

if(length(which((AtZmCS.table$Gene == AtZm.LECIF.PhyloP.table$Gene) == FALSE)) == 0){ print("Ready to join")}

AtZm.LECIF.PhyloP.CS.table <- as_tibble(cbind(AtZm.LECIF.PhyloP.table, AtZmCS.table[,4]))

if(length(which((AtZmPhastCons.table$Gene == AtZm.LECIF.PhyloP.CS.table$Gene) == FALSE)) == 0){ print("Ready to join")}

AtZm.LECIF.PhyloP.CS.PhastCons.table <- as_tibble(cbind(AtZm.LECIF.PhyloP.CS.table, AtZmPhastCons.table[,4]))

which((AtZmCNEs.table$Gene %in% AtZmLECIF.table$Gene) == FALSE)
AtZmCNEs.table <- AtZmCNEs.table[-5445,]
if(length(which((AtZmCNEs.table$Gene == AtZm.LECIF.PhyloP.CS.PhastCons.table$Gene) == FALSE)) == 0){ print("Ready to join")}

AtZm.LECIF.PhyloP.CS.PhastCons.CNEs.table <- as_tibble(cbind(AtZm.LECIF.PhyloP.CS.PhastCons.table, AtZmCNEs.table[,4:6]))

if(length(which((AtZmGWAS.table$Gene == AtZm.LECIF.PhyloP.CS.PhastCons.CNEs.table$Gene) == FALSE)) == 0){ print("Ready to join")}

AtZm.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table <- as_tibble(cbind(AtZm.LECIF.PhyloP.CS.PhastCons.CNEs.table, AtZmGWAS.table[,4:5]))

## Custom table

data2 <- as_tibble(rbind(AtOs.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table,AtZm.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table))
data2 <- AtZm.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table
sticky_style <- list(backgroundColor = "#e8e8e8")
html.object <- htmltools::browsable(
  tagList(
    tags$button(
      tagList(fontawesome::fa("download"), "Download Table"),
      onclick = "Reactable.downloadDataCSV('GenesTableFUNCO', 'GenesTableFUNCO.csv')"
    ),
    
    reactable(data2, searchable = T, filterable = T, pagination = TRUE, paginationType = "jump", defaultPageSize = 100, defaultColDef = colDef(vAlign = "center", headerVAlign = "center", align = "center", minWidth = 80, footer = function(values, name) { htmltools::div(name, style = list(fontWeight = 600))}), 
              bordered = T, striped = TRUE, highlight = TRUE, resizable = TRUE, wrap = FALSE, sortable = TRUE, elementId = "GenesTableFUNCO", defaultSorted = list(LECIF_Mean = "desc"), height = 800,
              columns = list(
                Species = colDef(sticky = "left",
                                 style = sticky_style,
                                 headerStyle = sticky_style),
                Comparison = colDef(sticky = "left",
                                    style = sticky_style,
                                    headerStyle = sticky_style),
                Gene = colDef(sticky = "left",
                              style = sticky_style,
                              headerStyle = sticky_style),
                LECIF_Quartile = colDef(
                  style = function(value){
                    if(value == "Q1"){
                      color <- "#309143"
                    } else if(value == "Q2"){
                      color <- "#8ace7e"
                    } else if(value == "Q3"){
                      color <- "#ff684c"
                    } else if(value == "Q4"){
                      color <- "#b60a1c"
                    }
                    list(color = color, fontWeight = "bold")
                  }
                ),
                LECIF_Trend = colDef(cell = function(value, index) {
                  sparkline(data2$LECIF_Trend[[index]])
                }),
                PhyloP_Quartile = colDef(
                  style = function(value){
                    if(value == "Q1"){
                      color <- "#309143"
                    } else if(value == "Q2"){
                      color <- "#8ace7e"
                    } else if(value == "Q3"){
                      color <- "#ff684c"
                    } else if(value == "Q4"){
                      color <- "#b60a1c"
                    }
                    list(color = color, fontWeight = "bold")
                  }
                ),
                PhyloP_Trend = colDef(cell = function(value, index) {
                  sparkline(data2$PhyloP_Trend[[index]])
                }),
                PhastCons_Element = colDef(cell = function(value) {
                  if (length(grep("Yes", value)) > 0) "\u2714\ufe0f Yes" else "\u274c No"
                }),
                CNEs_100kb = colDef(cell = function(value) {
                  if (length(grep("Yes", value)) > 0) "\u2714\ufe0f Yes" else "\u274c No"
                }),
                GWAS.SNPs_Presence = colDef(cell = function(value) {
                  if (length(grep("Yes", value)) > 0) "\u2714\ufe0f Yes" else "\u274c No"
                })
              ))
  )
)

save_html(html.object, "H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Arabidopsis-Zea/SearchGenes_At-Zm.html")
save_html(html.object, "H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/SearchGenes_Arabidopsis.html")


#### 4.3 Oryza-Arabidopsis ####

### Format ###

## LECIF

Os_At.LECIF <- read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Oryza-Arabidopsis/Genes.Os-At.sorted.binned.LECIFscore.bed", header = F, sep = "\t")
Os_At.LECIF$V8[which(Os_At.LECIF$V8 == ".")] <- gsub(".", 0, fixed = T, Os_At.LECIF$V8[which(Os_At.LECIF$V8 == ".")])
Os_At.LECIF$V8 <- as.numeric(Os_At.LECIF$V8)
Os_At.LECIF <- Os_At.LECIF[,c(1,2,3,4,8)]

## PhyloP

Os_At.PhyloP <- read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Oryza-Arabidopsis/Genes.Os-At.sorted.binned.PhyloPscore.bed", header = F, sep = "\t")
Os_At.PhyloP$V8 <- as.numeric(Os_At.PhyloP$V8)
Os_At.PhyloP <- Os_At.PhyloP[,c(1,2,3,4,8)]
Os_At.PhyloP$V8[which(is.na(Os_At.PhyloP$V8) == TRUE)] <- 0

## hiHMM Chromatin States (CS)

Os_At.CS <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Oryza-Arabidopsis/Genes.Os-At.sorted.binned.CS.bed", header = F, sep = "\t")
Os_At.CS <- Os_At.CS[,c(1,2,3,4,8)]

## PhastCons

Os_At.PhastCons <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Oryza-Arabidopsis/Genes.Os-At.sorted.binned.PhastCons.bed", header = F, sep = "\t")
Os_At.PhastCons <- Os_At.PhastCons[,c(1,2,3,4,7)]

## CNEs 

Os_At.CNEs <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Oryza-Arabidopsis/Genes.Os-At.sorted.CNEs.bed", header = F, sep = "\t")
Os_At.CNEs$CNEs_100kb <- "No"
Os_At.CNEs$CNEs_100kb[which(Os_At.CNEs$V14 <= 100000)] <- "Yes"
Os_At.CNEs$CNEsNearest_Distance <- Os_At.CNEs$V14
Os_At.CNEs$CNEsID <- NA
Os_At.CNEs$CNEsID <- Os_At.CNEs$V8
Os_At.CNEs$Gene <- Os_At.CNEs$V4
Os_At.CNEs <- Os_At.CNEs[,c(18,17,15,16)]

## SNPs-GWAS

Os_At.GWAS <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Oryza-Arabidopsis/Genes.Os-At.sorted.binned.SNPsGWAS.bed", header = F, sep = "\t")
Os_At.GWAS$GWAS.SNPs_Presence <- Os_At.GWAS$V8
Os_At.GWAS$GWAS.SNPs_Location <- paste0(Os_At.GWAS$V5,":",Os_At.GWAS$V6,":",Os_At.GWAS$V7)
Os_At.GWAS <- Os_At.GWAS[,c(1,2,3,4,9,10)]

### First group by interval*gene and mean n LECIF hits ###

## LECIF

Os_At.LECIF$PASTED <- paste0(Os_At.LECIF$V1, ":", Os_At.LECIF$V2, ":", Os_At.LECIF$V3) 
data <- Os_At.LECIF %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4),LECIFscore = list(V8))
data$LECIFscoreM <- unlist(lapply(data$LECIFscore, FUN = mean), use.names = F)

## PhyloP

Os_At.PhyloP$PASTED <- paste0(Os_At.PhyloP$V1, ":", Os_At.PhyloP$V2, ":", Os_At.PhyloP$V3)
data <- Os_At.PhyloP %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4),PhyloPscore = list(V8))
data$PhyloPscoreM <- unlist(lapply(data$PhyloPscore, FUN = mean), use.names = F)

## hiHMM Chromatin States (CS)

Os_At.CS$PASTED <- paste0(Os_At.CS$V1, ":", Os_At.CS$V2, ":", Os_At.CS$V3)
data <- Os_At.CS %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4),CS_Composition = paste(V8, collapse = ","))

## PhastCons 

Os_At.PhastCons$PASTED <- paste0(Os_At.PhastCons$V1, ":", Os_At.PhastCons$V2, ":", Os_At.PhastCons$V3)
Os_At.PhastCons$V7 <- gsub("phastCons_predicted", "Yes", fixed = T, x = Os_At.PhastCons$V7)
Os_At.PhastCons$V7 <- gsub(".", "No", fixed = T, x = Os_At.PhastCons$V7)
data <- Os_At.PhastCons %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4), PhastCons_Element = paste(V7, collapse = ","))

## CNEs
## Already grouped by gene intervals

## GWAS

Os_At.GWAS$PASTED <- paste0(Os_At.GWAS$V1, ":", Os_At.GWAS$V2, ":", Os_At.GWAS$V3)
Os_At.GWAS$GWAS.SNPs_Presence <- gsub("0", "No", fixed = T, x = Os_At.GWAS$GWAS.SNPs_Presence)
Os_At.GWAS$GWAS.SNPs_Presence <- gsub("1", "Yes", fixed = T, x = Os_At.GWAS$GWAS.SNPs_Presence)
data <- Os_At.GWAS %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4), GWAS.SNPs_Presence = paste(GWAS.SNPs_Presence, collapse = ","), GWAS.SNPs_Location = paste(GWAS.SNPs_Location, collapse = ","))

### Second group by gene only and mean five bins ###

## LECIF

data2 <- data %>%
  group_by(Gene) %>%
  summarise(LECIF_Trend = list(LECIFscoreM))

data2$LECIF_Mean <- unlist(lapply(data2$LECIF_Trend, mean), use.names = F)
data2$LECIF_Mean <- round(data2$LECIF_Mean, digits = 3)
data2$LECIF_Quartile <- NA
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0 & percent_rank(x = data2$LECIF_Mean) < 0.25)] <- "Q4"
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0.25 & percent_rank(x = data2$LECIF_Mean) < 0.5)] <- "Q3"
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0.5 & percent_rank(x = data2$LECIF_Mean) < 0.75)] <- "Q2"
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0.75 & percent_rank(x = data2$LECIF_Mean) <= 1)] <- "Q1"

data2$Species <- rep("Oryza sativa", nrow(data2))
data2$Comparison <- rep("Oryza-Arabidopsis", nrow(data2))
data2 <- data2[,c(5,6,1,3,4,2)]

## PhyloP

data2 <- data %>%
  group_by(Gene) %>%
  summarise(PhyloP_Trend = list(PhyloPscoreM))

data2$PhyloP_Mean <- unlist(lapply(data2$PhyloP_Trend, mean), use.names = F)
data2$PhyloP_Mean <- round(data2$PhyloP_Mean, digits = 3)
data2$PhyloP_Quartile <- NA
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0 & percent_rank(x = data2$PhyloP_Mean) < 0.25)] <- "Q4"
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0.25 & percent_rank(x = data2$PhyloP_Mean) < 0.5)] <- "Q3"
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0.5 & percent_rank(x = data2$PhyloP_Mean) < 0.75)] <- "Q2"
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0.75 & percent_rank(x = data2$PhyloP_Mean) <= 1)] <- "Q1"

data2$Species <- rep("Oryza sativa", nrow(data2))
data2$Comparison <- rep("Oryza-Arabidopsis", nrow(data2))
data2 <- data2[,c(5,6,1,3,4,2)]

## hiHMM Chromatin States (CS)

data2 <- data %>%
  group_by(Gene) %>%
  summarise(CS_Composition = paste(CS_Composition, collapse = ","))

data2$Species <- rep("Oryza sativa", nrow(data2))
data2$Comparison <- rep("Oryza-Arabidopsis", nrow(data2))
data2 <- data2[,c(3,4,1,2)]

## PhastCons

data2 <- data %>%
  group_by(Gene) %>%
  summarise(PhastCons_Element = paste(PhastCons_Element, collapse = ","))

data2$Species <- rep("Oryza sativa", nrow(data2))
data2$Comparison <- rep("Oryza-Arabidopsis", nrow(data2))
data2 <- data2[,c(3,4,1,2)]

## CNEs

data2 <- Os_At.CNEs %>%
  group_by(Gene) %>%
  summarise(CNEs_100kb = paste(CNEs_100kb, collapse = ","),
            CNEsNearest_Distance = unique(CNEsNearest_Distance),
            CNEsID = paste(CNEsID, collapse = ","))

data2$Species <- rep("Oryza sativa", nrow(data2))
data2$Comparison <- rep("Oryza-Arabidopsis", nrow(data2))
data2 <- data2[,c(5,6,1,2,3,4)]

## GWAS

data2 <- data %>%
  group_by(Gene) %>%
  summarise(GWAS.SNPs_Presence = paste(GWAS.SNPs_Presence, collapse = ","), GWAS.SNPs_Location = paste(GWAS.SNPs_Location, collapse = ","))

data2$GWAS.SNPs_Location <- gsub(".:-1:-1", "", data2$GWAS.SNPs_Location)
data2$GWAS.SNPs_Location <- gsub(",", "", data2$GWAS.SNPs_Location)
data2$Species <- rep("Oryza sativa", nrow(data2))
data2$Comparison <- rep("Oryza-Arabidopsis", nrow(data2))
data2 <- data2[,c(4,5,1,2,3)]

### Export & joint section ###

OsAtLECIF.table <- data2
OsAtLECIF.table.jic <- data2

OsAtPhyloP.table <- data2
OsAtPhyloP.table.jic <- data2

OsAtCS.table <- data2
OsAtCS.table.jic <- data2

OsAtPhastCons.table <- data2
OsAtPhastCons.table.jic <- data2

OsAtCNEs.table <- data2
OsAtCNEs.table.jic <- data2

OsAtGWAS.table <- data2
OsAtGWAS.table.jic <- data2

if(length(which((OsAtLECIF.table$Gene == OsAtPhyloP.table$Gene) == FALSE)) == 0){ print("Ready to join")}

OsAt.LECIF.PhyloP.table <- as_tibble(cbind(OsAtLECIF.table, OsAtPhyloP.table[, 4:6]))

if(length(which((OsAtCS.table$Gene == OsAt.LECIF.PhyloP.table$Gene) == FALSE)) == 0){ print("Ready to join")}

OsAt.LECIF.PhyloP.CS.table <- as_tibble(cbind(OsAt.LECIF.PhyloP.table, OsAtCS.table[,4]))

if(length(which((OsAtPhastCons.table$Gene == OsAt.LECIF.PhyloP.CS.table$Gene) == FALSE)) == 0){ print("Ready to join")}

OsAt.LECIF.PhyloP.CS.PhastCons.table <- as_tibble(cbind(OsAt.LECIF.PhyloP.CS.table, OsAtPhastCons.table[,4]))

#OsAtCNEs.table <- OsAtCNEs.table[-5445,]
if(length(which((OsAtCNEs.table$Gene == OsAt.LECIF.PhyloP.CS.PhastCons.table$Gene) == FALSE)) == 0){ print("Ready to join")}

OsAt.LECIF.PhyloP.CS.PhastCons.CNEs.table <- as_tibble(cbind(OsAt.LECIF.PhyloP.CS.PhastCons.table, OsAtCNEs.table[,4:6]))

if(length(which((OsAtGWAS.table$Gene == OsAt.LECIF.PhyloP.CS.PhastCons.CNEs.table$Gene) == FALSE)) == 0){ print("Ready to join")}

OsAt.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table <- as_tibble(cbind(OsAt.LECIF.PhyloP.CS.PhastCons.CNEs.table, OsAtGWAS.table[,4:5]))

## Custom table

data2 <- OsAt.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table
sticky_style <- list(backgroundColor = "#e8e8e8")
html.object <- htmltools::browsable(
  tagList(
    tags$button(
      tagList(fontawesome::fa("download"), "Download Table"),
      onclick = "Reactable.downloadDataCSV('GenesTableFUNCO', 'GenesTableFUNCO.csv')"
    ),
    
    reactable(data2, searchable = T, filterable = T, pagination = TRUE, paginationType = "jump", defaultPageSize = 100, defaultColDef = colDef(vAlign = "center", headerVAlign = "center", align = "center", minWidth = 80, footer = function(values, name) { htmltools::div(name, style = list(fontWeight = 600))}), 
              bordered = T, striped = TRUE, highlight = TRUE, resizable = TRUE, wrap = FALSE, sortable = TRUE, elementId = "GenesTableFUNCO", defaultSorted = list(LECIF_Mean = "desc"), height = 800,
              columns = list(
                Species = colDef(sticky = "left",
                                 style = sticky_style,
                                 headerStyle = sticky_style),
                Comparison = colDef(sticky = "left",
                                    style = sticky_style,
                                    headerStyle = sticky_style),
                Gene = colDef(sticky = "left",
                              style = sticky_style,
                              headerStyle = sticky_style),
                LECIF_Quartile = colDef(
                  style = function(value){
                    if(value == "Q1"){
                      color <- "#309143"
                    } else if(value == "Q2"){
                      color <- "#8ace7e"
                    } else if(value == "Q3"){
                      color <- "#ff684c"
                    } else if(value == "Q4"){
                      color <- "#b60a1c"
                    }
                    list(color = color, fontWeight = "bold")
                  }
                ),
                LECIF_Trend = colDef(cell = function(value, index) {
                  sparkline(data2$LECIF_Trend[[index]])
                }),
                PhyloP_Quartile = colDef(
                  style = function(value){
                    if(value == "Q1"){
                      color <- "#309143"
                    } else if(value == "Q2"){
                      color <- "#8ace7e"
                    } else if(value == "Q3"){
                      color <- "#ff684c"
                    } else if(value == "Q4"){
                      color <- "#b60a1c"
                    }
                    list(color = color, fontWeight = "bold")
                  }
                ),
                PhyloP_Trend = colDef(cell = function(value, index) {
                  sparkline(data2$PhyloP_Trend[[index]])
                }),
                PhastCons_Element = colDef(cell = function(value) {
                  if (length(grep("Yes", value)) > 0) "\u2714\ufe0f Yes" else "\u274c No"
                }),
                CNEs_100kb = colDef(cell = function(value) {
                  if (length(grep("Yes", value)) > 0) "\u2714\ufe0f Yes" else "\u274c No"
                }),
                GWAS.SNPs_Presence = colDef(cell = function(value) {
                  if (length(grep("Yes", value)) > 0) "\u2714\ufe0f Yes" else "\u274c No"
                })
              ))
  )
)

save_html(html.object, "H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Oryza-Arabidopsis/SearchGenes_Os-At.html")

#### 4.4 Oryza-Zea ####

### Format ###

## LECIF

Os_Zm.LECIF <- read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Oryza-Zea/Genes.Os-Zm.sorted.binned.LECIFscore.bed", header = F, sep = "\t")
Os_Zm.LECIF$V8[which(Os_Zm.LECIF$V8 == ".")] <- gsub(".", 0, fixed = T, Os_Zm.LECIF$V8[which(Os_Zm.LECIF$V8 == ".")])
Os_Zm.LECIF$V8 <- as.numeric(Os_Zm.LECIF$V8)
Os_Zm.LECIF <- Os_Zm.LECIF[,c(1,2,3,4,8)]

## PhyloP

Os_Zm.PhyloP <- read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Oryza-Zea/Genes.Os-Zm.sorted.binned.PhyloPscore.bed", header = F, sep = "\t")
Os_Zm.PhyloP$V8 <- as.numeric(Os_Zm.PhyloP$V8)
Os_Zm.PhyloP <- Os_Zm.PhyloP[,c(1,2,3,4,8)]
Os_Zm.PhyloP$V8[which(is.na(Os_Zm.PhyloP$V8) == TRUE)] <- 0

## hiHMM Chromatin States (CS)

Os_Zm.CS <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Oryza-Zea/Genes.Os-Zm.sorted.binned.CS.bed", header = F, sep = "\t")
Os_Zm.CS <- Os_Zm.CS[,c(1,2,3,4,8)]

## PhastCons

Os_Zm.PhastCons <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Oryza-Zea/Genes.Os-Zm.sorted.binned.PhastCons.bed", header = F, sep = "\t")
Os_Zm.PhastCons <- Os_Zm.PhastCons[,c(1,2,3,4,7)]

## CNEs 

Os_Zm.CNEs <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Oryza-Zea/Genes.Os-Zm.sorted.CNEs.bed", header = F, sep = "\t")
Os_Zm.CNEs$CNEs_100kb <- "No"
Os_Zm.CNEs$CNEs_100kb[which(Os_Zm.CNEs$V14 <= 100000)] <- "Yes"
Os_Zm.CNEs$CNEsNearest_Distance <- Os_Zm.CNEs$V14
Os_Zm.CNEs$CNEsID <- NA
Os_Zm.CNEs$CNEsID <- Os_Zm.CNEs$V8
Os_Zm.CNEs$Gene <- Os_Zm.CNEs$V4
Os_Zm.CNEs <- Os_Zm.CNEs[,c(18,17,15,16)]

## GWAS

Os_Zm.GWAS <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Oryza-Zea/Genes.Os-Zm.sorted.binned.SNPsGWAS.bed", header = F, sep = "\t")
Os_Zm.GWAS$GWAS.SNPs_Presence <- Os_Zm.GWAS$V8
Os_Zm.GWAS$GWAS.SNPs_Location <- paste0(Os_Zm.GWAS$V5,":",Os_Zm.GWAS$V6,":",Os_Zm.GWAS$V7)
Os_Zm.GWAS <- Os_Zm.GWAS[,c(1,2,3,4,9,10)]

### First group by interval*gene and mean n LECIF hits ###

## LECIF

Os_Zm.LECIF$PASTED <- paste0(Os_Zm.LECIF$V1, ":", Os_Zm.LECIF$V2, ":", Os_Zm.LECIF$V3) 
data <- Os_Zm.LECIF %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4),LECIFscore = list(V8))
data$LECIFscoreM <- unlist(lapply(data$LECIFscore, FUN = mean), use.names = F)

## PhyloP

Os_Zm.PhyloP$PASTED <- paste0(Os_Zm.PhyloP$V1, ":", Os_Zm.PhyloP$V2, ":", Os_Zm.PhyloP$V3)
data <- Os_Zm.PhyloP %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4),PhyloPscore = list(V8))
data$PhyloPscoreM <- unlist(lapply(data$PhyloPscore, FUN = mean), use.names = F)

## hiHMM Chromatin States (CS)

Os_Zm.CS$PASTED <- paste0(Os_Zm.CS$V1, ":", Os_Zm.CS$V2, ":", Os_Zm.CS$V3)
data <- Os_Zm.CS %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4),CS_Composition = paste(V8, collapse = ","))

## PhastCons 

Os_Zm.PhastCons$PASTED <- paste0(Os_Zm.PhastCons$V1, ":", Os_Zm.PhastCons$V2, ":", Os_Zm.PhastCons$V3)
Os_Zm.PhastCons$V7 <- gsub("phastCons_predicted", "Yes", fixed = T, x = Os_Zm.PhastCons$V7)
Os_Zm.PhastCons$V7 <- gsub(".", "No", fixed = T, x = Os_Zm.PhastCons$V7)
data <- Os_Zm.PhastCons %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4), PhastCons_Element = paste(V7, collapse = ","))

## CNEs
## Already grouped by gene intervals

## GWAS

Os_Zm.GWAS$PASTED <- paste0(Os_Zm.GWAS$V1, ":", Os_Zm.GWAS$V2, ":", Os_Zm.GWAS$V3)
Os_Zm.GWAS$GWAS.SNPs_Presence <- gsub("0", "No", fixed = T, x = Os_Zm.GWAS$GWAS.SNPs_Presence)
Os_Zm.GWAS$GWAS.SNPs_Presence <- gsub("1", "Yes", fixed = T, x = Os_Zm.GWAS$GWAS.SNPs_Presence)
data <- Os_Zm.GWAS %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4), GWAS.SNPs_Presence = paste(GWAS.SNPs_Presence, collapse = ","), GWAS.SNPs_Location = paste(GWAS.SNPs_Location, collapse = ","))

### Second group by gene only and mean five bins ###

## LECIF

data2 <- data %>%
  group_by(Gene) %>%
  summarise(LECIF_Trend = list(LECIFscoreM))

data2$LECIF_Mean <- unlist(lapply(data2$LECIF_Trend, mean), use.names = F)
data2$LECIF_Mean <- round(data2$LECIF_Mean, digits = 3)
data2$LECIF_Quartile <- NA
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0 & percent_rank(x = data2$LECIF_Mean) < 0.25)] <- "Q4"
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0.25 & percent_rank(x = data2$LECIF_Mean) < 0.5)] <- "Q3"
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0.5 & percent_rank(x = data2$LECIF_Mean) < 0.75)] <- "Q2"
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0.75 & percent_rank(x = data2$LECIF_Mean) <= 1)] <- "Q1"

data2$Species <- rep("Oryza sativa", nrow(data2))
data2$Comparison <- rep("Oryza-Zea", nrow(data2))
data2 <- data2[,c(5,6,1,3,4,2)]

## PhyloP

data2 <- data %>%
  group_by(Gene) %>%
  summarise(PhyloP_Trend = list(PhyloPscoreM))

data2$PhyloP_Mean <- unlist(lapply(data2$PhyloP_Trend, mean), use.names = F)
data2$PhyloP_Mean <- round(data2$PhyloP_Mean, digits = 3)
data2$PhyloP_Quartile <- NA
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0 & percent_rank(x = data2$PhyloP_Mean) < 0.25)] <- "Q4"
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0.25 & percent_rank(x = data2$PhyloP_Mean) < 0.5)] <- "Q3"
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0.5 & percent_rank(x = data2$PhyloP_Mean) < 0.75)] <- "Q2"
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0.75 & percent_rank(x = data2$PhyloP_Mean) <= 1)] <- "Q1"

data2$Species <- rep("Oryza sativa", nrow(data2))
data2$Comparison <- rep("Oryza-Zea", nrow(data2))
data2 <- data2[,c(5,6,1,3,4,2)]

## hiHMM Chromatin States (CS)

data2 <- data %>%
  group_by(Gene) %>%
  summarise(CS_Composition = paste(CS_Composition, collapse = ","))

data2$Species <- rep("Oryza sativa", nrow(data2))
data2$Comparison <- rep("Oryza-Zea", nrow(data2))
data2 <- data2[,c(3,4,1,2)]

## PhastCons

data2 <- data %>%
  group_by(Gene) %>%
  summarise(PhastCons_Element = paste(PhastCons_Element, collapse = ","))

data2$Species <- rep("Oryza sativa", nrow(data2))
data2$Comparison <- rep("Oryza-Zea", nrow(data2))
data2 <- data2[,c(3,4,1,2)]

## CNEs

data2 <- Os_Zm.CNEs %>%
  group_by(Gene) %>%
  summarise(CNEs_100kb = paste(CNEs_100kb, collapse = ","),
            CNEsNearest_Distance = unique(CNEsNearest_Distance),
            CNEsID = paste(CNEsID, collapse = ","))

data2$Species <- rep("Oryza sativa", nrow(data2))
data2$Comparison <- rep("Oryza-Zea", nrow(data2))
data2 <- data2[,c(5,6,1,2,3,4)]

## GWAS

data2 <- data %>%
  group_by(Gene) %>%
  summarise(GWAS.SNPs_Presence = paste(GWAS.SNPs_Presence, collapse = ","), GWAS.SNPs_Location = paste(GWAS.SNPs_Location, collapse = ","))

data2$GWAS.SNPs_Location <- gsub(".:-1:-1", "", data2$GWAS.SNPs_Location)
data2$GWAS.SNPs_Location <- gsub(",", "", data2$GWAS.SNPs_Location)
data2$Species <- rep("Oryza sativa", nrow(data2))
data2$Comparison <- rep("Oryza-Zea", nrow(data2))
data2 <- data2[,c(4,5,1,2,3)]

### Export & joint section ###

OsZmLECIF.table <- data2
OsZmLECIF.table.jic <- data2

OsZmPhyloP.table <- data2
OsZmPhyloP.table.jic <- data2

OsZmCS.table <- data2
OsZmCS.table.jic <- data2

OsZmPhastCons.table <- data2
OsZmPhastCons.table.jic <- data2

OsZmCNEs.table <- data2
OsZmCNEs.table.jic <- data2

OsZmGWAS.table <- data2
OsZmGWAS.table.jic <- data2

if(length(which((OsZmLECIF.table$Gene == OsZmPhyloP.table$Gene) == FALSE)) == 0){ print("Ready to join")}

OsZm.LECIF.PhyloP.table <- as_tibble(cbind(OsZmLECIF.table, OsZmPhyloP.table[, 4:6]))

if(length(which((OsZmCS.table$Gene == OsZm.LECIF.PhyloP.table$Gene) == FALSE)) == 0){ print("Ready to join")}

OsZm.LECIF.PhyloP.CS.table <- as_tibble(cbind(OsZm.LECIF.PhyloP.table, OsZmCS.table[,4]))

if(length(which((OsZmPhastCons.table$Gene == OsZm.LECIF.PhyloP.CS.table$Gene) == FALSE)) == 0){ print("Ready to join")}

OsZm.LECIF.PhyloP.CS.PhastCons.table <- as_tibble(cbind(OsZm.LECIF.PhyloP.CS.table, OsZmPhastCons.table[,4]))

which((OsZmCNEs.table$Gene %in% OsZmLECIF.table$Gene) == FALSE)
#OsZmCNEs.table <- OsZmCNEs.table[-5445,]
if(length(which((OsZmCNEs.table$Gene == OsZm.LECIF.PhyloP.CS.PhastCons.table$Gene) == FALSE)) == 0){ print("Ready to join")}

OsZm.LECIF.PhyloP.CS.PhastCons.CNEs.table <- as_tibble(cbind(OsZm.LECIF.PhyloP.CS.PhastCons.table, OsZmCNEs.table[,4:6]))

if(length(which((OsZmGWAS.table$Gene == OsZm.LECIF.PhyloP.CS.PhastCons.CNEs.table$Gene) == FALSE)) == 0){ print("Ready to join")}

OsZm.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table <- as_tibble(cbind(OsZm.LECIF.PhyloP.CS.PhastCons.CNEs.table, OsZmGWAS.table[,4:5]))

## Custom table

data2 <- as_tibble(rbind(OsAt.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table,OsZm.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table, AtOs.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table, AtZm.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table))
data2 <- OsZm.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table
sticky_style <- list(backgroundColor = "#e8e8e8")
html.object <- htmltools::browsable(
  tagList(
    tags$button(
      tagList(fontawesome::fa("download"), "Download Table"),
      onclick = "Reactable.downloadDataCSV('GenesTableFUNCO', 'GenesTableFUNCO.csv')"
    ),
    
    reactable(data2, searchable = T, filterable = T, pagination = TRUE, paginationType = "jump", defaultPageSize = 100, defaultColDef = colDef(vAlign = "center", headerVAlign = "center", align = "center", minWidth = 80, footer = function(values, name) { htmltools::div(name, style = list(fontWeight = 600))}), 
              bordered = T, striped = TRUE, highlight = TRUE, resizable = TRUE, wrap = FALSE, sortable = TRUE, elementId = "GenesTableFUNCO", defaultSorted = list(LECIF_Mean = "desc"), height = 800,
              columns = list(
                Species = colDef(sticky = "left",
                                 style = sticky_style,
                                 headerStyle = sticky_style),
                Comparison = colDef(sticky = "left",
                                    style = sticky_style,
                                    headerStyle = sticky_style),
                Gene = colDef(sticky = "left",
                              style = sticky_style,
                              headerStyle = sticky_style),
                LECIF_Quartile = colDef(
                  style = function(value){
                    if(value == "Q1"){
                      color <- "#309143"
                    } else if(value == "Q2"){
                      color <- "#8ace7e"
                    } else if(value == "Q3"){
                      color <- "#ff684c"
                    } else if(value == "Q4"){
                      color <- "#b60a1c"
                    }
                    list(color = color, fontWeight = "bold")
                  }
                ),
                LECIF_Trend = colDef(cell = function(value, index) {
                  sparkline(data2$LECIF_Trend[[index]])
                }),
                PhyloP_Quartile = colDef(
                  style = function(value){
                    if(value == "Q1"){
                      color <- "#309143"
                    } else if(value == "Q2"){
                      color <- "#8ace7e"
                    } else if(value == "Q3"){
                      color <- "#ff684c"
                    } else if(value == "Q4"){
                      color <- "#b60a1c"
                    }
                    list(color = color, fontWeight = "bold")
                  }
                ),
                PhyloP_Trend = colDef(cell = function(value, index) {
                  sparkline(data2$PhyloP_Trend[[index]])
                }),
                PhastCons_Element = colDef(cell = function(value) {
                  if (length(grep("Yes", value)) > 0) "\u2714\ufe0f Yes" else "\u274c No"
                }),
                CNEs_100kb = colDef(cell = function(value) {
                  if (length(grep("Yes", value)) > 0) "\u2714\ufe0f Yes" else "\u274c No"
                }),
                GWAS.SNPs_Presence = colDef(cell = function(value) {
                  if (length(grep("Yes", value)) > 0) "\u2714\ufe0f Yes" else "\u274c No"
                })
              ))
  )
)

save_html(html.object, "H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Oryza-Zea/SearchGenes_Os-Zm.html")

#### 4.5 Zea-Arabidopsis ####

### Format ###

## LECIF

Zm_At.LECIF <- read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Zea-Arabidopsis/Genes.Zm-At.sorted.binned.LECIFscore.bed", header = F, sep = "\t")
Zm_At.LECIF$V8[which(Zm_At.LECIF$V8 == ".")] <- gsub(".", 0, fixed = T, Zm_At.LECIF$V8[which(Zm_At.LECIF$V8 == ".")])
Zm_At.LECIF$V8 <- as.numeric(Zm_At.LECIF$V8)
Zm_At.LECIF <- Zm_At.LECIF[,c(1,2,3,4,8)]

## PhyloP

Zm_At.PhyloP <- read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Zea-Arabidopsis/Genes.Zm-At.sorted.binned.PhyloPscore.bed", header = F, sep = "\t")
Zm_At.PhyloP$V8 <- as.numeric(Zm_At.PhyloP$V8)
Zm_At.PhyloP <- Zm_At.PhyloP[,c(1,2,3,4,8)]
Zm_At.PhyloP$V8[which(is.na(Zm_At.PhyloP$V8) == TRUE)] <- 0

## hiHMM Chromatin States (CS)

Zm_At.CS <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Zea-Arabidopsis/Genes.Zm-At.sorted.binned.CS.bed", header = F, sep = "\t")
Zm_At.CS <- Zm_At.CS[,c(1,2,3,4,8)]

## PhastCons

Zm_At.PhastCons <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Zea-Arabidopsis/Genes.Zm-At.sorted.binned.PhastCons.bed", header = F, sep = "\t")
Zm_At.PhastCons <- Zm_At.PhastCons[,c(1,2,3,4,7)]

## CNEs 

Zm_At.CNEs <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Zea-Arabidopsis/Genes.Zm-At.sorted.CNEs.bed", header = F, sep = "\t")
Zm_At.CNEs$CNEs_100kb <- "No"
Zm_At.CNEs$CNEs_100kb[which(Zm_At.CNEs$V14 <= 100000)] <- "Yes"
Zm_At.CNEs$CNEsNearest_Distance <- Zm_At.CNEs$V14
Zm_At.CNEs$CNEsID <- NA
Zm_At.CNEs$CNEsID <- Zm_At.CNEs$V8
Zm_At.CNEs$Gene <- Zm_At.CNEs$V4
Zm_At.CNEs <- Zm_At.CNEs[,c(18,17,15,16)]

## SNPs-GWAS

Zm_At.GWAS <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Zea-Arabidopsis/Genes.Zm-At.sorted.SNPsGWAS.bed", header = F, sep = "\t")
Zm_At.GWAS$GWAS.SNPs_Presence <- Zm_At.GWAS$V8
Zm_At.GWAS$GWAS.SNPs_Location <- paste0(Zm_At.GWAS$V5,":",Zm_At.GWAS$V6,":",Zm_At.GWAS$V7)
Zm_At.GWAS <- Zm_At.GWAS[,c(1,2,3,4,9,10)]

### First group by interval*gene and mean n LECIF hits ###

## LECIF

Zm_At.LECIF$PASTED <- paste0(Zm_At.LECIF$V1, ":", Zm_At.LECIF$V2, ":", Zm_At.LECIF$V3) 
data <- Zm_At.LECIF %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4),LECIFscore = list(V8))
data$LECIFscoreM <- unlist(lapply(data$LECIFscore, FUN = mean), use.names = F)

## PhyloP

Zm_At.PhyloP$PASTED <- paste0(Zm_At.PhyloP$V1, ":", Zm_At.PhyloP$V2, ":", Zm_At.PhyloP$V3)
data <- Zm_At.PhyloP %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4),PhyloPscore = list(V8))
data$PhyloPscoreM <- unlist(lapply(data$PhyloPscore, FUN = mean), use.names = F)

## hiHMM Chromatin States (CS)

Zm_At.CS$PASTED <- paste0(Zm_At.CS$V1, ":", Zm_At.CS$V2, ":", Zm_At.CS$V3)
data <- Zm_At.CS %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4),CS_Composition = paste(V8, collapse = ","))

## PhastCons 

Zm_At.PhastCons$PASTED <- paste0(Zm_At.PhastCons$V1, ":", Zm_At.PhastCons$V2, ":", Zm_At.PhastCons$V3)
Zm_At.PhastCons$V7 <- gsub("phastCons_predicted", "Yes", fixed = T, x = Zm_At.PhastCons$V7)
Zm_At.PhastCons$V7 <- gsub(".", "No", fixed = T, x = Zm_At.PhastCons$V7)
data <- Zm_At.PhastCons %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4), PhastCons_Element = paste(V7, collapse = ","))

## CNEs
## Already grouped by gene intervals

## GWAS
## Already grouped by gene intervals

### Second group by gene only and mean five bins ###

## LECIF

data2 <- data %>%
  group_by(Gene) %>%
  summarise(LECIF_Trend = list(LECIFscoreM))

data2$LECIF_Mean <- unlist(lapply(data2$LECIF_Trend, mean), use.names = F)
data2$LECIF_Mean <- round(data2$LECIF_Mean, digits = 3)
data2$LECIF_Quartile <- NA
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0 & percent_rank(x = data2$LECIF_Mean) < 0.25)] <- "Q4"
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0.25 & percent_rank(x = data2$LECIF_Mean) < 0.5)] <- "Q3"
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0.5 & percent_rank(x = data2$LECIF_Mean) < 0.75)] <- "Q2"
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0.75 & percent_rank(x = data2$LECIF_Mean) <= 1)] <- "Q1"

data2$Species <- rep("Zea mays", nrow(data2))
data2$Comparison <- rep("Zea-Arabidopsis", nrow(data2))
data2 <- data2[,c(5,6,1,3,4,2)]

## PhyloP

data2 <- data %>%
  group_by(Gene) %>%
  summarise(PhyloP_Trend = list(PhyloPscoreM))

data2$PhyloP_Mean <- unlist(lapply(data2$PhyloP_Trend, mean), use.names = F)
data2$PhyloP_Mean <- round(data2$PhyloP_Mean, digits = 3)
data2$PhyloP_Quartile <- NA
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0 & percent_rank(x = data2$PhyloP_Mean) < 0.25)] <- "Q4"
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0.25 & percent_rank(x = data2$PhyloP_Mean) < 0.5)] <- "Q3"
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0.5 & percent_rank(x = data2$PhyloP_Mean) < 0.75)] <- "Q2"
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0.75 & percent_rank(x = data2$PhyloP_Mean) <= 1)] <- "Q1"

data2$Species <- rep("Zea mays", nrow(data2))
data2$Comparison <- rep("Zea-Arabidopsis", nrow(data2))
data2 <- data2[,c(5,6,1,3,4,2)]

## hiHMM Chromatin States (CS)

data2 <- data %>%
  group_by(Gene) %>%
  summarise(CS_Composition = paste(CS_Composition, collapse = ","))

data2$Species <- rep("Zea mays", nrow(data2))
data2$Comparison <- rep("Zea-Arabidopsis", nrow(data2))
data2 <- data2[,c(3,4,1,2)]

## PhastCons

data2 <- data %>%
  group_by(Gene) %>%
  summarise(PhastCons_Element = paste(PhastCons_Element, collapse = ","))

data2$Species <- rep("Zea mays", nrow(data2))
data2$Comparison <- rep("Zea-Arabidopsis", nrow(data2))
data2 <- data2[,c(3,4,1,2)]

## CNEs

data2 <- Zm_At.CNEs %>%
  group_by(Gene) %>%
  summarise(CNEs_100kb = paste(CNEs_100kb, collapse = ","),
            CNEsNearest_Distance = unique(CNEsNearest_Distance),
            CNEsID = paste(CNEsID, collapse = ","))

data2$Species <- rep("Zea mays", nrow(data2))
data2$Comparison <- rep("Zea-Arabidopsis", nrow(data2))
data2 <- data2[,c(5,6,1,2,3,4)]

## GWAS

Zm_At.GWAS$PASTED <- paste0(Zm_At.GWAS$V1, ":", Zm_At.GWAS$V2, ":", Zm_At.GWAS$V3)
Zm_At.GWAS$GWAS.SNPs_Presence <- gsub("0", "No", fixed = T, x = Zm_At.GWAS$GWAS.SNPs_Presence)
Zm_At.GWAS$GWAS.SNPs_Presence <- gsub("1", "Yes", fixed = T, x = Zm_At.GWAS$GWAS.SNPs_Presence)
colnames(Zm_At.GWAS) <- c("V1", "V2", "V3", "Gene", "GWAS.SNPs_Presence", "GWAS.SNPs_Location", "PASTED")
data2 <- Zm_At.GWAS %>%
  group_by(Gene) %>%
  summarise(GWAS.SNPs_Presence = paste(GWAS.SNPs_Presence, collapse = ","), GWAS.SNPs_Location = paste(GWAS.SNPs_Location, collapse = ","))

data2$GWAS.SNPs_Location <- gsub(".:-1:-1", "", data2$GWAS.SNPs_Location)
data2$GWAS.SNPs_Location <- gsub(",", "", data2$GWAS.SNPs_Location)
data2$Species <- rep("Zea mays", nrow(data2))
data2$Comparison <- rep("Zea-Arabidopsis", nrow(data2))
data2 <- data2[,c(4,5,1,2,3)]

### Export & joint section ###

ZmAtLECIF.table <- data2
ZmAtLECIF.table.jic <- data2

ZmAtPhyloP.table <- data2
ZmAtPhyloP.table.jic <- data2

ZmAtCS.table <- data2
ZmAtCS.table.jic <- data2

ZmAtPhastCons.table <- data2
ZmAtPhastCons.table.jic <- data2

ZmAtCNEs.table <- data2
ZmAtCNEs.table.jic <- data2

ZmAtGWAS.table <- data2
ZmAtGWAS.table.jic <- data2

if(length(which((ZmAtLECIF.table$Gene == ZmAtPhyloP.table$Gene) == FALSE)) == 0){ print("Ready to join")}

ZmAt.LECIF.PhyloP.table <- as_tibble(cbind(ZmAtLECIF.table, ZmAtPhyloP.table[, 4:6]))

if(length(which((ZmAtCS.table$Gene == ZmAt.LECIF.PhyloP.table$Gene) == FALSE)) == 0){ print("Ready to join")}

ZmAt.LECIF.PhyloP.CS.table <- as_tibble(cbind(ZmAt.LECIF.PhyloP.table, ZmAtCS.table[,4]))

if(length(which((ZmAtPhastCons.table$Gene == ZmAt.LECIF.PhyloP.CS.table$Gene) == FALSE)) == 0){ print("Ready to join")}

ZmAt.LECIF.PhyloP.CS.PhastCons.table <- as_tibble(cbind(ZmAt.LECIF.PhyloP.CS.table, ZmAtPhastCons.table[,4]))

#ZmAtCNEs.table <- ZmAtCNEs.table[-5445,]
if(length(which((ZmAtCNEs.table$Gene == ZmAt.LECIF.PhyloP.CS.PhastCons.table$Gene) == FALSE)) == 0){ print("Ready to join")}

ZmAt.LECIF.PhyloP.CS.PhastCons.CNEs.table <- as_tibble(cbind(ZmAt.LECIF.PhyloP.CS.PhastCons.table, ZmAtCNEs.table[,4:6]))

if(length(which((ZmAtGWAS.table$Gene == ZmAt.LECIF.PhyloP.CS.PhastCons.CNEs.table$Gene) == FALSE)) == 0){ print("Ready to join")}

ZmAt.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table <- as_tibble(cbind(ZmAt.LECIF.PhyloP.CS.PhastCons.CNEs.table, ZmAtGWAS.table[,4:5]))

## Custom table

data2 <- ZmAt.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table
sticky_style <- list(backgroundColor = "#e8e8e8")
html.object <- htmltools::browsable(
  tagList(
    tags$button(
      tagList(fontawesome::fa("download"), "Download Table"),
      onclick = "Reactable.downloadDataCSV('GenesTableFUNCO', 'GenesTableFUNCO.csv')"
    ),
    
    reactable(data2, searchable = T, filterable = T, pagination = TRUE, paginationType = "jump", defaultPageSize = 100, defaultColDef = colDef(vAlign = "center", headerVAlign = "center", align = "center", minWidth = 80, footer = function(values, name) { htmltools::div(name, style = list(fontWeight = 600))}), 
              bordered = T, striped = TRUE, highlight = TRUE, resizable = TRUE, wrap = FALSE, sortable = TRUE, elementId = "GenesTableFUNCO", defaultSorted = list(LECIF_Mean = "desc"), height = 800,
              columns = list(
                Species = colDef(sticky = "left",
                                 style = sticky_style,
                                 headerStyle = sticky_style),
                Comparison = colDef(sticky = "left",
                                    style = sticky_style,
                                    headerStyle = sticky_style),
                Gene = colDef(sticky = "left",
                              style = sticky_style,
                              headerStyle = sticky_style),
                LECIF_Quartile = colDef(
                  style = function(value){
                    if(value == "Q1"){
                      color <- "#309143"
                    } else if(value == "Q2"){
                      color <- "#8ace7e"
                    } else if(value == "Q3"){
                      color <- "#ff684c"
                    } else if(value == "Q4"){
                      color <- "#b60a1c"
                    }
                    list(color = color, fontWeight = "bold")
                  }
                ),
                LECIF_Trend = colDef(cell = function(value, index) {
                  sparkline(data2$LECIF_Trend[[index]])
                }),
                PhyloP_Quartile = colDef(
                  style = function(value){
                    if(value == "Q1"){
                      color <- "#309143"
                    } else if(value == "Q2"){
                      color <- "#8ace7e"
                    } else if(value == "Q3"){
                      color <- "#ff684c"
                    } else if(value == "Q4"){
                      color <- "#b60a1c"
                    }
                    list(color = color, fontWeight = "bold")
                  }
                ),
                PhyloP_Trend = colDef(cell = function(value, index) {
                  sparkline(data2$PhyloP_Trend[[index]])
                }),
                PhastCons_Element = colDef(cell = function(value) {
                  if (length(grep("Yes", value)) > 0) "\u2714\ufe0f Yes" else "\u274c No"
                }),
                CNEs_100kb = colDef(cell = function(value) {
                  if (length(grep("Yes", value)) > 0) "\u2714\ufe0f Yes" else "\u274c No"
                }),
                GWAS.SNPs_Presence = colDef(cell = function(value) {
                  if (length(grep("Yes", value)) > 0) "\u2714\ufe0f Yes" else "\u274c No"
                })
              ))
  )
)

save_html(html.object, "H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Zea-Arabidopsis/SearchGenes_Zm-At.html")

#### 4.6 Zea-Oryza ####

### Format ###

## LECIF

Zm_Os.LECIF <- read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Zea-Oryza/Genes.Zm-Os.sorted.binned.LECIFscore.bed", header = F, sep = "\t")
Zm_Os.LECIF$V8[which(Zm_Os.LECIF$V8 == ".")] <- gsub(".", 0, fixed = T, Zm_Os.LECIF$V8[which(Zm_Os.LECIF$V8 == ".")])
Zm_Os.LECIF$V8 <- as.numeric(Zm_Os.LECIF$V8)
Zm_Os.LECIF <- Zm_Os.LECIF[,c(1,2,3,4,8)]

## PhyloP

Zm_Os.PhyloP <- read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Zea-Oryza/Genes.Zm-Os.sorted.binned.PhyloPscore.bed", header = F, sep = "\t")
Zm_Os.PhyloP$V8 <- as.numeric(Zm_Os.PhyloP$V8)
Zm_Os.PhyloP <- Zm_Os.PhyloP[,c(1,2,3,4,8)]
Zm_Os.PhyloP$V8[which(is.na(Zm_Os.PhyloP$V8) == TRUE)] <- 0

## hiHMM Chromatin States (CS)

Zm_Os.CS <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Zea-Oryza/Genes.Zm-Os.sorted.binned.CS.bed", header = F, sep = "\t")
Zm_Os.CS <- Zm_Os.CS[,c(1,2,3,4,8)]

## PhastCons

Zm_Os.PhastCons <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Zea-Oryza/Genes.Zm-Os.sorted.binned.PhastCons.bed", header = F, sep = "\t")
Zm_Os.PhastCons <- Zm_Os.PhastCons[,c(1,2,3,4,7)]

## CNEs 

Zm_Os.CNEs <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Zea-Oryza/Genes.Zm-Os.sorted.CNEs.bed", header = F, sep = "\t")
Zm_Os.CNEs$CNEs_100kb <- "No"
Zm_Os.CNEs$CNEs_100kb[which(Zm_Os.CNEs$V14 <= 100000)] <- "Yes"
Zm_Os.CNEs$CNEsNearest_Distance <- Zm_Os.CNEs$V14
Zm_Os.CNEs$CNEsID <- NA
Zm_Os.CNEs$CNEsID <- Zm_Os.CNEs$V8
Zm_Os.CNEs$Gene <- Zm_Os.CNEs$V4
Zm_Os.CNEs <- Zm_Os.CNEs[,c(18,17,15,16)]

## GWAS

Zm_Os.GWAS <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Zea-Oryza/Genes.Zm-Os.sorted.SNPsGWAS.bed", header = F, sep = "\t")
Zm_Os.GWAS$GWAS.SNPs_Presence <- Zm_Os.GWAS$V8
Zm_Os.GWAS$GWAS.SNPs_Location <- paste0(Zm_Os.GWAS$V5,":",Zm_Os.GWAS$V6,":",Zm_Os.GWAS$V7)
Zm_Os.GWAS <- Zm_Os.GWAS[,c(1,2,3,4,9,10)]

### First group by interval*gene and mean n LECIF hits ###

## LECIF

Zm_Os.LECIF$PASTED <- paste0(Zm_Os.LECIF$V1, ":", Zm_Os.LECIF$V2, ":", Zm_Os.LECIF$V3) 
data <- Zm_Os.LECIF %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4),LECIFscore = list(V8))
data$LECIFscoreM <- unlist(lapply(data$LECIFscore, FUN = mean), use.names = F)

## PhyloP

Zm_Os.PhyloP$PASTED <- paste0(Zm_Os.PhyloP$V1, ":", Zm_Os.PhyloP$V2, ":", Zm_Os.PhyloP$V3)
data <- Zm_Os.PhyloP %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4),PhyloPscore = list(V8))
data$PhyloPscoreM <- unlist(lapply(data$PhyloPscore, FUN = mean), use.names = F)

## hiHMM Chromatin States (CS)

Zm_Os.CS$PASTED <- paste0(Zm_Os.CS$V1, ":", Zm_Os.CS$V2, ":", Zm_Os.CS$V3)
data <- Zm_Os.CS %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4),CS_Composition = paste(V8, collapse = ","))

## PhastCons 

Zm_Os.PhastCons$PASTED <- paste0(Zm_Os.PhastCons$V1, ":", Zm_Os.PhastCons$V2, ":", Zm_Os.PhastCons$V3)
Zm_Os.PhastCons$V7 <- gsub("phastCons_predicted", "Yes", fixed = T, x = Zm_Os.PhastCons$V7)
Zm_Os.PhastCons$V7 <- gsub(".", "No", fixed = T, x = Zm_Os.PhastCons$V7)
data <- Zm_Os.PhastCons %>%
  group_by(PASTED) %>%
  summarise(Gene = unique(V4), PhastCons_Element = paste(V7, collapse = ","))

## CNEs
## Already grouped by gene intervals

## GWAS
## Already grouped by gene intervals

### Second group by gene only and mean five bins ###

## LECIF

data2 <- data %>%
  group_by(Gene) %>%
  summarise(LECIF_Trend = list(LECIFscoreM))

data2$LECIF_Mean <- unlist(lapply(data2$LECIF_Trend, mean), use.names = F)
data2$LECIF_Mean <- round(data2$LECIF_Mean, digits = 3)
data2$LECIF_Quartile <- NA
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0 & percent_rank(x = data2$LECIF_Mean) < 0.25)] <- "Q4"
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0.25 & percent_rank(x = data2$LECIF_Mean) < 0.5)] <- "Q3"
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0.5 & percent_rank(x = data2$LECIF_Mean) < 0.75)] <- "Q2"
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0.75 & percent_rank(x = data2$LECIF_Mean) <= 1)] <- "Q1"

data2$Species <- rep("Zea mays", nrow(data2))
data2$Comparison <- rep("Zea-Oryza", nrow(data2))
data2 <- data2[,c(5,6,1,3,4,2)]

## PhyloP

data2 <- data %>%
  group_by(Gene) %>%
  summarise(PhyloP_Trend = list(PhyloPscoreM))

data2$PhyloP_Mean <- unlist(lapply(data2$PhyloP_Trend, mean), use.names = F)
data2$PhyloP_Mean <- round(data2$PhyloP_Mean, digits = 3)
data2$PhyloP_Quartile <- NA
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0 & percent_rank(x = data2$PhyloP_Mean) < 0.25)] <- "Q4"
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0.25 & percent_rank(x = data2$PhyloP_Mean) < 0.5)] <- "Q3"
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0.5 & percent_rank(x = data2$PhyloP_Mean) < 0.75)] <- "Q2"
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0.75 & percent_rank(x = data2$PhyloP_Mean) <= 1)] <- "Q1"

data2$Species <- rep("Zea mays", nrow(data2))
data2$Comparison <- rep("Zea-Oryza", nrow(data2))
data2 <- data2[,c(5,6,1,3,4,2)]

## hiHMM Chromatin States (CS)

data2 <- data %>%
  group_by(Gene) %>%
  summarise(CS_Composition = paste(CS_Composition, collapse = ","))

data2$Species <- rep("Zea mays", nrow(data2))
data2$Comparison <- rep("Zea-Oryza", nrow(data2))
data2 <- data2[,c(3,4,1,2)]

## PhastCons

data2 <- data %>%
  group_by(Gene) %>%
  summarise(PhastCons_Element = paste(PhastCons_Element, collapse = ","))

data2$Species <- rep("Zea mays", nrow(data2))
data2$Comparison <- rep("Zea-Oryza", nrow(data2))
data2 <- data2[,c(3,4,1,2)]

## CNEs

data2 <- Zm_Os.CNEs %>%
  group_by(Gene) %>%
  summarise(CNEs_100kb = paste(CNEs_100kb, collapse = ","),
            CNEsNearest_Distance = unique(CNEsNearest_Distance),
            CNEsID = paste(CNEsID, collapse = ","))

data2$Species <- rep("Zea mays", nrow(data2))
data2$Comparison <- rep("Zea-Oryza", nrow(data2))
data2 <- data2[,c(5,6,1,2,3,4)]

## GWAS

Zm_Os.GWAS$PASTED <- paste0(Zm_Os.GWAS$V1, ":", Zm_Os.GWAS$V2, ":", Zm_Os.GWAS$V3)
Zm_Os.GWAS$GWAS.SNPs_Presence <- gsub("0", "No", fixed = T, x = Zm_Os.GWAS$GWAS.SNPs_Presence)
Zm_Os.GWAS$GWAS.SNPs_Presence <- gsub("1", "Yes", fixed = T, x = Zm_Os.GWAS$GWAS.SNPs_Presence)
colnames(Zm_Os.GWAS) <- c("V1", "V2", "V3", "Gene", "GWAS.SNPs_Presence", "GWAS.SNPs_Location", "PASTED")
data2 <- Zm_Os.GWAS %>%
  group_by(Gene) %>%
  summarise(GWAS.SNPs_Presence = paste(GWAS.SNPs_Presence, collapse = ","), GWAS.SNPs_Location = paste(GWAS.SNPs_Location, collapse = ","))

data2$GWAS.SNPs_Location <- gsub(".:-1:-1", "", data2$GWAS.SNPs_Location)
data2$GWAS.SNPs_Location <- gsub(",", "", data2$GWAS.SNPs_Location)
data2$Species <- rep("Zea mays", nrow(data2))
data2$Comparison <- rep("Zea-Oryza", nrow(data2))
data2 <- data2[,c(4,5,1,2,3)]

### Export & joint section ###

ZmOsLECIF.table <- data2
ZmOsLECIF.table.jic <- data2

ZmOsPhyloP.table <- data2
ZmOsPhyloP.table.jic <- data2

ZmOsCS.table <- data2
ZmOsCS.table.jic <- data2

ZmOsPhastCons.table <- data2
ZmOsPhastCons.table.jic <- data2

ZmOsCNEs.table <- data2
ZmOsCNEs.table.jic <- data2

ZmOsGWAS.table <- data2
ZmOsGWAS.table.jic <- data2

if(length(which((ZmOsLECIF.table$Gene == ZmOsPhyloP.table$Gene) == FALSE)) == 0){ print("Ready to join")}

ZmOs.LECIF.PhyloP.table <- as_tibble(cbind(ZmOsLECIF.table, ZmOsPhyloP.table[, 4:6]))

if(length(which((ZmOsCS.table$Gene == ZmOs.LECIF.PhyloP.table$Gene) == FALSE)) == 0){ print("Ready to join")}

ZmOs.LECIF.PhyloP.CS.table <- as_tibble(cbind(ZmOs.LECIF.PhyloP.table, ZmOsCS.table[,4]))

if(length(which((ZmOsPhastCons.table$Gene == ZmOs.LECIF.PhyloP.CS.table$Gene) == FALSE)) == 0){ print("Ready to join")}

ZmOs.LECIF.PhyloP.CS.PhastCons.table <- as_tibble(cbind(ZmOs.LECIF.PhyloP.CS.table, ZmOsPhastCons.table[,4]))

which((ZmOsCNEs.table$Gene %in% ZmOsLECIF.table$Gene) == FALSE)
#ZmOsCNEs.table <- ZmOsCNEs.table[-5445,]
if(length(which((ZmOsCNEs.table$Gene == ZmOs.LECIF.PhyloP.CS.PhastCons.table$Gene) == FALSE)) == 0){ print("Ready to join")}

ZmOs.LECIF.PhyloP.CS.PhastCons.CNEs.table <- as_tibble(cbind(ZmOs.LECIF.PhyloP.CS.PhastCons.table, ZmOsCNEs.table[,4:6]))

if(length(which((ZmOsGWAS.table$Gene == ZmOs.LECIF.PhyloP.CS.PhastCons.CNEs.table$Gene) == FALSE)) == 0){ print("Ready to join")}

ZmOs.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table <- as_tibble(cbind(ZmOs.LECIF.PhyloP.CS.PhastCons.CNEs.table, ZmOsGWAS.table[,4:5]))

## Custom table

data2 <- as_tibble(rbind(ZmAt.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table,ZmOs.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table,OsAt.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table,OsZm.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table, AtOs.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table, AtZm.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table))
data2 <- ZmOs.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table
sticky_style <- list(backgroundColor = "#e8e8e8")
html.object <- htmltools::browsable(
  tagList(
    tags$button(
      tagList(fontawesome::fa("download"), "Download Table"),
      onclick = "Reactable.downloadDataCSV('GenesTableFUNCO', 'GenesTableFUNCO.csv')"
    ),
    
    reactable(data2, searchable = T, filterable = T, pagination = TRUE, paginationType = "jump", defaultPageSize = 100, defaultColDef = colDef(vAlign = "center", headerVAlign = "center", align = "center", minWidth = 80, footer = function(values, name) { htmltools::div(name, style = list(fontWeight = 600))}), 
              bordered = T, striped = TRUE, highlight = TRUE, resizable = TRUE, wrap = FALSE, sortable = TRUE, elementId = "GenesTableFUNCO", defaultSorted = list(LECIF_Mean = "desc"), height = 800,
              columns = list(
                Species = colDef(sticky = "left",
                                 style = sticky_style,
                                 headerStyle = sticky_style),
                Comparison = colDef(sticky = "left",
                                    style = sticky_style,
                                    headerStyle = sticky_style),
                Gene = colDef(sticky = "left",
                              style = sticky_style,
                              headerStyle = sticky_style),
                LECIF_Quartile = colDef(
                  style = function(value){
                    if(value == "Q1"){
                      color <- "#309143"
                    } else if(value == "Q2"){
                      color <- "#8ace7e"
                    } else if(value == "Q3"){
                      color <- "#ff684c"
                    } else if(value == "Q4"){
                      color <- "#b60a1c"
                    }
                    list(color = color, fontWeight = "bold")
                  }
                ),
                LECIF_Trend = colDef(cell = function(value, index) {
                  sparkline(data2$LECIF_Trend[[index]])
                }),
                PhyloP_Quartile = colDef(
                  style = function(value){
                    if(value == "Q1"){
                      color <- "#309143"
                    } else if(value == "Q2"){
                      color <- "#8ace7e"
                    } else if(value == "Q3"){
                      color <- "#ff684c"
                    } else if(value == "Q4"){
                      color <- "#b60a1c"
                    }
                    list(color = color, fontWeight = "bold")
                  }
                ),
                PhyloP_Trend = colDef(cell = function(value, index) {
                  sparkline(data2$PhyloP_Trend[[index]])
                }),
                PhastCons_Element = colDef(cell = function(value) {
                  if (length(grep("Yes", value)) > 0) "\u2714\ufe0f Yes" else "\u274c No"
                }),
                CNEs_100kb = colDef(cell = function(value) {
                  if (length(grep("Yes", value)) > 0) "\u2714\ufe0f Yes" else "\u274c No"
                }),
                GWAS.SNPs_Presence = colDef(cell = function(value) {
                  if (length(grep("Yes", value)) > 0) "\u2714\ufe0f Yes" else "\u274c No"
                })
              ))
  )
)
save_html(html.object, "H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/Zea-Oryza/SearchGenes_Zm-Os.html")
save_html(html.object, "H:/Taiotactical/Taiotactical_LECIF/3.Applications/Tables/Tables/SearchGenes.html")

#### 2. At Super-enhancers Tables #### 

#### 2.1 Arabidopsis-Oryza ####

### Format ###

## LECIF

At_Os.LECIF <- read.delim("H:/Taiotactical/Taiotactical_LECIF/Tables_Superenhancers/At-Os.SuperEnhancersID.sorted.binned.LECIFscore.bed", header = F, sep = "\t")
At_Os.LECIF$V8[which(At_Os.LECIF$V8 == ".")] <- gsub(".", 0, fixed = T, At_Os.LECIF$V8[which(At_Os.LECIF$V8 == ".")])
At_Os.LECIF$V8 <- as.numeric(At_Os.LECIF$V8)
At_Os.LECIF <- At_Os.LECIF[,c(1,2,3,4,8)]

## PhyloP

At_Os.PhyloP <- read.delim("H:/Taiotactical/Taiotactical_LECIF/Tables_Superenhancers/SuperEnhancersID.sorted.binned.PhyloP.bed", header = F, sep = "\t")
At_Os.PhyloP$V8 <- as.numeric(At_Os.PhyloP$V8)
At_Os.PhyloP <- At_Os.PhyloP[,c(1,2,3,4,8)]
At_Os.PhyloP$V8[which(is.na(At_Os.PhyloP$V8) == TRUE)] <- 0

## hiHMM Chromatin States (CS)

At_Os.CS <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/Tables_Superenhancers/SuperEnhancersID.sorted.binned.CS.bed", header = F, sep = "\t")
At_Os.CS <- At_Os.CS[,c(1,2,3,4,8)]

## PhastCons

At_Os.PhastCons <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/Tables_Superenhancers/SuperEnhancersID.sorted.binned.PhastCons.bed", header = F, sep = "\t")
At_Os.PhastCons <- At_Os.PhastCons[,c(1,2,3,4,7)]

## CNEs 

At_Os.CNEs <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/Tables_Superenhancers/At-Os.SuperEnhancersID.sorted.CNEs.bed", header = F, sep = "\t")
At_Os.CNEs$CNEs_100kb <- "No"
At_Os.CNEs$CNEs_100kb[which(At_Os.CNEs$V14 <= 100000)] <- "Yes"
At_Os.CNEs$CNEsNearest_Distance <- At_Os.CNEs$V14
At_Os.CNEs$CNEsID <- NA
At_Os.CNEs$CNEsID <- At_Os.CNEs$V8
At_Os.CNEs$SuperEnhancerID <- At_Os.CNEs$V4
At_Os.CNEs <- At_Os.CNEs[,c(18,17,15,16)]

## SNPs-GWAS

At_Os.GWAS <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/Tables_Superenhancers/SuperEnhancersID.sorted.binned.SNPsGWAS.bed", header = F, sep = "\t")
At_Os.GWAS$GWAS.SNPs_Presence <- At_Os.GWAS$V8
At_Os.GWAS$GWAS.SNPs_Location <- paste0(At_Os.GWAS$V5,":",At_Os.GWAS$V6,":",At_Os.GWAS$V7)
At_Os.GWAS <- At_Os.GWAS[,c(1,2,3,4,9,10)]

## Genes

At_Os.Genes <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/Tables_Superenhancers/SuperEnhancers_NearestGenes.bed", header = F, sep = "\t")
At_Os.Genes <- At_Os.Genes[,c(4,8,9)]
colnames(At_Os.Genes) <- c("SuperEnhancerID", "GeneID", "Distance")

### First group by interval*SuperEnhancer ###

## Genes
## Already grouped by SuperEnhancer intervals

## LECIF

At_Os.LECIF$PASTED <- paste0(At_Os.LECIF$V1, ":", At_Os.LECIF$V2, ":", At_Os.LECIF$V3) 
data <- At_Os.LECIF %>%
  group_by(PASTED) %>%
  summarise(SuperEnhancerID = unique(V4),LECIFscore = list(V8))
data$LECIFscoreM <- unlist(lapply(data$LECIFscore, FUN = mean), use.names = F)

## PhyloP

At_Os.PhyloP$PASTED <- paste0(At_Os.PhyloP$V1, ":", At_Os.PhyloP$V2, ":", At_Os.PhyloP$V3)
data <- At_Os.PhyloP %>%
  group_by(PASTED) %>%
  summarise(SuperEnhancerID = unique(V4),PhyloPscore = list(V8))
data$PhyloPscoreM <- unlist(lapply(data$PhyloPscore, FUN = mean), use.names = F)

## hiHMM Chromatin States (CS)

At_Os.CS$PASTED <- paste0(At_Os.CS$V1, ":", At_Os.CS$V2, ":", At_Os.CS$V3)
data <- At_Os.CS %>%
  group_by(PASTED) %>%
  summarise(SuperEnhancerID = unique(V4),CS_Composition = paste(V8, collapse = ","))

## PhastCons 

At_Os.PhastCons$PASTED <- paste0(At_Os.PhastCons$V1, ":", At_Os.PhastCons$V2, ":", At_Os.PhastCons$V3)
At_Os.PhastCons$V7 <- gsub("phastCons_predicted", "Yes", fixed = T, x = At_Os.PhastCons$V7)
At_Os.PhastCons$V7 <- gsub(".", "No", fixed = T, x = At_Os.PhastCons$V7)
data <- At_Os.PhastCons %>%
  group_by(PASTED) %>%
  summarise(SuperEnhancerID = unique(V4), PhastCons_Element = paste(V7, collapse = ","))

## CNEs
## Already grouped by gene intervals

## GWAS

At_Os.GWAS$PASTED <- paste0(At_Os.GWAS$V1, ":", At_Os.GWAS$V2, ":", At_Os.GWAS$V3)
At_Os.GWAS$GWAS.SNPs_Presence <- gsub("0", "No", fixed = T, x = At_Os.GWAS$GWAS.SNPs_Presence)
At_Os.GWAS$GWAS.SNPs_Presence <- gsub("1", "Yes", fixed = T, x = At_Os.GWAS$GWAS.SNPs_Presence)
data <- At_Os.GWAS %>%
  group_by(PASTED) %>%
  summarise(SuperEnhancerID = unique(V4), GWAS.SNPs_Presence = paste(GWAS.SNPs_Presence, collapse = ","), GWAS.SNPs_Location = paste(GWAS.SNPs_Location, collapse = ","))

### Second group by gene only and mean five bins ###

## Genes

data2 <- At_Os.Genes %>%
  group_by(SuperEnhancerID) %>%
  summarise(GeneID = paste(GeneID, collapse = ","), GeneDistance = paste(Distance, collapse = ","))
data2$Species <- rep("Arabidopsis thaliana", nrow(data2))
data2$Comparison <- rep("Arabidopsis-Oryza", nrow(data2))
data2 <- data2[,c(5,6,1,3,4,2)]

## LECIF

data2 <- data %>%
  group_by(SuperEnhancerID) %>%
  summarise(LECIF_Trend = list(LECIFscoreM))

data2$LECIF_Mean <- unlist(lapply(data2$LECIF_Trend, mean), use.names = F)
data2$LECIF_Mean <- round(data2$LECIF_Mean, digits = 3)
data2$LECIF_Quartile <- NA
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0 & percent_rank(x = data2$LECIF_Mean) < 0.25)] <- "Q4"
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0.25 & percent_rank(x = data2$LECIF_Mean) < 0.5)] <- "Q3"
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0.5 & percent_rank(x = data2$LECIF_Mean) < 0.75)] <- "Q2"
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0.75 & percent_rank(x = data2$LECIF_Mean) <= 1)] <- "Q1"

data2$Species <- rep("Arabidopsis thaliana", nrow(data2))
data2$Comparison <- rep("Arabidopsis-Oryza", nrow(data2))
data2 <- data2[,c(4,5,1,2,3)]

## PhyloP

data2 <- data %>%
  group_by(Gene) %>%
  summarise(PhyloP_Trend = list(PhyloPscoreM))

data2$PhyloP_Mean <- unlist(lapply(data2$PhyloP_Trend, mean), use.names = F)
data2$PhyloP_Mean <- round(data2$PhyloP_Mean, digits = 3)
data2$PhyloP_Quartile <- NA
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0 & percent_rank(x = data2$PhyloP_Mean) < 0.25)] <- "Q4"
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0.25 & percent_rank(x = data2$PhyloP_Mean) < 0.5)] <- "Q3"
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0.5 & percent_rank(x = data2$PhyloP_Mean) < 0.75)] <- "Q2"
data2$PhyloP_Quartile[which(percent_rank(x = data2$PhyloP_Mean) >= 0.75 & percent_rank(x = data2$PhyloP_Mean) <= 1)] <- "Q1"

data2$Species <- rep("Arabidopsis thaliana", nrow(data2))
data2$Comparison <- rep("Arabidopsis-Oryza", nrow(data2))
data2 <- data2[,c(5,6,1,3,4,2)]

## hiHMM Chromatin States (CS)

data2 <- data %>%
  group_by(SuperEnhancerID) %>%
  summarise(CS_Composition = paste(CS_Composition, collapse = ","))

data2$Species <- rep("Arabidopsis thaliana", nrow(data2))
data2$Comparison <- rep("Arabidopsis-Oryza", nrow(data2))
data2 <- data2[,c(3,4,1,2)]

## PhastCons

data2 <- data %>%
  group_by(SuperEnhancerID) %>%
  summarise(PhastCons_Element = paste(PhastCons_Element, collapse = ","))

data2$Species <- rep("Arabidopsis thaliana", nrow(data2))
data2$Comparison <- rep("Arabidopsis-Oryza", nrow(data2))
data2 <- data2[,c(3,4,1,2)]

## CNEs

data2 <- At_Os.CNEs %>%
  group_by(SuperEnhancerID) %>%
  summarise(CNEs_100kb = paste(CNEs_100kb, collapse = ","),
            CNEsNearest_Distance = unique(CNEsNearest_Distance),
            CNEsID = paste(CNEsID, collapse = ","))

data2$Species <- rep("Arabidopsis thaliana", nrow(data2))
data2$Comparison <- rep("Arabidopsis-Oryza", nrow(data2))
data2 <- data2[,c(5,6,1,2,3,4)]

## GWAS

data2 <- data %>%
  group_by(SuperEnhancerID) %>%
  summarise(GWAS.SNPs_Presence = paste(GWAS.SNPs_Presence, collapse = ","), GWAS.SNPs_Location = paste(GWAS.SNPs_Location, collapse = ","))

data2$GWAS.SNPs_Location <- gsub(".:-1:-1", "", data2$GWAS.SNPs_Location)
data2$GWAS.SNPs_Location <- gsub(",", "", data2$GWAS.SNPs_Location)
data2$Species <- rep("Arabidopsis thaliana", nrow(data2))
data2$Comparison <- rep("Arabidopsis-Oryza", nrow(data2))
data2 <- data2[,c(4,5,1,2,3)]

### Export & joint section ###

AtOsGenes.table <- data2
AtOsGenes.table.jic <- data2

AtOsLECIF.table <- data2
AtOsLECIF.table.jic <- data2

AtOsPhyloP.table <- data2
AtOsPhyloP.table.jic <- data2

AtOsCS.table <- data2
AtOsCS.table.jic <- data2

AtOsPhastCons.table <- data2
AtOsPhastCons.table.jic <- data2

AtOsCNEs.table <- data2
AtOsCNEs.table.jic <- data2

AtOsGWAS.table <- data2
AtOsGWAS.table.jic <- data2

if(length(which((AtOsGenes.table$SuperEnhancerID == AtOsLECIF.table$SuperEnhancerID) == FALSE)) == 0){ print("Ready to join")}

AtOs.Genes.LECIF.table <- as_tibble(cbind(AtOsGenes.table, AtOsLECIF.table[, 4:6]))

if(length(which((AtOs.Genes.LECIF.table$SuperEnhancerID == AtOsPhyloP.table$SuperEnhancerID) == FALSE)) == 0){ print("Ready to join")}

AtOs.Genes.LECIF.PhyloP.table <- as_tibble(cbind(AtOs.Genes.LECIF.table, AtOsPhyloP.table[, 4:6]))

if(length(which((AtOs.Genes.LECIF.PhyloP.table$SuperEnhancerID == AtOsCS.table$SuperEnhancerID) == FALSE)) == 0){ print("Ready to join")}

AtOs.Genes.LECIF.PhyloP.CS.table <- as_tibble(cbind(AtOs.Genes.LECIF.PhyloP.table, AtOsCS.table[,4]))

if(length(which((AtOsPhastCons.table$SuperEnhancerID == AtOs.Genes.LECIF.PhyloP.CS.table$SuperEnhancerID) == FALSE)) == 0){ print("Ready to join")}

AtOs.Genes.LECIF.PhyloP.CS.PhastCons.table <- as_tibble(cbind(AtOs.Genes.LECIF.PhyloP.CS.table, AtOsPhastCons.table[,4]))

if(length(which((AtOsCNEs.table$SuperEnhancerID == AtOs.Genes.LECIF.PhyloP.CS.PhastCons.table$SuperEnhancerID) == FALSE)) == 0){ print("Ready to join")}

AtOs.Genes.LECIF.PhyloP.CS.PhastCons.CNEs.table <- as_tibble(cbind(AtOs.Genes.LECIF.PhyloP.CS.PhastCons.table, AtOsCNEs.table[,4:6]))

if(length(which((AtOsGWAS.table$SuperEnhancerID == AtOs.Genes.LECIF.PhyloP.CS.PhastCons.CNEs.table$SuperEnhancerID) == FALSE)) == 0){ print("Ready to join")}

AtOs.Genes.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table <- as_tibble(cbind(AtOs.Genes.LECIF.PhyloP.CS.PhastCons.CNEs.table, AtOsGWAS.table[,4:5]))

#### 2.2 Arabidopsis-Zea ####

### Format ###

## LECIF

At_Zm.LECIF <- read.delim("H:/Taiotactical/Taiotactical_LECIF/Tables_Superenhancers/At-Zm.SuperEnhancersID.sorted.binned.LECIFscore.bed", header = F, sep = "\t")
At_Zm.LECIF$V8[which(At_Zm.LECIF$V8 == ".")] <- gsub(".", 0, fixed = T, At_Zm.LECIF$V8[which(At_Zm.LECIF$V8 == ".")])
At_Zm.LECIF$V8 <- as.numeric(At_Zm.LECIF$V8)
At_Zm.LECIF <- At_Zm.LECIF[,c(1,2,3,4,8)]

## CNEs 

At_Zm.CNEs <-  read.delim("H:/Taiotactical/Taiotactical_LECIF/Tables_Superenhancers/At-Zm.SuperEnhancersID.sorted.CNEs.bed", header = F, sep = "\t")
At_Zm.CNEs$CNEs_100kb <- "No"
At_Zm.CNEs$CNEs_100kb[which(At_Zm.CNEs$V14 <= 100000)] <- "Yes"
At_Zm.CNEs$CNEsNearest_Distance <- At_Zm.CNEs$V14
At_Zm.CNEs$CNEsID <- NA
At_Zm.CNEs$CNEsID <- At_Zm.CNEs$V8
At_Zm.CNEs$SuperEnhancerID <- At_Zm.CNEs$V4
At_Zm.CNEs <- At_Zm.CNEs[,c(18,17,15,16)]

### First group by interval*gene and mean n LECIF hits ###

## LECIF

At_Zm.LECIF$PASTED <- paste0(At_Zm.LECIF$V1, ":", At_Zm.LECIF$V2, ":", At_Zm.LECIF$V3) 
data <- At_Zm.LECIF %>%
  group_by(PASTED) %>%
  summarise(SuperEnhancerID = unique(V4),LECIFscore = list(V8))
data$LECIFscoreM <- unlist(lapply(data$LECIFscore, FUN = mean), use.names = F)

## CNEs
## Already grouped by gene intervals

### Second group by gene only and mean five bins ###

## LECIF

data2 <- data %>%
  group_by(SuperEnhancerID) %>%
  summarise(LECIF_Trend = list(LECIFscoreM))

data2$LECIF_Mean <- unlist(lapply(data2$LECIF_Trend, mean), use.names = F)
data2$LECIF_Mean <- round(data2$LECIF_Mean, digits = 3)
data2$LECIF_Quartile <- NA
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0 & percent_rank(x = data2$LECIF_Mean) < 0.25)] <- "Q4"
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0.25 & percent_rank(x = data2$LECIF_Mean) < 0.5)] <- "Q3"
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0.5 & percent_rank(x = data2$LECIF_Mean) < 0.75)] <- "Q2"
data2$LECIF_Quartile[which(percent_rank(x = data2$LECIF_Mean) >= 0.75 & percent_rank(x = data2$LECIF_Mean) <= 1)] <- "Q1"

data2$Species <- rep("Arabidopsis thaliana", nrow(data2))
data2$Comparison <- rep("Arabidopsis-Zea", nrow(data2))
data2 <- data2[,c(5,6,1,3,4,2)]

## CNEs

data2 <- At_Zm.CNEs %>%
  group_by(SuperEnhancerID) %>%
  summarise(CNEs_100kb = paste(CNEs_100kb, collapse = ","),
            CNEsNearest_Distance = unique(CNEsNearest_Distance),
            CNEsID = paste(CNEsID, collapse = ","))

data2$Species <- rep("Arabidopsis thaliana", nrow(data2))
data2$Comparison <- rep("Arabidopsis-Zea", nrow(data2))
data2 <- data2[,c(5,6,1,2,3,4)]

### Export & joint section ###

AtZmLECIF.table <- data2
AtZmLECIF.table.jic <- data2

AtZmCNEs.table <- data2
AtZmCNEs.table.jic <- data2

if(length(which((AtOsGenes.table$SuperEnhancerID == AtZmLECIF.table$SuperEnhancerID) == FALSE)) == 0){ print("Ready to join")}

AtZm.Genes.LECIF.table <- as_tibble(cbind(AtOsGenes.table, AtZmLECIF.table[, 4:6]))

if(length(which((AtZm.Genes.LECIF.table$SuperEnhancerID == AtOsPhyloP.table$SuperEnhancerID) == FALSE)) == 0){ print("Ready to join")}

AtZm.Genes.LECIF.PhyloP.table <- as_tibble(cbind(AtZm.Genes.LECIF.table, AtOsPhyloP.table[, 4:6]))

if(length(which((AtZm.Genes.LECIF.PhyloP.table$SuperEnhancerID == AtOsCS.table$SuperEnhancerID) == FALSE)) == 0){ print("Ready to join")}

AtZm.Genes.LECIF.PhyloP.CS.table <- as_tibble(cbind(AtZm.Genes.LECIF.PhyloP.table, AtOsCS.table[,4]))

if(length(which((AtOsPhastCons.table$SuperEnhancerID == AtZm.Genes.LECIF.PhyloP.CS.table$SuperEnhancerID) == FALSE)) == 0){ print("Ready to join")}

AtZm.Genes.LECIF.PhyloP.CS.PhastCons.table <- as_tibble(cbind(AtZm.Genes.LECIF.PhyloP.CS.table, AtOsPhastCons.table[,4]))

if(length(which((AtZmCNEs.table$SuperEnhancerID == AtZm.Genes.LECIF.PhyloP.CS.PhastCons.table$SuperEnhancerID) == FALSE)) == 0){ print("Ready to join")}

AtZm.Genes.LECIF.PhyloP.CS.PhastCons.CNEs.table <- as_tibble(cbind(AtZm.Genes.LECIF.PhyloP.CS.PhastCons.table, AtZmCNEs.table[,4:6]))

if(length(which((AtOsGWAS.table$SuperEnhancerID == AtZm.Genes.LECIF.PhyloP.CS.PhastCons.CNEs.table$SuperEnhancerID) == FALSE)) == 0){ print("Ready to join")}

AtZm.Genes.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table <- as_tibble(cbind(AtZm.Genes.LECIF.PhyloP.CS.PhastCons.CNEs.table, AtOsGWAS.table[,4:5]))

## Custom table

AtZm.Genes.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table$Comparison <- gsub("Arabidopsis-Oryza", "Arabidopsis-Zea", AtZm.Genes.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table$Comparison, fixed = T)
data2 <- as_tibble(rbind(AtOs.Genes.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table, AtZm.Genes.LECIF.PhyloP.CS.PhastCons.CNEs.GWAS.table))
sticky_style <- list(backgroundColor = "#e8e8e8")
html.object <- htmltools::browsable(
  tagList(
    tags$button(
      tagList(fontawesome::fa("download"), "Download Table"),
      onclick = "Reactable.downloadDataCSV('SuperEnhancersTableFUNCO', 'SuperEnhancersTableFUNCO.csv')"
    ),
    
    reactable(data2, searchable = T, filterable = T, pagination = TRUE, paginationType = "jump", defaultPageSize = 100, defaultColDef = colDef(vAlign = "center", headerVAlign = "center", align = "center", minWidth = 80, footer = function(values, name) { htmltools::div(name, style = list(fontWeight = 600))}), 
              bordered = T, striped = TRUE, highlight = TRUE, resizable = TRUE, wrap = FALSE, sortable = TRUE, elementId = "SuperEnhancersTableFUNCO", defaultSorted = list(LECIF_Mean = "desc"), height = 800,
              columns = list(
                Species = colDef(sticky = "left",
                                 style = sticky_style,
                                 headerStyle = sticky_style),
                Comparison = colDef(sticky = "left",
                                    style = sticky_style,
                                    headerStyle = sticky_style),
                SuperEnhancerID = colDef(sticky = "left",
                              style = sticky_style,
                              headerStyle = sticky_style),
                LECIF_Quartile = colDef(
                  style = function(value){
                    if(value == "Q1"){
                      color <- "#309143"
                    } else if(value == "Q2"){
                      color <- "#8ace7e"
                    } else if(value == "Q3"){
                      color <- "#ff684c"
                    } else if(value == "Q4"){
                      color <- "#b60a1c"
                    }
                    list(color = color, fontWeight = "bold")
                  }
                ),
                LECIF_Trend = colDef(cell = function(value, index) {
                  sparkline(data2$LECIF_Trend[[index]])
                }),
                PhyloP_Quartile = colDef(
                  style = function(value){
                    if(value == "Q1"){
                      color <- "#309143"
                    } else if(value == "Q2"){
                      color <- "#8ace7e"
                    } else if(value == "Q3"){
                      color <- "#ff684c"
                    } else if(value == "Q4"){
                      color <- "#b60a1c"
                    }
                    list(color = color, fontWeight = "bold")
                  }
                ),
                PhyloP_Trend = colDef(cell = function(value, index) {
                  sparkline(data2$PhyloP_Trend[[index]])
                }),
                PhastCons_Element = colDef(cell = function(value) {
                  if (length(grep("Yes", value)) > 0) "\u2714\ufe0f Yes" else "\u274c No"
                }),
                CNEs_100kb = colDef(cell = function(value) {
                  if (length(grep("Yes", value)) > 0) "\u2714\ufe0f Yes" else "\u274c No"
                }),
                GWAS.SNPs_Presence = colDef(cell = function(value) {
                  if (length(grep("Yes", value)) > 0) "\u2714\ufe0f Yes" else "\u274c No"
                })
              ))
  )
)

save_html(html.object, "H:/Taiotactical/Taiotactical_LECIF/Tables_Superenhancers/SearchSuperEnhancers.html")


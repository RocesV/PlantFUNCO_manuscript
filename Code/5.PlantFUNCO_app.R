###################################
# TAIOTACTICAL- Overlap Enrichmet #
###################################

#### 0. Load libs and pkgs ####

library(BiocManager)
options(repos = BiocManager::repositories())
library(shiny)
library(shinyjs)
library(LOLA)
library(pheatmap)
library(RColorBrewer)
library(GenomicFeatures)
library(cowplot)
library(shinycssloaders)
library(rsconnect)

#### 1. Shiny app | UI ####

actionLink <- function(inputId, ...) {
  tags$a(href='javascript:void',
         id=inputId,
         class='action-button',
         ...)
}
options(spinner.color="#0275D8", spinner.color.background="#ffffff", spinner.size=2)
ui <- fluidPage(
  # Title Panel
  titlePanel("Genomic Overlap Enrichment"),
  actionButton("example", "Load Example"),
  actionButton("reset", "Clear"),
  actionButton("custom", "Load your sets and GO!"),
  downloadButton("export", "Export PDF"),
  hr(),
  # Main Panel
  fluidRow(
    column(width = 3, # Left Panel
           wellPanel( # Welcome Welcome section
             h4("Purpose"),
             textOutput("welcome1"), textOutput("welcome2"),
             textOutput("welcome3"), textOutput("welcome4"),
             textOutput("welcome5"), textOutput("welcome6"),
             textOutput("welcome7"), textOutput("welcome8"),
             textOutput("welcome9"),
             HTML("<br>"),
             HTML("<h5><strong>1) Select species, upload genomics regions and universe background</strong></h5>"),
             HTML("<h5><strong>2) Click GO! button to see results (it can take a while)</strong></h5>"),
             HTML("<h5><strong>3) Check results and export using PDF button </strong></h5>"),
             HTML("<h5><strong>Click Clear button to refresh or Load Example button to see results integrated with  Martin et al., 2021. Genome Biology (DOI: https://doi.org/10.1186/s13059-020-02258-y) data.</strong></h5>"),
             h5("This small application is focused in facilitating interpretation of functional genomics tracks generated (inter-species chromatin states and LECIF functional conservation score) and comparative genomics tracks used (PhyloP) in this study taking into account user-specified set of genomic regions.", style = "text-align: justify;"),
             h5("First of all, Users need to select one of the three species available (Arabidopsis, Oryza or Zea). Users can upload multiple Genomic Regions and 1 Universe Background corresponding to the species selected. Finally,  enrichment of the genomic regions of interest are tested for significant overlap against four databases (Chromatin States, LECIF vs Species1, LECIF vs Species2 and PhyloP) taking into account overlap distribution of the Universe background. Four heatmaps (1 heatmap x database) are plotted with OddsRatio results only highlighting significant overlaps.", style = "text-align: justify;"),
             HTML("<br>"),
             tags$small(paste0(
               "NOTE: LECIF and PhyloP are Score-based tracks, thus, regions were divided into equal sized bins (Bin5 = genomics regions with highest scores and vice versa for Bin1).",
               "Using the same set of genomic regions of interest, results could change depending of the Universe Background. See LOLA package for more info: http://code.databio.org/LOLA/"
             ))
           ),
           wellPanel( # Input section || NEED VALIDATE
             selectInput(inputId = "species", label = "Species", choices = c("Arabidopsis", "Oryza", "Zea")),
             fileInput("uploadgr", "Users Input BED files", buttonLabel = "Upload Genomic Regions", multiple = TRUE),
             fileInput("uploadub", NULL, buttonLabel = "Upload Universe Background", multiple = FALSE),
             tags$small(paste(
               "NOTE: Example data is only available for Arabidopsis thaliana and consists in differential expression and splicing regulation of abiotic, biotic and tissues genes/events characterized in this cool article: Martin et al., 2021. Genome Biology (DOI: https://doi.org/10.1186/s13059-020-02258-y). Universe Background are all genes of Arabidopsis thaliana genome.",
               "BED example for each Species:", 
               "cat Arabidopsis_thaliana.bed (only from Chr1 to Chr5)",
               "Chr1 200 450",
               "Chr2 250 730",
               "cat Oryza_sativa.bed (only from Chr1 to Chr12 and ChrSy/ChrUn)",
               "Chr1 200 450",
               "Chr2 250 730",
               "cat Zea_mays.file (only from 1 to 10)",
               "1 200 450",
               "2 250 730"
             ))
           )
    ),
    column(9, # Right panel and output section
           tableOutput("files1"),
           tableOutput("files2"),
           shinycssloaders::withSpinner(
           plotOutput("overlap"), type = 6)
             )),

)

#### 2. Shiny app | SERVER ####

server <- function(input, output, session) {
  #### general options ####
  options(shiny.maxRequestSize = 100 * 1024^2)
  #### Welcome Welcome section ####
  output$welcome1 <- renderText("--------------------------------------------------------")
  output$welcome2 <- renderText("Welcome, Welcome | Thaks for using <3")
  output$welcome3 <- renderText("--------------------------------------------------------")
  output$welcome4 <- renderText(".")
  output$welcome5 <- renderText(".")
  output$welcome6 <- renderText("( )_( )")
  output$welcome7 <- renderText("(='.'=)")
  output$welcome8 <- renderText("(^)_(^) [Roces, V.]  ")
  output$welcome9 <- renderText("")
  #### reactive expressions #### 
  v <- reactiveValues(uploadgr = NULL, uploadub = NULL)
  
  observeEvent(input$reset, {
    v$uploadgr <- NULL
    v$uploadub <- NULL
    v$OverlapPlot <- NULL
    reset("uploadgr")
    reset("uploadub")
    reset("OverlapPlot")
  })
  
  observeEvent(input$example, {
    v$species <- "Arabidopsis"
    #### Genomic Regions of Interest - EXAMPLE ####
    v$uploadgr <- data.frame(name = c("DiffGenes.abiotic", "DiffGenes.biotic", "DiffGenes.tissues", "DiffSplicing.abiotic", "DiffSplicing.biotic", "DiffSplicing.tissues"), 
                             datapath = c("./example/ArabidopsisDiffGenes.abiotic.bed", "./example/ArabidopsisDiffGenes.biotic.bed", "./example/ArabidopsisDiffGenes.tissues.bed",
                                          "./example/ArabidopsisDiffSplicing.abiotic.bed", "./example/ArabidopsisDiffSplicing.biotic.bed", "./example/ArabidopsisDiffSplicing.tissues.bed"))
    upload.genomicRegions <- list()
    for(i in 1:nrow(v$uploadgr)){
      upload.genomicRegions[[i]] <- read.table(file = v$uploadgr[i, 'datapath'], header = F, sep = "\t")
      upload.genomicRegions[[i]] <- upload.genomicRegions[[i]][,c(1,2,3)]
      colnames(upload.genomicRegions[[i]]) <- c("chr", "start", "end")
      upload.genomicRegions[[i]] <- makeGRangesFromDataFrame(upload.genomicRegions[[i]])
    }
    UserSets <- GRangesList(upload.genomicRegions)
    names(UserSets) <- v$uploadgr[,'name']
    rm(upload.genomicRegions)
    #### Universe Background - EXAMPLE ####
    v$uploadub <- data.frame(name = c("Universe"), datapath = c("./example/ArabidopsisGenome.background.bed"))
    upload.universeBackground <- read.table(file = v$uploadub[1, 'datapath'], header = F, sep = "\t")
    upload.universeBackground <- upload.universeBackground[,c(1,2,3)]
    colnames(upload.universeBackground) <- c("chr", "start", "end")
    Universe <- makeGRangesFromDataFrame(upload.universeBackground)
    rm(upload.universeBackground)
    
    #### Computation shit - EXAMPLE ####
    
    regionResults.CShiHMM <- lapply(UserSets, FUN = function(x){
      Results <- runLOLA(x, Universe, regionDB[[v$species]]$CShiHMM, cores = 1)
    })
    regionResults.LECIF1 <- lapply(UserSets, FUN = function(x){
      Results <- runLOLA(x, Universe, regionDB[[v$species]]$LECIF1, cores = 1)
    })
    regionResults.LECIF2 <- lapply(UserSets, FUN = function(x){
      Results <- runLOLA(x, Universe, regionDB[[v$species]]$LECIF2, cores = 1)
    })
    regionResults.PhyloP <- lapply(UserSets, FUN = function(x){
      Results <- runLOLA(x, Universe, regionDB[[v$species]]$PhyloP, cores = 1)
    })
    
    #### Aesthethic shit - EXAMPLE ####
    
    ### CShiHMM
    
    for(i in 1:length(regionResults.CShiHMM)){
      name <- names(regionResults.CShiHMM[i])
      regionResults.CShiHMM[[i]]$file <- rep(name, nrow(regionResults.CShiHMM[[i]]))
      regionResults.CShiHMM[[i]] <- regionResults.CShiHMM[[i]][,qValue:=NULL]
      regionResults.CShiHMM[[i]]$oddsRatio[which(regionResults.CShiHMM[[i]]$pValueLog < 1.30103)] <- 0
      regionResults.CShiHMM[[i]]$oddsRatio[which(regionResults.CShiHMM[[i]]$oddsRatio > 10)] <- 10
    }
    
    regionResults.CShiHMM <- do.call(rbind, regionResults.CShiHMM)
    df.CShiHMM <- data.frame(id = regionResults.CShiHMM$file, variable = regionResults.CShiHMM$description, value = regionResults.CShiHMM$oddsRatio)
    df.CShiHMM <- reshape(df.CShiHMM, idvar = "id", v.names = "value",  timevar = "variable", direction = "wide")
    colnames(df.CShiHMM) <- gsub("value.", "", colnames(df.CShiHMM), fixed = T)
    rownames(df.CShiHMM) <- df.CShiHMM$id
    df.CShiHMM <- df.CShiHMM[,-1]
    col.CShiHMM <- colorRampPalette(brewer.pal(7, "Greens"))(250)
    df.CShiHMM <- as.matrix(df.CShiHMM)
    title.CShiHMM <- "Inter-Specific Chromatin States"
    
    ### LECIF1
    
    for(i in 1:length(regionResults.LECIF1)){
      name <- names(regionResults.LECIF1[i])
      regionResults.LECIF1[[i]]$file <- rep(name, nrow(regionResults.LECIF1[[i]]))
      regionResults.LECIF1[[i]] <- regionResults.LECIF1[[i]][,qValue:=NULL]
      regionResults.LECIF1[[i]]$oddsRatio[which(regionResults.LECIF1[[i]]$pValueLog < 1.30103)] <- 0
      regionResults.LECIF1[[i]]$oddsRatio[which(regionResults.LECIF1[[i]]$oddsRatio > 10)] <- 10
    }
    regionResults.LECIF1 <- do.call(rbind, regionResults.LECIF1)
    df.LECIF1 <- data.frame(id = regionResults.LECIF1$file, variable = regionResults.LECIF1$description, value = regionResults.LECIF1$oddsRatio)
    df.LECIF1 <- reshape(df.LECIF1, idvar = "id", v.names = "value",  timevar = "variable", direction = "wide")
    colnames(df.LECIF1) <- gsub("value.", "", colnames(df.LECIF1), fixed = T)
    rownames(df.LECIF1) <- df.LECIF1$id
    df.LECIF1 <- df.LECIF1[,-1]
    col.LECIF1 <- colorRampPalette(brewer.pal(7, "Purples"))(250)
    df.LECIF1 <- as.matrix(df.LECIF1)
    title.LECIF1 <- "LECIF score vs Oryza"
    
    ### LECIF2
    
    for(i in 1:length(regionResults.LECIF2)){
      name <- names(regionResults.LECIF2[i])
      regionResults.LECIF2[[i]]$file <- rep(name, nrow(regionResults.LECIF2[[i]]))
      regionResults.LECIF2[[i]] <- regionResults.LECIF2[[i]][,qValue:=NULL]
      regionResults.LECIF2[[i]]$oddsRatio[which(regionResults.LECIF2[[i]]$pValueLog < 1.30103)] <- 0
      regionResults.LECIF2[[i]]$oddsRatio[which(regionResults.LECIF2[[i]]$oddsRatio > 10)] <- 10
    }
    regionResults.LECIF2 <- do.call(rbind, regionResults.LECIF2)
    df.LECIF2 <- data.frame(id = regionResults.LECIF2$file, variable = regionResults.LECIF2$description, value = regionResults.LECIF2$oddsRatio)
    df.LECIF2 <- reshape(df.LECIF2, idvar = "id", v.names = "value",  timevar = "variable", direction = "wide")
    colnames(df.LECIF2) <- gsub("value.", "", colnames(df.LECIF2), fixed = T)
    rownames(df.LECIF2) <- df.LECIF2$id
    df.LECIF2 <- df.LECIF2[,-1]
    col.LECIF2 <- colorRampPalette(brewer.pal(7, "Purples"))(250)
    df.LECIF2 <- as.matrix(df.LECIF2)
    title.LECIF2 <- "LECIF score vs Zea"
    
    ### PhyloP
    
    for(i in 1:length(regionResults.PhyloP)){
      name <- names(regionResults.PhyloP[i])
      regionResults.PhyloP[[i]]$file <- rep(name, nrow(regionResults.PhyloP[[i]]))
      regionResults.PhyloP[[i]] <- regionResults.PhyloP[[i]][,qValue:=NULL]
      regionResults.PhyloP[[i]]$oddsRatio[which(regionResults.PhyloP[[i]]$pValueLog < 1.30103)] <- 0
      regionResults.PhyloP[[i]]$oddsRatio[which(regionResults.PhyloP[[i]]$oddsRatio > 10)] <- 10
    }
    regionResults.PhyloP <- do.call(rbind, regionResults.PhyloP)
    df.PhyloP <- data.frame(id = regionResults.PhyloP$file, variable = regionResults.PhyloP$description, value = regionResults.PhyloP$oddsRatio)
    df.PhyloP <- reshape(df.PhyloP, idvar = "id", v.names = "value",  timevar = "variable", direction = "wide")
    colnames(df.PhyloP) <- gsub("value.", "", colnames(df.PhyloP), fixed = T)
    rownames(df.PhyloP) <- df.PhyloP$id
    df.PhyloP <- df.PhyloP[,-1]
    col.PhyloP <- colorRampPalette(brewer.pal(7, "Blues"))(250)
    df.PhyloP <- as.matrix(df.PhyloP)
    title.PhyloP <- "PhyloP score"
    
    #### Plots zone - EXAMPLE ####
    
    p.CShiHMM <- pheatmap(df.CShiHMM, scale = "none", color = col.CShiHMM,  cluster_rows = F, cluster_cols = F,  clustering_method = "ward.D2", cellwidth = 35, cellheight = 35,main = title.CShiHMM, border_color = "white", silent = T)
    p.LECIF1 <- pheatmap(df.LECIF1, scale = "none", color = col.LECIF1,  cluster_rows = F, cluster_cols = F,  clustering_method = "ward.D2", cellwidth = 35, cellheight = 35,main = title.LECIF1, border_color = "white", silent = T)
    p.LECIF2 <- pheatmap(df.LECIF2, scale = "none", color = col.LECIF2,  cluster_rows = F, cluster_cols = F,  clustering_method = "ward.D2", cellwidth = 35, cellheight = 35,main = title.LECIF2, border_color = "white", silent = T)
    p.PhyloP <- pheatmap(df.PhyloP, scale = "none", color = col.PhyloP,  cluster_rows = F, cluster_cols = F,  clustering_method = "ward.D2", cellwidth = 35, cellheight = 35,main = title.PhyloP, border_color = "white", silent = T)
    v$OverlapPlot <- ggdraw() +
      draw_plot(p.CShiHMM$gtable, x = 0, y = 0.61, width = 1, height = .4) +
      draw_plot(p.LECIF1$gtable, x = 0, y = 0, width = 0.35, height = .60) +
      draw_plot(p.LECIF2$gtable, x = 0.33, y = 0, width = 0.35, height = .60) +
      draw_plot(p.PhyloP$gtable, x = 0.66, y = 0, width = 0.35, height = .60)
    
  })
  
  observeEvent(input$custom, {
    v$species <- input$species
    #### Genomic Regions of Interest - CUSTOM ####
    v$uploadgr <- input$uploadgr
    upload.genomicRegions <- list()
    for(i in 1:nrow(input$uploadgr)){
      upload.genomicRegions[[i]] <- read.table(file = input$uploadgr[i, 'datapath'], header = F, sep = "\t")
      upload.genomicRegions[[i]] <- upload.genomicRegions[[i]][,c(1,2,3)]
      colnames(upload.genomicRegions[[i]]) <- c("chr", "start", "end")
      upload.genomicRegions[[i]] <- makeGRangesFromDataFrame(upload.genomicRegions[[i]])
    }
    UserSets <- GRangesList(upload.genomicRegions)
    names(UserSets) <- input$uploadgr[,'name']
    rm(upload.genomicRegions)
    #### Universe Background - CUSTOM ####
    v$uploadub <- input$uploadub
    upload.universeBackground <- read.table(file = input$uploadub[1, 'datapath'], header = F, sep = "\t")
    upload.universeBackground <- upload.universeBackground[,c(1,2,3)]
    colnames(upload.universeBackground) <- c("chr", "start", "end")
    Universe <- makeGRangesFromDataFrame(upload.universeBackground)
    rm(upload.universeBackground)
    
    #### Computation shit - CUSTOM ####
    
    regionResults.CShiHMM <- lapply(UserSets, FUN = function(x){
      Results <- runLOLA(x, Universe, regionDB[[v$species]]$CShiHMM, cores = 1)
    })
    regionResults.LECIF1 <- lapply(UserSets, FUN = function(x){
      Results <- runLOLA(x, Universe, regionDB[[v$species]]$LECIF1, cores = 1)
    })
    regionResults.LECIF2 <- lapply(UserSets, FUN = function(x){
      Results <- runLOLA(x, Universe, regionDB[[v$species]]$LECIF2, cores = 1)
    })
    regionResults.PhyloP <- lapply(UserSets, FUN = function(x){
      Results <- runLOLA(x, Universe, regionDB[[v$species]]$PhyloP, cores = 1)
    })
    
    #### Aesthethic shit - CUSTOM ####
    
    ### CShiHMM
    
    for(i in 1:length(regionResults.CShiHMM)){
      name <- names(regionResults.CShiHMM[i])
      regionResults.CShiHMM[[i]]$file <- rep(name, nrow(regionResults.CShiHMM[[i]]))
      regionResults.CShiHMM[[i]] <- regionResults.CShiHMM[[i]][,qValue:=NULL]
      regionResults.CShiHMM[[i]]$oddsRatio[which(regionResults.CShiHMM[[i]]$pValueLog < 1.30103)] <- 0
      regionResults.CShiHMM[[i]]$oddsRatio[which(regionResults.CShiHMM[[i]]$oddsRatio > 10)] <- 10
    }
    
    regionResults.CShiHMM <- do.call(rbind, regionResults.CShiHMM)
    df.CShiHMM <- data.frame(id = regionResults.CShiHMM$file, variable = regionResults.CShiHMM$description, value = regionResults.CShiHMM$oddsRatio)
    df.CShiHMM <- reshape(df.CShiHMM, idvar = "id", v.names = "value",  timevar = "variable", direction = "wide")
    colnames(df.CShiHMM) <- gsub("value.", "", colnames(df.CShiHMM), fixed = T)
    rownames(df.CShiHMM) <- df.CShiHMM$id
    df.CShiHMM <- df.CShiHMM[,-1]
    col.CShiHMM <- colorRampPalette(brewer.pal(7, "Greens"))(250)
    df.CShiHMM <- as.matrix(df.CShiHMM)
    title.CShiHMM <- "Inter-Specific Chromatin States"
    
    ### LECIF1
    
    for(i in 1:length(regionResults.LECIF1)){
      name <- names(regionResults.LECIF1[i])
      regionResults.LECIF1[[i]]$file <- rep(name, nrow(regionResults.LECIF1[[i]]))
      regionResults.LECIF1[[i]] <- regionResults.LECIF1[[i]][,qValue:=NULL]
      regionResults.LECIF1[[i]]$oddsRatio[which(regionResults.LECIF1[[i]]$pValueLog < 1.30103)] <- 0
      regionResults.LECIF1[[i]]$oddsRatio[which(regionResults.LECIF1[[i]]$oddsRatio > 10)] <- 10
    }
    regionResults.LECIF1 <- do.call(rbind, regionResults.LECIF1)
    df.LECIF1 <- data.frame(id = regionResults.LECIF1$file, variable = regionResults.LECIF1$description, value = regionResults.LECIF1$oddsRatio)
    df.LECIF1 <- reshape(df.LECIF1, idvar = "id", v.names = "value",  timevar = "variable", direction = "wide")
    colnames(df.LECIF1) <- gsub("value.", "", colnames(df.LECIF1), fixed = T)
    rownames(df.LECIF1) <- df.LECIF1$id
    df.LECIF1 <- df.LECIF1[,-1]
    col.LECIF1 <- colorRampPalette(brewer.pal(7, "Purples"))(250)
    df.LECIF1 <- as.matrix(df.LECIF1)
    if(v$species == "Arabidopsis"){
      title.LECIF1 <- "LECIF score vs Oryza"
    }else if(v$species == "Oryza" | v$species == "Zea"){
      title.LECIF1 <- "LECIF score vs Arabidopsis"
    }
    
    ### LECIF2
    
    for(i in 1:length(regionResults.LECIF2)){
      name <- names(regionResults.LECIF2[i])
      regionResults.LECIF2[[i]]$file <- rep(name, nrow(regionResults.LECIF2[[i]]))
      regionResults.LECIF2[[i]] <- regionResults.LECIF2[[i]][,qValue:=NULL]
      regionResults.LECIF2[[i]]$oddsRatio[which(regionResults.LECIF2[[i]]$pValueLog < 1.30103)] <- 0
      regionResults.LECIF2[[i]]$oddsRatio[which(regionResults.LECIF2[[i]]$oddsRatio > 10)] <- 10
    }
    regionResults.LECIF2 <- do.call(rbind, regionResults.LECIF2)
    df.LECIF2 <- data.frame(id = regionResults.LECIF2$file, variable = regionResults.LECIF2$description, value = regionResults.LECIF2$oddsRatio)
    df.LECIF2 <- reshape(df.LECIF2, idvar = "id", v.names = "value",  timevar = "variable", direction = "wide")
    colnames(df.LECIF2) <- gsub("value.", "", colnames(df.LECIF2), fixed = T)
    rownames(df.LECIF2) <- df.LECIF2$id
    df.LECIF2 <- df.LECIF2[,-1]
    col.LECIF2 <- colorRampPalette(brewer.pal(7, "Purples"))(250)
    df.LECIF2 <- as.matrix(df.LECIF2)
    if(v$species == "Arabidopsis" | v$species == "Oryza"){
      title.LECIF2 <- "LECIF score vs Zea"
    }else if(v$species == "Zea"){
      title.LECIF2 <- "LECIF score vs Oryza"
    }
    
    
    ### PhyloP
    
    for(i in 1:length(regionResults.PhyloP)){
      name <- names(regionResults.PhyloP[i])
      regionResults.PhyloP[[i]]$file <- rep(name, nrow(regionResults.PhyloP[[i]]))
      regionResults.PhyloP[[i]] <- regionResults.PhyloP[[i]][,qValue:=NULL]
      regionResults.PhyloP[[i]]$oddsRatio[which(regionResults.PhyloP[[i]]$pValueLog < 1.30103)] <- 0
      regionResults.PhyloP[[i]]$oddsRatio[which(regionResults.PhyloP[[i]]$oddsRatio > 10)] <- 10
    }
    regionResults.PhyloP <- do.call(rbind, regionResults.PhyloP)
    df.PhyloP <- data.frame(id = regionResults.PhyloP$file, variable = regionResults.PhyloP$description, value = regionResults.PhyloP$oddsRatio)
    df.PhyloP <- reshape(df.PhyloP, idvar = "id", v.names = "value",  timevar = "variable", direction = "wide")
    colnames(df.PhyloP) <- gsub("value.", "", colnames(df.PhyloP), fixed = T)
    rownames(df.PhyloP) <- df.PhyloP$id
    df.PhyloP <- df.PhyloP[,-1]
    col.PhyloP <- colorRampPalette(brewer.pal(7, "Blues"))(250)
    df.PhyloP <- as.matrix(df.PhyloP)
    title.PhyloP <- "PhyloP score"
    
    #### Plots zone - CUSTOM ####
    
    p.CShiHMM <- pheatmap(df.CShiHMM, scale = "none", color = col.CShiHMM,  cluster_rows = F, cluster_cols = F,  clustering_method = "ward.D2", cellwidth = 35, cellheight = 35,main = title.CShiHMM, border_color = "white", silent = T)
    p.LECIF1 <- pheatmap(df.LECIF1, scale = "none", color = col.LECIF1,  cluster_rows = F, cluster_cols = F,  clustering_method = "ward.D2", cellwidth = 35, cellheight = 35,main = title.LECIF1, border_color = "white", silent = T)
    p.LECIF2 <- pheatmap(df.LECIF2, scale = "none", color = col.LECIF2,  cluster_rows = F, cluster_cols = F,  clustering_method = "ward.D2", cellwidth = 35, cellheight = 35,main = title.LECIF2, border_color = "white", silent = T)
    p.PhyloP <- pheatmap(df.PhyloP, scale = "none", color = col.PhyloP,  cluster_rows = F, cluster_cols = F,  clustering_method = "ward.D2", cellwidth = 35, cellheight = 35,main = title.PhyloP, border_color = "white", silent = T)
    v$OverlapPlot <- ggdraw() +
      draw_plot(p.CShiHMM$gtable, x = 0, y = 0.61, width = 1, height = .4) +
      draw_plot(p.LECIF1$gtable, x = 0, y = 0, width = 0.35, height = .60) +
      draw_plot(p.LECIF2$gtable, x = 0.33, y = 0, width = 0.35, height = .60) +
      draw_plot(p.PhyloP$gtable, x = 0.66, y = 0, width = 0.35, height = .60)
    
  })
  
  #### Outputs section ####
  ## Genomic Regions Tables ##
  output$files1 <- renderTable({ 
    # Reset button
    if (is.null(v$uploadgr)) return()
    # Plot table
    v$uploadgr
  })
  # Universe Background Table
  output$files2 <- renderTable({
    # Reset button
    if (is.null(v$uploadub)){
      return()
    }else if(!is.null(v$uploadub)){
      # Plot table
      v$uploadub  
    }
  })
  ## Heatmaps and export ##
  output$overlap <- renderPlot({
    # Reset button
    if (is.null(v$OverlapPlot)){
      return()
    }else if(!is.null(v$OverlapPlot)){
      # Plot heatmap comp
      v$OverlapPlot  
    }
  }, height = 800)
  # Export heatmap
  output$export <- downloadHandler(
      filename = "PlantFUNCOOverlapEnrichment.pdf",
      content = function(file) {
        pdf(file, onefile = TRUE, width = 16.5, height = 11.7)
        print(v$OverlapPlot)
        dev.off()
      }
  )
  
}

shinyApp(ui, server)







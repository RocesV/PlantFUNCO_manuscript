###################################
# TAIOTACTICAL-CNEseek MS SCRIPT #
###################################

# If running from my docker before executing this include NonNestedFilter.pl in $PATH
# 0. Load libs

library(CNEr)
library(rtracklayer)

#workDir
assemblyDir <- "/home/Transfer/0.WGA/CNEs/query/"
#resDir
axtDir <- "/home/Transfer/0.WGA/CNEs/alignments/"
#storeDir
cnesDir <- "/home/Transfer/0.WGA/CNEs/CNEs/"

# for loop two by two CNEs seeker
Species <- list.files(assemblyDir)
Combinations <- expand.grid(Species, Species)
Combinations <- Combinations[-which(Combinations[,1] == Combinations[,2]),]
Combinations <- paste0(Combinations[,1], sep = "-", Combinations[,2])

for(i in Combinations){
  # global variables and inputs needed
  Dir1to2 <- i
  Dir2to1 <- paste0(strsplit(Dir1to2, split = "-", fixed = T)[[1]][2], sep = "-", strsplit(Dir1to2, split = "-", fixed = T)[[1]][1])
  One <- strsplit(Dir1to2, split = "-", fixed = T)[[1]][1]
  Two <- strsplit(Dir1to2, split = "-", fixed = T)[[1]][2]
  ExonsOne <- read.table(file.path(paste0(assemblyDir,One,"/Exons_",One, ".bed")), col.names = c("chr", "start", "end", "strand")) 
  ExonsTwo <- read.table(file.path(paste0(assemblyDir,Two,"/Exons_",Two, ".bed")), col.names = c("chr", "start", "end", "strand"))
  RepeatsOne <- read.table(file.path(paste0(assemblyDir,One,"/RepeatMasker_",One, ".bed")), col.names = c("chr", "start", "end"))  
  RepeatsTwo <- read.table(file.path(paste0(assemblyDir,Two,"/RepeatMasker_",Two, ".bed")), col.names = c("chr", "start", "end"))
  SizesOne <- file.path(paste0(assemblyDir,One, "/", One, ".fa.sizes"))
  SizesTwo <- file.path(paste0(assemblyDir,Two, "/", Two, ".fa.sizes"))
  AssemblyOne <- file.path(paste0(assemblyDir,One, "/", One, ".fa.2bit"))
  AssemblyTwo <- file.path(paste0(assemblyDir,Two, "/", Two, ".fa.2bit"))
  axtFilesDir1to2 <- file.path(paste0(axtDir,Dir1to2,"/",gsub(pattern = "-", replacement = ".", Dir1to2, fixed = T),".net.axt"))
  axtFilesDir2to1 <- file.path(paste0(axtDir,Dir2to1,"/",gsub(pattern = "-", replacement = ".", Dir2to1, fixed = T),".net.axt"))
  # filter information to limit search to CNEs: Exons and Repeats
  Filter_ExonsOne <- GenomicRanges::makeGRangesFromDataFrame(ExonsOne)
  Filter_ExonsTwo <- GenomicRanges::makeGRangesFromDataFrame(ExonsTwo)
  Filter_RepeatsOne <- GenomicRanges::makeGRangesFromDataFrame(RepeatsOne)
  Filter_RepeatsTwo <- GenomicRanges::makeGRangesFromDataFrame(RepeatsTwo)
  FilterOne <- unlist(GenomicRanges::GRangesList(Filter_ExonsOne, Filter_RepeatsOne))
  FilterTwo <- unlist(GenomicRanges::GRangesList(Filter_ExonsTwo, Filter_RepeatsTwo))
  # create CNE class: hit values modified for each species according to additional whole genome duplications in Ren Ren et al 2018, Molecular Plant
  cutoffOne <- 4L
  cutoffTwo <- 4L 
  if(One == "Zea_mays"){
    cutoffOne <- 8L
  } else if(Two == "Zea_mays"){
    cutoffTwo <- 8L
  }
  cneObject <- CNE(assembly1Fn = AssemblyOne, assembly2Fn = AssemblyTwo, axt12Fn = axtFilesDir1to2, axt21Fn = axtFilesDir2to1, cutoffs1 = cutoffOne, cutoffs2 = cutoffTwo)
  # CNEs identification
  identities <- c(45L)
  windows <- c(50L)
  cneList <- ceScan(x = cneObject, tFilter = FilterOne, qFilter = FilterTwo, window = windows, identity = identities)
  cneMergedList <- lapply(cneList, cneMerge)
  if(length(cneMergedList[['45_50']]@CNEMerged) > 0){ 
    cneFinalList <- lapply(cneMergedList, blatCNE) # remember change to bowtie2 aligner in cneObject 
  } else {
      cat(paste0("NO 50 bp CNEs DETECTED FOR ", One, " and ", Two))
    # recursive conditional stopping the loop if not more CNEs seeker is needed
      Combinations <- Combinations[-which(Combinations %in% c(Dir1to2, Dir2to1))]
      if(length(Combinations) == 0){
        cat("Getting out of CNEs-seeker loop because no more combinations between species are available ...")
        break
        }
    }
  # CNE storage and query (SQLite + BedGraph borwser files)
  cat(paste0("\n Creating SQLdb for ", One, " and ", Two, " CNEs \n"))
  saveCNEToSQLite(x = cneFinalList[[1]], dbName = file.path(paste0(cnesDir, Dir1to2, "_", Dir2to1, "_CNEdb.sqlite")))
  setwd(cnesDir)
  makeCNEDensity(cneFinalList[["45_50"]]@CNEFinal, threshold = "45_50", genomeFirst = One, genomeSecond = Two)
  # recursive conditional stopping the loop if not more CNEs seeker is needed
  Combinations <- Combinations[-which(Combinations %in% c(Dir1to2, Dir2to1))]
  if(length(Combinations) == 0){
    cat("Getting out of CNEs-seeker loop because no more combinations between species are available ...")
    break
  }
}

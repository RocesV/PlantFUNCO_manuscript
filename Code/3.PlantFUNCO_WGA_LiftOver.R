###################################
# TAIOTACTICAL WGA-Lift MS SCRIPT #
###################################

# If running from my docker before executing this include NonNestedFilter.pl in $PATH
# 0. Load libs

library(CNEr)

## Modified functions for LiftOvers

my.system <- function(cmd, echo=TRUE, intern=FALSE, ...){
  if (echo){
    message(cmd)
  }
  res <- system(cmd, intern=intern, ...)
  if (!intern){
    stopifnot(res == 0)
  }
  return(res)
}

lastz_mod <- function (assemblyTarget, assemblyQuery, outputDir = ".", chrsTarget = NULL, 
                       chrsQuery = NULL, distance = c("far", "medium", "near"), 
                       binary = "lastz", mc.cores = getOption("mc.cores", 2L), 
                       echoCommand = FALSE) 
{
  distance <- match.arg(distance)
  if (!all(file.exists(c(assemblyTarget, assemblyQuery)))) {
    stop(assemblyTarget, " and ", assemblyQuery, " must exist!")
  }
  if (tools::file_ext(assemblyTarget) != "2bit" || tools::file_ext(assemblyQuery) != 
      "2bit") {
    stop("The assembly must be in .2bit format!")
  }
  if (!dir.exists(outputDir)) {
    dir.create(outputDir, recursive = TRUE)
  }
  matrixFile <- tempfile(fileext = ".lastzMatrix")
  if (!isTRUE(echoCommand)) {
    on.exit(unlink(matrixFile))
  }
  write.table(scoringMatrix(distance), file = matrixFile, 
              quote = FALSE, sep = " ", row.names = FALSE, col.names = TRUE)
  lastzOptions <- list(near = paste0("C=0 E=150 H=0 K=4500 L=3000 M=254 O=600 T=2 Y=15000 Q=", 
                                     matrixFile), medium = paste0("C=0 E=30 H=0 K=3000 L=3000 M=50 O=400 T=1 Y=9400 Q=", 
                                                                  matrixFile), far = paste0("C=0 E=30 H=2000 K=2200 L=6000 M=50 O=400 T=2 Y=3400 Q=", 
                                                                                            matrixFile))
  if (is.null(chrsTarget)) {
    chrsTarget <- seqnames(seqinfo(rtracklayer::TwoBitFile(assemblyTarget)))
  }
  else {
    stopifnot(all(chrsTarget %in% seqnames(seqinfo(rtracklayer::TwoBitFile(assemblyTarget)))))
  }
  if (is.null(chrsQuery)) {
    chrsQuery <- seqnames(seqinfo(rtracklayer::TwoBitFile(assemblyQuery)))
  }
  else {
    stopifnot(all(chrsQuery %in% seqnames(seqinfo(rtracklayer::TwoBitFile(assemblyQuery)))))
  }
  outputToReturn <- c()
  format <- "lav"
  runLastz <- function(chrTarget, chrQuery, assemblyTarget, 
                        assemblyQuery, format = "lav") {
    output <- file.path(outputDir, paste0(gsub("|", "\\|", 
                                               chrTarget, fixed = TRUE), ".", sub("\\..*$", "", 
                                                                                  basename(assemblyTarget)), "-", gsub("|", "\\|", 
                                                                                                                       chrQuery, fixed = TRUE), ".", sub("\\..*$", "", 
                                                                                                                                                         basename(assemblyQuery)), ".", format))
    if (identical(paste(assemblyTarget, chrTarget), paste(assemblyQuery, 
                                                          chrQuery))) {
      cmd <- paste0(binary, " ", assemblyTarget, "/", 
                    gsub("|", "\\|", chrTarget, fixed = TRUE), " ", 
                    "--self --nomirror", " ", lastzOptions[[distance]], 
                    " --format=", format, " --output=", output, 
                    " --markend --allocate:traceback=1.99G")
    }
    else {
      cmd <- paste0(binary, " ", assemblyTarget, "/", 
                    gsub("|", "\\|", chrTarget, fixed = TRUE), " ", 
                    assemblyQuery, "/", gsub("|", "\\|", chrQuery, 
                                             fixed = TRUE), " ", lastzOptions[[distance]], 
                    " --format=", format, " --output=", output, 
                    " --markend --allocate:traceback=1.99G")
    }
    if (echoCommand) {
      output <- cmd
    }
    else {
      if (!file.exists(output)) {
        my.system(cmd)
      }
      else {
        warning("The output ", output, " already exists! Skipping..")
      }
    }
    return(output)
  }
  combinations <- expand.grid(chrsTarget, chrsQuery, stringsAsFactors = FALSE)
  mc.cores <- min(mc.cores, nrow(combinations))
  ans <- parallel::mcmapply(runLastz, combinations[[1]], combinations[[2]], 
                  MoreArgs = list(assemblyTarget = assemblyTarget, assemblyQuery = assemblyQuery, 
                                  format = format), mc.cores = mc.cores)
  invisible(unname(ans))
}


#workDir
assemblyQueryDir <- "/home/Transfer/Liftover/query/"
assemblyTargetDir <- "/home/Transfer/Liftover/target/"
#resDir
axtDir <- "/home/Transfer/Liftover/alignments/"


for(i in list.files(assemblyTargetDir)){
  for(j in list.files(assemblyQueryDir)){
      dir.create(paste0(axtDir, i, sep = "-", j))
      outputDir <- paste0(axtDir, i, sep = "-", j)
      assemblyTarget <- paste0(assemblyTargetDir, i, sep = "/", i, ".fa.2bit")
      assemblyQuery  <- paste0(assemblyQueryDir, j, sep = "/", j, ".fa.2bit")
      distance <- "near"
      # Determination if it is a fragmented assembly or not
      nlinesTarget <- as.numeric(strsplit(system(paste0("wc -l ",assemblyTargetDir, i, sep = "/", i, ".fa.sizes"), intern = T), split = " ")[[1]][1])
      nlinesQuery <- as.numeric(strsplit(system(paste0("wc -l ", assemblyQueryDir, j, sep = "/", j, ".fa.sizes"), intern = T), split = " ")[[1]][1])
      if(nlinesQuery < 150 & nlinesTarget < 150){
        cat("Not fragmented genomes detected. Using lastz alignment ...")
        lavs <- lastz_mod(assemblyTarget = assemblyTarget, assemblyQuery = assemblyQuery, 
                      outputDir =     outputDir, chrsTarget = NULL, chrsQuery = NULL, distance = distance, 
                      mc.cores = 6, binary = "/home/rocesv/lastz_install/lastz_32") #aligner
        #lav files to psl files
        psls <- lavToPsl(lavs, removeLav = FALSE, 
                         binary = "/home/rocesv/HillerGenomeUtils/GenomeAlignmentTools/kent/bin/lavToPsl")
      }else if(nlinesQuery > 150 | nlinesTarget > 150){
        cat("Fragmented genomes detected. Using last alignment ...")
        ## Build the lastdb index
        system2(command="lastdb", args=c("-c", file.path(paste0(outputDir,sep = "/", i)),
                                         file.path(paste0(assemblyTargetDir, i, sep = "/", i, ".fa"))))
        ## Run last aligner
        lastal(db=file.path(paste0(outputDir, sep = "/", i)),
               queryFn=file.path(paste0(assemblyQueryDir, j, sep = "/", j, ".fa")),
               outputFn=file.path(paste0(outputDir, i, sep = "-", j, ".maf")),
               distance=distance, binary="lastal", mc.cores=6L)
        
        ## maf to psl 
        psls <- file.path(paste0(outputDir, i, sep = "-", j, ".psl"))
        system2(command="maf-convert", args=c("psl", 
                                              file.path(paste0(outputDir, i, sep = "-", j, ".maf")),
                                              ">", psls))
      }
      
      
      #### Chaining and processing - If two matching alignments are close enough, they are joined into one fragment in a process called chaining (axtChain). Then these chain files are sorted and combined into one file (chainMergeSort).
      
      cat("Chaining and processing ...")
      chains <- axtChain(psls = psls, assemblyTarget = assemblyTarget, assemblyQuery = assemblyQuery,
                         distance = distance, removePsl = FALSE, binary = "/home/rocesv/HillerGenomeUtils/GenomeAlignmentTools/kent/bin/axtChain")
      allChain <- chainMergeSort(chains = chains, assemblyTarget = assemblyTarget, assemblyQuery = assemblyQuery, removeChains = FALSE, 
                                 binary = "/home/rocesv/HillerGenomeUtils/GenomeAlignmentTools/kent/bin/chainMergeSort", 
                                 allChain = file.path(outputDir, paste0(sub("\\.fa.2bit$", "", basename(assemblyTarget),
                                                                            ignore.case=TRUE), ".", sub("\\.fa.2bit$", "", basename(assemblyQuery), ignore.case=TRUE), ".all.chain")))
      
      
      
      
      #### Netting - Filter chains to remove chains that do not have a chance of being netted by *chainPreNet*. During the alignment, every genomic fragment can match with several others, and certainly we want to keep the best one by *chainNet*. Then we add the synteny information with *netSyntenic*.
      
      cat("Netting and filtering chains ...")
      allPreChain <- chainPreNet(file.path(outputDir, paste0(sub("\\.fa.2bit$", "", basename(assemblyTarget),
                                                                 ignore.case=TRUE), ".", sub("\\.fa.2bit$", "", basename(assemblyQuery), ignore.case=TRUE), ".all.chain")), 
                                 assemblyTarget, assemblyQuery, removeAllChain = FALSE, 
                                 binary ="/home/rocesv/HillerGenomeUtils/GenomeAlignmentTools/kent/bin/chainPreNet",
                                 allPreChain=file.path(outputDir, paste0(sub("\\.fa.2bit$", "", basename(assemblyTarget), ignore.case = TRUE), ".", 
                                                                         sub("\\.fa.2bit$", "", basename(assemblyQuery), ignore.case = TRUE),".all.pre.chain")))     
      ### Bash Docker
      tNibDir <- paste0("-tNibDir=",assemblyTarget)
      qNibDir <- paste0("-qNibDir=",assemblyQuery)
      system2(command="/home/rocesv/HillerGenomeUtils/GenomeAlignmentTools/bin/chainNet", 
              args=c(file.path(outputDir, paste0(sub("\\.fa.2bit$", "", basename(assemblyTarget),
                                                     ignore.case=TRUE), ".", sub("\\.fa.2bit$", "", basename(assemblyQuery), ignore.case=TRUE), ".all.pre.chain")),
                     file.path(paste0(assemblyTargetDir, i, sep = "/", i, ".fa.sizes")),
                     file.path(paste0(assemblyQueryDir, j, sep = "/", j, ".fa.sizes")),
                     file.path(paste0(outputDir,sep="/", i, ".net")),
                     file.path(paste0(outputDir,sep="/", j, ".net")),
                     "-rescore",
                     tNibDir,
                     qNibDir,
                     "-linearGap=medium"
              ))
      
      system2(command="netSyntenic",
              args = c(file.path(paste0(outputDir,sep="/", i, ".net")),
                       file.path(paste0(outputDir,sep="/", i, ".noClass.net"))
              ))
      
      #### axtNet - Create *.net.axt* file from the previous *net* and *chain* files.
      cat("Passing to .net.axt format ...")
      netSyntenicFile <- file.path(paste0(outputDir,sep="/", i, ".noClass.net"))
      netToAxt(netSyntenicFile, allPreChain, assemblyTarget, assemblyQuery, axtFile=file.path(outputDir, paste0(sub("\\.fa.2bit$", "", basename(assemblyTarget),
                                                                                                                    ignore.case = TRUE), ".", sub("\\.fa.2bit$", "", basename(assemblyQuery), ignore.case = TRUE), ".net.axt")), removeFiles=FALSE,
               binaryNetToAxt="/home/rocesv/HillerGenomeUtils/GenomeAlignmentTools/kent/bin/netToAxt", 
               binaryAxtSort="/home/rocesv/HillerGenomeUtils/GenomeAlignmentTools/kent/bin/axtSort")     
    }}


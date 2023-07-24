###################################
# TAIOTACTICAL WGA-CNEs MS SCRIPT #
###################################

# If running from my docker before executing this include NonNestedFilter.pl in $PATH
# 0. Load libs

library(CNEr)

#workDir
assemblyDir <- "/home/Transfer/CNEs/query/"
#resDir
axtDir <- "/home/Transfer/CNEs/alignments/"

#for loop - bidirectional lastz alignment between all species in query
for(i in list.files(assemblyDir)){
  for(j in list.files(assemblyDir)){
    if(i == j){print("Pass to next iteration")
    }else if(i != j){
      dir.create(paste0(axtDir, i, sep = "-", j))
      outputDir <- paste0(axtDir, i, sep = "-", j)
      assemblyTarget <- paste0(assemblyDir, i, sep = "/", i, ".fa.2bit")
      assemblyQuery  <- paste0(assemblyDir, j, sep = "/", j, ".fa.2bit")
      if(i == "Oryza_sativa" & j == "Zea_mays" | i == "Zea_mays" & j == "Oryza_sativa"){
        distance <- "medium"
      }else {
        distance <- "far"
      }
      if(i == "Arabidopsis_thaliana" & j == "Oryza_sativa"){
        next
      }else if(i == "Oryza_sativa" & j == "Arabidopsis_thaliana"){
        next
      }else if(i == "Oryza_sativa" & j == "Zea_mays"){
        next
      }else if(i == "Arabidopsis_thaliana" & j == "Zea_mays"){
        next
      }
      # Determination if it is a fragmented assembly or not
      nlinesTarget <- as.numeric(strsplit(system(paste0("wc -l ",assemblyDir, i, sep = "/", i, ".fa.sizes"), intern = T), split = " ")[[1]][1])
      nlinesQuery <- as.numeric(strsplit(system(paste0("wc -l ", assemblyDir, j, sep = "/", j, ".fa.sizes"), intern = T), split = " ")[[1]][1])
      if(nlinesQuery < 150 & nlinesTarget < 150){
        cat("Not fragmented genomes detected. Using lastz alignment ...")
        lavs <- lastz(assemblyTarget = assemblyTarget, assemblyQuery = assemblyQuery, 
                      outputDir =     outputDir, chrsTarget = NULL, chrsQuery = NULL, distance = distance, 
                      mc.cores = 6, binary = "/home/rocesv/lastz_install/lastz_32") #aligner
        #lav files to psl files
        psls <- lavToPsl(lavs, removeLav = FALSE, 
                         binary = "/home/rocesv/HillerGenomeUtils/GenomeAlignmentTools/kent/bin/lavToPsl")
      }else if(nlinesQuery > 150 | nlinesTarget > 150){
        cat("Fragmented genomes detected. Using last alignment ...")
        ## Build the lastdb index
        system2(command="lastdb", args=c("-c", file.path(paste0(outputDir,sep = "/", i)),
                                         file.path(paste0(assemblyDir, i, sep = "/", i, ".fa"))))
        ## Run last aligner
        lastal(db=file.path(paste0(outputDir, sep = "/", i)),
               queryFn=file.path(paste0(assemblyDir, j, sep = "/", j, ".fa")),
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
      
      
      #### Improve chains by Hiller group -- two rounds of senstivie alignment in order to improve non coding alignments | RepeatFiller & chainCleaner
      ### Bash Docker: careful with / in path with OutputDir
      
      cat("RepeatFiller improvement of chains ...")
      system2(command="/home/rocesv/HillerGenomeUtils/GenomeAlignmentTools/src/RepeatFiller.py", 
             args=c("-c", file.path(outputDir, paste0(sub("\\.fa.2bit$", "", basename(assemblyTarget),
                                                          ignore.case=TRUE), ".", sub("\\.fa.2bit$", "", basename(assemblyQuery), ignore.case=TRUE), ".all.chain")),
                    "-T2", assemblyTarget,
                    "-Q2", assemblyQuery,
                    "-o", file.path(outputDir, paste0(sub("\\.fa.2bit$", "", basename(assemblyTarget),
                                                          ignore.case=TRUE), ".", sub("\\.fa.2bit$", "", basename(assemblyQuery), ignore.case=TRUE), ".repeatfiller.all.chain"))
                    ))
      
      cat("ChainNet before ChainCleaner improvement of chains ...")
      system2(command = "chainNet",
              args = c("-minScore=0",
                       file.path(outputDir, paste0(sub("\\.fa.2bit$", "", basename(assemblyTarget),
                                                       ignore.case=TRUE), ".", sub("\\.fa.2bit$", "", basename(assemblyQuery), ignore.case=TRUE), ".repeatfiller.all.chain")),
                       file.path(paste0(assemblyDir, i, sep = "/", i, ".fa.sizes")),
                       file.path(paste0(assemblyDir, j, sep = "/", j, ".fa.sizes")),
                       "stdout /dev/null | NetFilterNonNested.perl /dev/stdin -minScore1 3000 >",
                       file.path(paste0(outputDir, "/tmp.chainCleaner.XXFXsAd.net"))
              ))
      
      
      cat("ChainCleaner improvement of chains ...")
      in.net <- paste0("-net=", file.path(paste0(outputDir, "/tmp.chainCleaner.XXFXsAd.net")))
      system2(command="chainCleaner", 
              args=c(file.path(outputDir, paste0(sub("\\.fa.2bit$", "", basename(assemblyTarget),
                                                     ignore.case=TRUE), ".", sub("\\.fa.2bit$", "", basename(assemblyQuery), ignore.case=TRUE), ".repeatfiller.all.chain")),
                     assemblyTarget,
                     assemblyQuery,
                     file.path(outputDir, paste0(sub("\\.fa.2bit$", "", basename(assemblyTarget),
                                                     ignore.case=TRUE), ".", sub("\\.fa.2bit$", "", basename(assemblyQuery), ignore.case=TRUE), ".chaincleaner.all.chain")),
                     file.path(outputDir, paste0("/chainCleaner.bed")),
                     in.net,
                     "-linearGap=medium"
                     ))
      
      #### Netting - Filter chains to remove chains that do not have a chance of being netted by *chainPreNet*. During the alignment, every genomic fragment can match with several others, and certainly we want to keep the best one by *chainNet*. Then we add the synteny information with *netSyntenic*.
      
      cat("Netting and filtering chains ...")
      allPreChain <- chainPreNet(file.path(outputDir, paste0(sub("\\.fa.2bit$", "", basename(assemblyTarget),
                                                                 ignore.case=TRUE), ".", sub("\\.fa.2bit$", "", basename(assemblyQuery), ignore.case=TRUE), ".chaincleaner.all.chain")), 
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
                     file.path(paste0(assemblyDir, i, sep = "/", i, ".fa.sizes")),
                     file.path(paste0(assemblyDir, j, sep = "/", j, ".fa.sizes")),
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
}


## data folder D:/Dropbox/Dropbox/FBNNet/FBNNet2/output/Experiments_data/leukeamia/differential_cutOffInduction_1_majority_dot_6
## files2<-'D:/Dropbox/Dropbox/FBNNet/ChildhoodLeukeamiaDataFile/GSE2677_RAW'
## files3<-'D:\\Dropbox\\Dropbox\\FBNNet\\FBNNet2\\study\\leukeamia_data'
## files1<-c('D:\\Dropbox\\Dropbox\\FBNNet\\ChildhoodLeukeamiaDataFile\\GSE2677_RAW','G:\\Dropbox\\Dropbox\\FBNNet\\Genome_Data\\GSE13670_RAW\\GSE13670_RAW','G:\\Dropbox\\Dropbox\\FBNNet\\Genome_Data\\GSE20489_RAW\\GSE20489_RAW','G:\\Dropbox\\Dropbox\\FBNNet\\Genome_Data\\GSE42088_RAW\\GSE42088_RAW'，'G:\\Dropbox\\Dropbox\\FBNNet\\Genome_Data\\GSE54992_RAW\\GSE54992_RAW'，'G:\\Dropbox\\Dropbox\\FBNNet\\Genome
## Data\\GSE57194_RAW\\GSE57194_RAW')

## targetsamples <- c("B-ALL-13", "B-ALL-17", "B-ALL-24", "B-ALL-31", "B-ALL-32", "B-ALL-33", "B-ALL-37", "B-ALL-38", "B-ALL-40", "B-ALL-43", "T-ALL-2", "T-ALL-20",
##                    "T-ALL-25")
## targetsamples <- c("B-ALL-13", "B-ALL-17", "B-ALL-24", "B-ALL-31", "B-ALL-32", "B-ALL-33", "B-ALL-37", "B-ALL-38", "B-ALL-40", "B-ALL-43")
## targetsamples <- c("T-ALL-2", "T-ALL-20","T-ALL-25")
## targetsamples <- c("B-ALL-Adult", "B-ALL-IV_EtOH_40", "B-ALL-IV_GC_40", "C-Line-C7R1dim_high_GC",
##                    "C-Line-CEMC1_ratGR_GC", "HD1-STS-1", "HD2-RPK-1", "R-Line-C7R1dim_low_GC",
##                    "R-Line-CEMC1_GC", "R-Line-PreB_EtOH", "R-Line-PreB_GC", "S-Line-C7H2_GC",        
##                    "S-Line-PreB_GC" )

#'@export
differentiallyExpressionStudy <- function(cellDirectory, 
                                  sortedtimeseries = NULL, 
                                  useGCRMA = FALSE, 
                                  cutOffInduction = 0.7, 
                                  cutOffRepression = 0.7, 
                                  majority = 7,   
                                  targetsamples = c("B-ALL-13", 
                                                     "B-ALL-17", 
                                                     "B-ALL-24", 
                                                     "B-ALL-31", 
                                                     "B-ALL-32", 
                                                     "B-ALL-33", 
                                                     "B-ALL-37", 
                                                     "B-ALL-38", 
                                                     "B-ALL-40", 
                                                     "B-ALL-43", 
                                                     "T-ALL-2", 
                                                     "T-ALL-20",
                                                     "T-ALL-25")) {
  
  # read affy files and normalized by RMA
  if (is.null(sortedtimeseries)) {
    sortedtimeseries2 <- convertAffyRawDataIntoNormalizedStructureData(cellDirectory, useGCRMA = useGCRMA)
  }
  sortedtimeseries <- sortedtimeseries2[targetsamples]
  
  cond <- sapply(sortedtimeseries, function(entry) !is.null(entry))
  sortedtimeseries <- sortedtimeseries[cond]
  

  print("###############################################our differentially experiments ###########################3")
  diffgenes_RMA <- identifyDifferentiallyExpressedGenes(sortedtimeseries, cutOffInduction = cutOffInduction, cutOffRepression = cutOffRepression, 
                                                        majority = majority)
  
  commonGeneSet1i <- diffgenes_RMA$DifferentialExpression[[1]]$Induced_ProbeID ##6to0
  commonGeneSet1r <- diffgenes_RMA$DifferentialExpression[[1]]$Repressed__ProbeID##6to0
  commonGeneSet2i <- diffgenes_RMA$DifferentialExpression[[2]]$Induced_ProbeID ##24to0
  commonGeneSet2r <- diffgenes_RMA$DifferentialExpression[[2]]$Repressed__ProbeID##24to0
  commonGeneSet3i <- diffgenes_RMA$DifferentialExpression[[3]]$Induced_ProbeID##24to6
  commonGeneSet3r <- diffgenes_RMA$DifferentialExpression[[3]]$Repressed__ProbeID##24to6
  
  print("############################")
  print("6To0 induced via FBNNET bioinformatics:")
  Genes6To0induced_fbnnet <- unlist(mapProbesetNames(commonGeneSet1i))
  deduplicated6To0induced_fbnnet <- unique(Genes6To0induced_fbnnet)
  print(deduplicated6To0induced_fbnnet)
  print(paste("In total =", length(commonGeneSet1i), " unqiue =", length(deduplicated6To0induced_fbnnet), sep = ""))
  
  print("############################")
  
  print("6To0 repressed via FBNNET bioinformatics:")
  Genes6To0repressed_fbnnet <- unlist(mapProbesetNames(commonGeneSet1r))
  deduplicated6To0repressed_fbnnet <- unique(Genes6To0repressed_fbnnet)
  print(deduplicated6To0repressed_fbnnet)
  print(paste("In total =", length(commonGeneSet1r), " unqiue =", length(deduplicated6To0repressed_fbnnet), sep = ""))
  
  
  print("############################")
  print("24To0 induced via FBNNET bioinformatics:")
  Genes24To0induced_fbnnet <- unlist(mapProbesetNames(commonGeneSet2i))
  deduplicated24To0induced_fbnnet <- unique(Genes24To0induced_fbnnet)
  print(deduplicated24To0induced_fbnnet)
  print(paste("In total =", length(commonGeneSet2i), " unqiue =", length(deduplicated24To0induced_fbnnet), sep = ""))
  
  print("############################")
  
  print("24To0 repressed via FBNNET bioinformatics:")
  Genes24To0repressed_fbnnet <- unlist(mapProbesetNames(commonGeneSet2r))
  deduplicated24To0repressed_fbnnet <- unique(Genes24To0repressed_fbnnet)
  print(deduplicated24To0repressed_fbnnet)
  print(paste("In total =", length(commonGeneSet2r), " unqiue =", length(deduplicated24To0repressed_fbnnet), sep = ""))
  
  
  print("############################")
  print("24To6 induced via FBNNET bioinformatics:")
  Genes24To6induced_fbnnet <- unlist(mapProbesetNames(commonGeneSet3i))
  deduplicated24To6induced_fbnnet <- unique(Genes24To6induced_fbnnet)
  print(deduplicated24To6induced_fbnnet)
  print(paste("In total =", length(commonGeneSet3i), " unqiue =", length(deduplicated24To6induced_fbnnet), sep = ""))
  
  print("############################")
  
  print("24To6 repressed via FBNNET bioinformatics:")
  Genes24To6repressed_fbnnet <- unlist(mapProbesetNames(commonGeneSet3r))
  deduplicated24To6repressed_fbnnet <- unique(Genes24To6repressed_fbnnet)
  print(deduplicated24To6repressed_fbnnet)
  print(paste("In total =", length(commonGeneSet3r), " unqiue =", length(deduplicated24To6repressed_fbnnet), sep = ""))
  
  ############################### output ################################################################################ todo: add the difference as well, plus the primise genes for
  ############################### further mining
  finalSet <- unique(c(deduplicated6To0induced_fbnnet, 
                       deduplicated6To0repressed_fbnnet, 
                       deduplicated24To0induced_fbnnet, 
                       deduplicated24To0repressed_fbnnet, 
                       deduplicated24To6induced_fbnnet, 
                       deduplicated24To6repressed_fbnnet))
  print(paste("Total unique genes cross all time points are:", length(finalSet), sep = ""))
  print(finalSet)
  finalProbesetSet <- unique(c(commonGeneSet1i, 
                               commonGeneSet1r, 
                               commonGeneSet2i, 
                               commonGeneSet2r, 
                               commonGeneSet3i, 
                               commonGeneSet3r))
  
  filtered_timeseries <- lapply(sortedtimeseries2, function(subdata) subdata[rownames(subdata) %in% finalProbesetSet, ])
  
  reproduceSchmidtStudy_res <- list()
  reproduceSchmidtStudy_res[["CombinedGeneSet"]] <- finalSet
  reproduceSchmidtStudy_res[["FBNNetGenes6To0induced"]] <- deduplicated6To0induced_fbnnet
  reproduceSchmidtStudy_res[["FBNNetGenes6To0repressed"]] <- deduplicated6To0repressed_fbnnet
  reproduceSchmidtStudy_res[["FBNNetGenes24To0induced"]] <- deduplicated24To0induced_fbnnet
  reproduceSchmidtStudy_res[["FBNNetGenes24To0repressed"]] <- deduplicated24To0repressed_fbnnet
  reproduceSchmidtStudy_res[["FBNNetGenes24To6Or8induced"]] <- deduplicated24To6induced_fbnnet
  reproduceSchmidtStudy_res[["FBNNetGenes24To6or8repressed"]] <- deduplicated24To6repressed_fbnnet
  reproduceSchmidtStudy_res[["timeseries_original"]] <- sortedtimeseries2
  reproduceSchmidtStudy_res[["filtered_timeseries"]] <- filtered_timeseries
  reproduceSchmidtStudy_res[["total_ProbesetSet"]] <- finalProbesetSet
  #save(reproduceSchmidtStudy_res, file = "temp/reproduceSchmidtStudy_temp.Rdata")
  reproduceSchmidtStudy_res
}

## data folder D:/Dropbox/Dropbox/FBNNet/FBNNet2/output/Experiments_data/leukeamia/differential_cutOffInduction_1_majority_dot_6
## files2<-'D:\\Dropbox\\Dropbox\\FBNNet\\ChildhoodLeukeamiaDataFile\\GSE2677_RAW'
## files3<-'D:\\Dropbox\\Dropbox\\FBNNet\\FBNNet2\\study\\leukeamia_data'
## files1<-c('D:\\Dropbox\\Dropbox\\FBNNet\\ChildhoodLeukeamiaDataFile\\GSE2677_RAW','G:\\Dropbox\\Dropbox\\FBNNet\\Genome_Data\\GSE13670_RAW\\GSE13670_RAW','G:\\Dropbox\\Dropbox\\FBNNet\\Genome_Data\\GSE20489_RAW\\GSE20489_RAW','G:\\Dropbox\\Dropbox\\FBNNet\\Genome_Data\\GSE42088_RAW\\GSE42088_RAW'，'G:\\Dropbox\\Dropbox\\FBNNet\\Genome_Data\\GSE54992_RAW\\GSE54992_RAW'，'G:\\Dropbox\\Dropbox\\FBNNet\\Genome
## Data\\GSE57194_RAW\\GSE57194_RAW')

#'@export
reproduceSchmidtStudy <- function(cellDirectory, 
                                  sortedtimeseries = NULL, 
                                  useGCRMA = FALSE, 
                                  cutOffInduction = 0.7, 
                                  cutOffRepression = 0.7, 
                                  majority = 7,
                                  targetsamples = c("B-ALL-13", 
                                                    "B-ALL-17", 
                                                    "B-ALL-24", 
                                                    "B-ALL-31", 
                                                    "B-ALL-32", 
                                                    "B-ALL-33", 
                                                    "B-ALL-37", 
                                                    "B-ALL-38", 
                                                    "B-ALL-40", 
                                                    "B-ALL-43", 
                                                    "T-ALL-2", 
                                                    "T-ALL-20",
                                                    "T-ALL-25")) {

    # read affy files and normalized by RMA
    if (is.null(sortedtimeseries)) {
        sortedtimeseries2 <- convertAffyRawDataIntoNormalizedStructureData(cellDirectory, useGCRMA = useGCRMA)
    }
    
    sortedtimeseries <- sortedtimeseries2[targetsamples]
    
    cond <- sapply(sortedtimeseries, function(entry) !is.null(entry))
    sortedtimeseries <- sortedtimeseries[cond]
    
    #originSchmidtEvalues <- read.delim("output/Schmidt2006/GSE2677_E-Values-Patients.csv", header = TRUE, sep = ",")
    # testMat <- originSchmidtEvalues[, 5:7]
    # rownames(testMat) <- originSchmidtEvalues[, 1]
    # colnames(testMat) <- c(0, 8, 24)
    
    # try to reproduce the result of the paper 'Identification of glucocorticoid-response genes in children with acute lymphoblastic leukemia', written by
    # Stefan Schmidt schmidtMvalues = read.table('study/Schmidt2006/GSE2677_M-Values-Patients.xls', header = TRUE)
    schmidtMvalues <- read.delim("study/Schmidt2006/GSE2677_M-Values-Patients.csv", header = TRUE, sep = ",")
    schmidtMvalues6To0 <- schmidtMvalues[, c(1, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28)]
    schmidtMvalues24To0 <- schmidtMvalues[, c(1, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29)]
    
    schmidtMvalues24To6Or8 <- read.delim("study/Schmidt2006/GSE2677_M-Values-Patients-6-24.csv", header = TRUE, sep = ",")
    schmidtMvalues24To6Or8 <- schmidtMvalues24To6Or8[, c(1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16)]
    
    # remove the id column
    schmidtMvalues6To0_temp <- sapply(schmidtMvalues6To0[, -1], as.numeric)
    schmidtMvalues24To0_temp <- sapply(schmidtMvalues24To0[, -1], as.numeric)
    schmidtMvalues24To6Or8_temp <- sapply(schmidtMvalues24To6Or8[, -1], as.numeric)
    
    rownames(schmidtMvalues6To0_temp) <- schmidtMvalues6To0[, 1]
    rownames(schmidtMvalues24To0_temp) <- schmidtMvalues24To0[, 1]
    rownames(schmidtMvalues24To6Or8_temp) <- schmidtMvalues24To6Or8[, 1]
    
    schmidtMvalues6To0 <- schmidtMvalues6To0_temp
    schmidtMvalues24To0 <- schmidtMvalues24To0_temp
    schmidtMvalues24To6Or8 <- schmidtMvalues24To6Or8_temp
    
    # experiment with the criteria of logRatio=0.7(+-), majority 6 of 13
    cat("Experiment with the criteria of cutOffInduction=", cutOffInduction, " & cutOffRepression=", cutOffRepression, " (+-), majority=", majority, " of 13")
    print(paste("Experiment with cutOffInduction=", cutOffInduction, " and the majority of ", majority, sep = ""))
    # induced 0 -> 6
    newRownames <- c()
    newMatinduced <- do.call(rbind, lapply(1:nrow(schmidtMvalues6To0), function(rowindex, logRatio, majority, schmidtMvalues6To0) {
        # get row vector
       rowname <- rownames(schmidtMvalues6To0)[rowindex]
        rowvector <- schmidtMvalues6To0[rowindex, ]
        totalSamples <- length(rowvector)
        numOfMatchCriteria <- length(which(as.numeric(rowvector) >= logRatio))
        
        if (numOfMatchCriteria >= majority) {
            return(c(rownames(schmidtMvalues6To0)[rowindex], schmidtMvalues6To0[rowindex, ], majority = majority))
        }else {
           return (NULL)
        }
    }, cutOffInduction, majority, schmidtMvalues6To0))
    
    print("6To0 induced:")
    Probeset6To0induced <- newMatinduced[, 1]
    if (is.null(Probeset6To0induced)) 
        Probeset6To0induced <- c()
    
    Genes6To0induced <- unlist(mapProbesetNames(Probeset6To0induced))
    deduplicated6To0induced <- unique(Genes6To0induced)
    print(deduplicated6To0induced)
    print(paste("In total =", length(Probeset6To0induced), " unqiue =", length(deduplicated6To0induced), sep = ""))
    
    # repressed
    newRownames <- c()
    newMatrepressed <- do.call(rbind, lapply(1:nrow(schmidtMvalues6To0), function(rowindex, logRatio, majority, schmidtMvalues6To0) {
        # get row vector
        rowvector <- schmidtMvalues6To0[rowindex, ]
        totalSamples <- length(rowvector)
        numOfMatchCriteria <- length(which(as.numeric(rowvector) <= (-1 * logRatio)))

        if (numOfMatchCriteria >= majority) {
            return(c(rownames(schmidtMvalues6To0)[rowindex], schmidtMvalues6To0[rowindex, ], majority = majority))
        }else {
          return (NULL)
        }
    }, cutOffRepression, majority, schmidtMvalues6To0))
    
    print(paste("Experiment with cutOffRepression=", cutOffRepression, " and the majority of ", majority, " out of 13", sep = ""))
    print("6To0 repressed:")
    Probeset6To0repressed <- newMatrepressed[, 1]
    if (is.null(Probeset6To0repressed)) 
        Probeset6To0repressed <- c()
    
    Genes6To0repressed <- unlist(mapProbesetNames(Probeset6To0repressed))
    deduplicated6To0repressed <- unique(Genes6To0repressed)
    print(deduplicated6To0repressed)
    print(paste("In total =", length(Probeset6To0repressed), " unqiue =", length(deduplicated6To0repressed), sep = ""))
    
    # 24-0 induced
    newRownames <- c()
    newMat2induced <- do.call(rbind, lapply(1:nrow(schmidtMvalues24To0), function(rowindex, logRatio, majority, schmidtMvalues24To0) {
        # get row vector
        rowvector <- schmidtMvalues24To0[rowindex, ]
        totalSamples <- length(rowvector)
        numOfMatchCriteria <- length(which(as.numeric(rowvector) >= logRatio))

        if (numOfMatchCriteria >= majority) {
            return(c(rownames(schmidtMvalues24To0)[rowindex], schmidtMvalues24To0[rowindex, ], majority = majority))
        }else {
          return (NULL)
        }
    }, cutOffInduction, majority, schmidtMvalues24To0))
    
    print("24To0 induced:")
    Probeset24To0induced <- newMat2induced[, 1]
    if (is.null(Probeset24To0induced)) 
        Probeset24To0induced <- c()
    
    Genes24To0induced <- unlist(mapProbesetNames(Probeset24To0induced))
    deduplicated24To0induced <- unique(Genes24To0induced)
    print(deduplicated24To0induced)
    print(paste("In total =", length(Probeset24To0induced), " unqiue =", length(deduplicated24To0induced), sep = ""))
    
    # repressed
    newRownames <- c()
    newMat2repressed <- do.call(rbind, lapply(1:nrow(schmidtMvalues24To0), function(rowindex, logRatio, majority, schmidtMvalues24To0) {
        # get row vector
        rowvector <- schmidtMvalues24To0[rowindex, ]
        totalSamples <- length(rowvector)
        numOfMatchCriteria <- length(which(as.numeric(rowvector) <= (-1 * logRatio)))
        
        if (numOfMatchCriteria >= majority) {
            return(c(rownames(schmidtMvalues24To0)[rowindex], schmidtMvalues24To0[rowindex, ], majority = majority))
        }else {
          return (NULL)
        }
    }, cutOffRepression, majority, schmidtMvalues24To0))
    
    print("24To0 repressed:")
    Probeset24To0repressed <- newMat2repressed[, 1]
    if (is.null(Probeset24To0repressed)) 
        Probeset24To0repressed <- c()
    
    Genes24To0repressed <- unlist(mapProbesetNames(Probeset24To0repressed))
    deduplicated24To0repressed <- unique(Genes24To0repressed)
    print(deduplicated24To0repressed)
    print(paste("In total =", length(Probeset24To0repressed), " unqiue =", length(deduplicated24To0repressed), sep = ""))
    
    # 24To6/8 schmidtMvalues24To6Or8 induced
    newRownames <- c()
    newMat3induced <- do.call(rbind, lapply(1:nrow(schmidtMvalues24To6Or8), function(rowindex, logRatio, majority, schmidtMvalues24To6Or8) {
        # get row vector
        rowvector <- schmidtMvalues24To6Or8[rowindex, ]
        totalSamples <- length(rowvector)
        numOfMatchCriteria <- length(which(as.numeric(rowvector) >= logRatio))
  
        if (numOfMatchCriteria >= majority) {
            return(c(rownames(schmidtMvalues24To6Or8)[rowindex], schmidtMvalues24To6Or8[rowindex, ], majority = majority))
        }else {
          return (NULL)
        }
    }, cutOffInduction, majority, schmidtMvalues24To6Or8))
    
    print("24To6/8 induced:")
    Probeset24To6Or8induced <- newMat3induced[, 1]
    if (is.null(Probeset24To6Or8induced)) 
        Probeset24To6Or8induced <- c()
    
    Genes24To6Or8induced <- unlist(mapProbesetNames(Probeset24To6Or8induced))
    deduplicated24To6Or8induced <- unique(Genes24To6Or8induced)
    print(deduplicated24To6Or8induced)
    print(paste("In total =", length(Probeset24To6Or8induced), " unqiue =", length(deduplicated24To6Or8induced), sep = ""))
    
    # repressed
    newRownames <- c()
    newMat3repressed <- do.call(rbind, lapply(1:nrow(schmidtMvalues24To6Or8), function(rowindex, logRatio, majority, schmidtMvalues24To6Or8) {
        # get row vector
        rowvector <- schmidtMvalues24To6Or8[rowindex, ]
        totalSamples <- length(rowvector)
        numOfMatchCriteria <- length(which(as.numeric(rowvector) <= (-1 * logRatio)))
        
        if (numOfMatchCriteria >= majority) {
            return(c(rownames(schmidtMvalues24To6Or8)[rowindex], schmidtMvalues24To6Or8[rowindex, ], majority = majority))
        }else {
          return (NULL)
        }
    }, cutOffRepression, majority, schmidtMvalues24To6Or8))
    
    print("24To6/8 repressed:")
    Probeset24To6or8repressed <- newMat3repressed[, 1]
    if (is.null(Probeset24To6or8repressed)) 
        Probeset24To6or8repressed <- c()
    
    Genes24To6or8repressed <- unlist(mapProbesetNames(Probeset24To6or8repressed))
    deduplicated24To6or8repressed <- unique(Genes24To6or8repressed)
    print(deduplicated24To6or8repressed)
    print(paste("In total =", length(Probeset24To6or8repressed), " unqiue =", length(deduplicated24To6or8repressed), sep = ""))
    
    ############################### output ################################################################################ todo: add the difference as well, plus the primise genes for
    ############################### further mining
    finalSet <- unique(c(deduplicated6To0induced, 
                         deduplicated6To0repressed, 
                         deduplicated24To0induced, 
                         deduplicated24To0repressed, 
                         deduplicated24To6Or8induced, 
                         deduplicated24To6or8repressed))
    print(paste("Total unique genes cross all time points are:", length(finalSet), sep = ""))
    print(finalSet)
    finalProbesetSet <- unique(c(Probeset6To0induced, Probeset6To0repressed, Probeset24To0induced, Probeset24To0repressed, Probeset24To6Or8induced, Probeset24To6or8repressed))
    filtered_timeseries <- lapply(sortedtimeseries2, function(subdata) subdata[rownames(subdata) %in% finalProbesetSet, ])
    
    reproduceSchmidtStudy_res <- list()
    reproduceSchmidtStudy_res[["CombinedGeneSet"]] <- finalSet
    reproduceSchmidtStudy_res[["SchmidtGenes6To0induced"]] <- deduplicated6To0induced
    reproduceSchmidtStudy_res[["SchmidtGenes6To0repressed"]] <- deduplicated6To0repressed
    reproduceSchmidtStudy_res[["SchmidtGenes24To0induced"]] <- deduplicated24To0induced
    reproduceSchmidtStudy_res[["SchmidtGenes24To0repressed"]] <- deduplicated24To0repressed
    reproduceSchmidtStudy_res[["SchmidtGenes24To6Or8induced"]] <- deduplicated24To6Or8induced
    reproduceSchmidtStudy_res[["SchmidtGenes24To6or8repressed"]] <- deduplicated24To6or8repressed
    reproduceSchmidtStudy_res[["timeseries_original"]] <- sortedtimeseries2
    reproduceSchmidtStudy_res[["filtered_timeseries"]] <- filtered_timeseries
    reproduceSchmidtStudy_res[["total_ProbesetSet"]] <- finalProbesetSet
    reproduceSchmidtStudy_res[["Schmidt6To0inducedMAT"]] <- newMatinduced
    reproduceSchmidtStudy_res[["Schmidt6To0repressedMAT"]] <- newMatrepressed
    reproduceSchmidtStudy_res[["Schmidt24To0inducedMAT"]] <- newMat2induced
    reproduceSchmidtStudy_res[["Schmidt24To0repressedMAT"]] <- newMat2repressed
    reproduceSchmidtStudy_res[["Schmidt24To6Or8inducedMAT"]] <- newMat3induced
    reproduceSchmidtStudy_res[["Schmidt24To6or8repressedMAT"]] <- newMat3repressed
    #save(reproduceSchmidtStudy_res, file = "temp/reproduceSchmidtStudy_temp.Rdata")
    reproduceSchmidtStudy_res
}

#'@export
leukeamia_study_with_differential_output <- function(schmidts_realanlysis_output, 
                                                 method = c("average", "kmeans", "edgeDetector", "scanStatistic"), 
                                                 maxK = 4, 
                                                 temporal = 2,
                                                 target_genes = c()) {
  sortedtimeseries_leukaemia <- schmidts_realanlysis_output$filtered_timeseries
  ## shouldn't have duplicate row names
  convertedgenedagta <- convertTimeseriesProbsetNameToGeneName(sortedtimeseries_leukaemia)$convert_data
  finalSet <- unique(schmidts_realanlysis_output$CombinedGeneSet)
  if( length(target_genes) > 0) {
    finalSet = target_genes
  }

  
  # remove duplicate
  convertedgenedagta <- lapply(convertedgenedagta, function(subdata)subdata[which(rownames(subdata) %in% finalSet), ])
  
    
    print(rownames(convertedgenedagta[[1]]))
    
    # timeseries_RMA_LogRatio_0.7<<-discreteTimeSeries(convertedgenedagta,method='average')
    
    if (method == "average") {
        timeseries_RMA_LogRatio_0.7 <- discreteTimeSeries(convertedgenedagta, method = "average")
    } else {
        require(BoolNet)
        timeseries_RMA_LogRatio_0.7 <- BoolNet::binarizeTimeSeries(convertedgenedagta, method = method)$binarizedMeasurements
    }
    
    # membexp=3 more fussy
    genes <- rownames(timeseries_RMA_LogRatio_0.7[[1]])
    # build all cubes for all clusters
    cubeLeukaemia_RMA_LogRatio_temp <- constructFBNCube(genes, genes, timeseries_RMA_LogRatio_0.7, maxK, temporal, TRUE)
    gc()
    #.rs.restartR()
    save(cubeLeukaemia_RMA_LogRatio_temp, file = "temp/cubeLeukaemia_RMA_LogRatio_temp.Rdata")
    totalNetworks_RMA_LogRatio_0.7 <- mineFBNNetwork(cubeLeukaemia_RMA_LogRatio_temp, useParallel = TRUE)
    
    leukeamia_study_with_schmidts_output_res <- list(schmidts_realanlysis_output = schmidts_realanlysis_output, cube = cubeLeukaemia_RMA_LogRatio_temp, network = totalNetworks_RMA_LogRatio_0.7, 
        timeseries = timeseries_RMA_LogRatio_0.7)
    #save(leukeamia_study_with_schmidts_output_res, file = "temp/leukeamia_study_with_schmidts_output_temp.Rdata")
    leukeamia_study_with_schmidts_output_res
}

#'@export
leukeamia_study_with_schmidts_output_cluster <- function(schmidts_realanlysis_output, 
                                                         method = c("average", "kmeans", "edgeDetector", "scanStatistic"), 
                                                         maxK = 4, 
                                                         temporal = 2, 
                                                         minElementInCluster = 10, 
                                                         maxElementInCluster = 20) {
    sortedtimeseries_leukaemia <- schmidts_realanlysis_output$filtered_timeseries
    convertedgenedagta <- convertTimeseriesProbsetNameToGeneName(sortedtimeseries_leukaemia)$convert_data
    
    finalSet <- unique(schmidts_realanlysis_output$CombinedGeneSet)
    print(finalSet)
    
    # remove duplicate
    convertedgenedagta <- lapply(convertedgenedagta, function(subdata) subdata[which(rownames(subdata) %in% finalSet), ])
    
    print(rownames(convertedgenedagta[[1]]))
    
    # timeseries_RMA_LogRatio_0.7<<-discreteTimeSeries(convertedgenedagta,method='average')
    
    if (method == "average") {
        timeseries_RMA_LogRatio_0.7 <- discreteTimeSeries(convertedgenedagta, method = "average")
    } else {
        require(BoolNet)
        timeseries_RMA_LogRatio_0.7 <- BoolNet::binarizeTimeSeries(convertedgenedagta, method = method)$binarizedMeasurements
    }
    
    totalNetworks_RMA_LogRatio_0.7 <- generateFBMNetwork(timeseries_data = timeseries_RMA_LogRatio_0.7, 
                                  maxK = 4, 
                                  max_deep_temporal = 2, 
                                  useParallel = TRUE,
                                  maxGenesForSingleCube = 10,
                                  parallel_on_group = TRUE)
    
    #getClusteredTimeseries_RMA_LogRatio_0.7 <- dividedDataIntoSubgroups(timeseries_RMA_LogRatio_0.7, maxElementInCluster)
    # build all cubes for all clusters build all cubes for all clusters
    #cubeLeukaemia_RMA_LogRatio_temp <- constructFBNCubeAndNetworkInClusters_combine(getClusteredTimeseries_RMA_LogRatio_0.7, maxK = maxK, temporal = temporal, useParallel = TRUE)
    # merge all networks? seperate by clusters?
    #networks_RMA_LogRatio_0.7 <- mergeClusterNetworks(cubeLeukaemia_RMA_LogRatio_temp)
    #totalNetworks_RMA_LogRatio_0.7 <- filterNetworkConnections(networks_RMA_LogRatio_0.7)
    
    #save(cubeLeukaemia_RMA_LogRatio_temp, file = "temp/cubeLeukaemia_RMA_LogRatio_temp_cluster.Rdata")
    leukeamia_study_with_schmidts_output_res <- list(schmidts_realanlysis_output = schmidts_realanlysis_output, network = totalNetworks_RMA_LogRatio_0.7, timeseries = timeseries_RMA_LogRatio_0.7)
    #save(leukeamia_study_with_schmidts_output_res, file = "temp/leukeamia_study_with_schmidts_output_temp.Rdata")
    leukeamia_study_with_schmidts_output_res
}

## files2<-'D:\\Dropbox\\Dropbox\\FBNNet\\ChildhoodLeukeamiaDataFile\\GSE2677_RAW'
## files1<-c('D:\\Dropbox\\Dropbox\\FBNNet\\ChildhoodLeukeamiaDataFile\\GSE2677_RAW','D:\\Dropbox\\Dropbox\\FBNNet\\Genome
## Data\\GSE13670_RAW\\GSE13670_RAW','D:\\Dropbox\\Dropbox\\FBNNet\\Genome
## Data\\GSE20489_RAW\\GSE20489_RAW','D:\\Dropbox\\Dropbox\\FBNNet\\Genome
## Data\\GSE42088_RAW\\GSE42088_RAW'，'D:\\Dropbox\\Dropbox\\FBNNet\\Genome
## Data\\GSE54992_RAW\\GSE54992_RAW'，'D:\\Dropbox\\Dropbox\\FBNNet\\Genome Data\\GSE57194_RAW\\GSE57194_RAW')
#'@export
leukeamia_study_cluster <- function(cellDirectory, 
                                    cutOffInduction = 0.7, 
                                    cutOffRepression = 0.7, 
                                    majority = 0.5, 
                                    sortedtimeseries = NULL, 
                                    useGCRMA = FALSE, 
    method = c("average", "kmeans", "edgeDetector", "scanStatistic")) {
    # read affy files and normalized by RMA
    if (is.null(sortedtimeseries)) {
        sortedtimeseries <- convertAffyRawDataIntoNormalizedStructureData(cellDirectory, useGCRMA = useGCRMA)
    }
    
    targetsamples <- c("B-ALL-13", "B-ALL-17", "B-ALL-24", "B-ALL-31", "B-ALL-32", "B-ALL-33", "B-ALL-37", "B-ALL-38", "B-ALL-40", "B-ALL-43", "T-ALL-2", "T-ALL-20", 
        "T-ALL-25")
    sortedtimeseries_leukaemia <- sortedtimeseries[targetsamples]
    cond <- sapply(sortedtimeseries_leukaemia, function(entry) !is.null(entry))
    sortedtimeseries_leukaemia <- sortedtimeseries_leukaemia[cond]
    
    # Generate cube 0.7 = 2.0 ^ 0.7
    probesets <- rownames(sortedtimeseries_leukaemia[[1]])
    probesetGeneNameMappings <- mapProbesetNames(probesets)
    
    diffgenes_RMA <- identifyDifferentiallyExpressedGenes(sortedtimeseries_leukaemia, cutOffInduction = cutOffInduction, cutOffRepression = cutOffRepression, 
        majority = majority, probesetGeneNameMappings = probesetGeneNameMappings)
    
    commonGeneSet1i <- diffgenes_RMA$DifferentialExpression[[1]]$Induced_ProbeID
    commonGeneSet1r <- diffgenes_RMA$DifferentialExpression[[1]]$Repressed__ProbeID
    commonGeneSet2i <- diffgenes_RMA$DifferentialExpression[[2]]$Induced_ProbeID
    commonGeneSet2r <- diffgenes_RMA$DifferentialExpression[[2]]$Repressed__ProbeID
    commonGeneSet3i <- diffgenes_RMA$DifferentialExpression[[3]]$Induced_ProbeID
    commonGeneSet3r <- diffgenes_RMA$DifferentialExpression[[3]]$Repressed__ProbeID
    
    finalSet <- unique(c(commonGeneSet1i, commonGeneSet1r, commonGeneSet2i, commonGeneSet2r, commonGeneSet3i, commonGeneSet3r))
    print(finalSet)
    
    # remove duplicate
    subsetgenedata <- lapply(sortedtimeseries, function(subdata) subdata[rownames(subdata) %in% finalSet, ])
    convertedgenedagta <- convertTimeseriesProbsetNameToGeneName(subsetgenedata)$convert_data
    print(rownames(convertedgenedagta[[1]]))
    
    if (method == "average") {
        timeseries_RMA_LogRatio_0.7 <- discreteTimeSeries(convertedgenedagta, method = "average")
    } else {
        require(BoolNet)
        timeseries_RMA_LogRatio_0.7 <- BoolNet::binarizeTimeSeries(convertedgenedagta, method = method)$binarizedMeasurements
    }
    
    # membexp=3 more fussy
    getClusteredTimeseries_RMA_LogRatio_0.7 <- clusterdDiscreteData(convertedgenedagta,timeseries_RMA_LogRatio_0.7, 30)
    # build all cubes for all clusters
    cubeLeukaemia_RMA_LogRatio_0.7 <- constructFBNCubeAndNetworkInClusters(getClusteredTimeseries_RMA_LogRatio_0.7, 4, 1, TRUE)
    # merge all networks? seperate by clusters?
    networks_RMA_LogRatio_0.7 <- mergeClusterNetworks(cubeLeukaemia_RMA_LogRatio_0.7)
    totalNetworks_RMA_LogRatio_0.7 <- filterNetworkConnections(networks_RMA_LogRatio_0.7)
    
    diffgenes <- list(Induction_0_to_6_or_8 = commonGeneSet1i, Repression_0_to_6_or_8 = commonGeneSet1r, Induction_0_to_24 = commonGeneSet2i, Repression_0_to_24 = commonGeneSet2r, 
        Induction_6_or_8_to_24 = commonGeneSet3i, Repression_6_or_8_to_24 = commonGeneSet3r)
    list(diff_genes = diffgenes, cube = timeseries_RMA_LogRatio_0.7, network = totalNetworks_RMA_LogRatio_0.7, timeseries = timeseries_RMA_LogRatio_0.7)
}


#' 
#' # use boolnet
#' #'@export
#' leukeamia_study <- function(cellDirectory, 
#'                             cutOffInduction = 1, 
#'                             cutOffRepression = 1, 
#'                             majority = 0.5, 
#'                             sortedtimeseries = NULL, 
#'                             useGCRMA = FALSE, 
#'                             method = c("average", 
#'     "kmeans", "edgeDetector", "scanStatistic")) {
#'     require(BoolNet)
#'     # read affy files and normalized by RMA
#'     if (is.null(sortedtimeseries)) {
#'         sortedtimeseries <- convertAffyRawDataIntoNormalizedStructureData(cellDirectory, useGCRMA = useGCRMA)
#'     }
#'     
#'     targetsamples <- c("B-ALL-13", "B-ALL-17", "B-ALL-24", "B-ALL-31", "B-ALL-32", "B-ALL-33", "B-ALL-37", "B-ALL-38", "B-ALL-40", "B-ALL-43", "T-ALL-2", "T-ALL-20", 
#'         "T-ALL-25")
#'     sortedtimeseries_leukaemia <- sortedtimeseries[targetsamples]
#'     cond <- sapply(sortedtimeseries_leukaemia, function(entry) !is.null(entry))
#'     sortedtimeseries_leukaemia <- sortedtimeseries_leukaemia[cond]
#'     
#'     # Generate cube 0.7 = 2.0 ^ 0.7
#'     cutOffInduction = cutOffInduction
#'     cutOffRepression = cutOffRepression
#'     majority <- majority
#'     probesets <- rownames(sortedtimeseries_leukaemia[[1]])
#'     probesetGeneNameMappings <- mapProbesetNames(probesets)
#'     
#'     diffgenes_RMA <- identifyDifferentiallyExpressedGenes(sortedtimeseries_leukaemia, cutOffInduction = cutOffInduction, cutOffRepression = cutOffRepression, 
#'         majority = majority, probesetGeneNameMappings = probesetGeneNameMappings)
#'     
#'     commonGeneSet1i <- diffgenes_RMA$DifferentialExpression[[1]]$Induced_ProbeID
#'     commonGeneSet1r <- diffgenes_RMA$DifferentialExpression[[1]]$Repressed__ProbeID
#'     commonGeneSet2i <- diffgenes_RMA$DifferentialExpression[[2]]$Induced_ProbeID
#'     commonGeneSet2r <- diffgenes_RMA$DifferentialExpression[[2]]$Repressed__ProbeID
#'     commonGeneSet3i <- diffgenes_RMA$DifferentialExpression[[3]]$Induced_ProbeID
#'     commonGeneSet3r <- diffgenes_RMA$DifferentialExpression[[3]]$Repressed__ProbeID
#'     
#'     finalSet <- unique(c(commonGeneSet1i, commonGeneSet1r, commonGeneSet2i, commonGeneSet2r, commonGeneSet3i, commonGeneSet3r))
#'     print(finalSet)
#'     
#'     # remove duplicate
#'     subsetgenedata <- lapply(sortedtimeseries, function(subdata) subdata[rownames(subdata) %in% finalSet, ])
#'     convertedgenedagta <- convertTimeseriesProbsetNameToGeneName(subsetgenedata)$convert_data
#'     print(rownames(convertedgenedagta[[1]]))
#'     
#'     # timeseries_RMA_LogRatio_0.7<<-discreteTimeSeries(convertedgenedagta,method='average')
#'     
#'     if (method == "average") {
#'         timeseries_RMA_LogRatio_0.7 <- discreteTimeSeries(convertedgenedagta, method = "average")
#'     } else {
#'         require(BoolNet)
#'         timeseries_RMA_LogRatio_0.7 <- BoolNet::binarizeTimeSeries(convertedgenedagta, method = method)$binarizedMeasurements
#'     }
#'     
#'     # membexp=3 more fussy
#'     genes <- rownames(timeseries_RMA_LogRatio_0.7[[1]])
#'     # build all cubes for all clusters
#'     cubeLeukaemia_RMA_LogRatio_0.7 <- constructFBNCube(genes, genes, timeseries_RMA_LogRatio_0.7, 4, 2, TRUE, TRUE)
#'     totalNetworks_RMA_LogRatio_0.7 <- mineFBNNetwork(cubeLeukaemia_RMA_LogRatio_0.7)
#'     
#'     diffgenes <- list(Induction_0_to_6_or_8 = commonGeneSet1i, Repression_0_to_6_or_8 = commonGeneSet1r, Induction_0_to_24 = commonGeneSet2i, Repression_0_to_24 = commonGeneSet2r, 
#'         Induction_6_or_8_to_24 = commonGeneSet3i, Repression_6_or_8_to_24 = commonGeneSet3r)
#'     list(diff_genes = diffgenes, cube = timeseries_RMA_LogRatio_0.7, network = totalNetworks_RMA_LogRatio_0.7, timeseries = timeseries_RMA_LogRatio_0.7)
#'     # merge all networks? seperate by clusters?  networks_RMA_LogRatio_0.7<<-mergeClusterNetworks(cubeLeukaemia_RMA_LogRatio_0.7)
#'     # totalNetworks_RMA_LogRatio_0.7<<-networks_RMA_LogRatio_0.7$TotalNetwork totalGenes_RMA_LogRatio_0.7<<-networks_RMA_LogRatio_0.7$TotalNetwork$genes
#'     # individualNetworks_RMA_LogRatio_0.7<<-networks_RMA_LogRatio_0.7$IndividualNetwork
#' }
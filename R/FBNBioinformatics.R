## if (!requireNamespace("BiocManager", quietly = TRUE))
##   install.packages("BiocManager")

## BiocManager::install("AnnotationDbi")
## BiocManager::install("hgu133plus2.db")


#' Functions related to Bioinformatics
#' 
#' A collection of functions that used to handle 
#' bioinformatics related tasks 
#' 
#' @param files The file folder path or directory that 
#' contains the affy files (*.cel)
#' @param useGCRMA optional, if true it use GCRMA method to normalize
#'  the affy data, otherwise, use the default RMA method.
#'   GCRMA, A bias-corrected RMA; FALSE, use RMA, which is based on
#'    Robust Multi-Chip' average
#' @param cdfname The name of CDF that is associated with the affy data. 
#' @param cellDirectory A directory that contains affy raw data
#' @param outPutEvalue optional to write the result into disk
#' @param rawData the raw data that is the output from the method 
#' getRelatedAffyRawData
#' @param ges_assess_no the GSE identity id, such as GSE35635
#' @author Leshi Chen, leshi, chen@lincolnuni.ac.nz, chenleshi@hotmail.com
#' @keywords Fundamental Boolean Network, Boolean Network, Genetic Regulatory Network
#'
#' @references Chen et al.(2018), Front. Physiol., 25 September 2018, 
#' (\href{https://doi.org/10.3389/fphys.2018.01328}{Front. Physiol.})
#' @references Bioconductor (2019), https://www.bioconductor.org/
#' @examples 
#' exprs <- downloadGenomeDataExprs(GSE35635)
#' exprs
#' @name FBNBioinformatics
NULL 

#' This method is used to construct time series data such as
#'  Leukeamia data from GSE2677 Experimental method
#' @rdname "FBNBioinformatics"
#' @export
convertAffyRawDataIntoNormalizedStructureData <- function(files, 
                                                          useGCRMA = FALSE, 
                                                          cdfname = "HG-U133_Plus_2",
                                                          isNew = TRUE) {
  futile.logger::flog.info(sprintf("Enter convertAffyRawDataIntoNormalizedStructureData zone: 
                                   files=%s, 
                                   useGCRMA=%s, 
                                   cdfname=%s",
                                   files,
                                   useGCRMA,
                                   cdfname))
  typeOfRMA = "RMA"
  if (useGCRMA) {
    typeOfRMA = "GCRMA"
  }
  
  dir.create(file.path(getwd(), "output"), showWarnings = FALSE)
  
  combinedTimeseries <- list()
  for (filename in files) {
    futile.logger::flog.info(paste0("Processing file: "), filename)
    output_normalized_file <- paste("output\\NormalizedEValue", "_", typeOfRMA, ".csv")
    if (!file.exists(output_normalized_file) || isNew) {
      rawdata <- getRelatedAffyRawData(filename, cdfname)
      ## normalize data either using GCRMA or RMA 
      nrawdata <- normalizeTimesereisRawData(rawdata, typeOfRMA)$NormalizedExpData
      write.csv(nrawdata, file = output_normalized_file)
    } else {
      nrawdata <- read.csv(file = output_normalized_file, 
                           header = TRUE, 
                           sep = ",", 
                           check.names = FALSE)
      nrawdata2 <- nrawdata[, -1]
      rownames(nrawdata2) <- nrawdata[, 1]
      nrawdata <- as.matrix(nrawdata2)
    }
    
    # convert into sample->time points i.e., group by samples and order by time points
    stimeseries <- convertIntoSampleTimeSeries(nrawdata)
    sortedtimeseries <- reorderSampleTimeSeries(stimeseries)
    
    # combine all time series from different source/files
    combinedTimeseries[[length(combinedTimeseries) + 1]] <- sortedtimeseries
  }
  futile.logger::flog.info("Leave convertAffyRawDataIntoNormalizedStructureData zone.")
  dissolve(combinedTimeseries)
}

## try http:// if https:// URLs are not supported, use follow to install required bioconductor packages source('https://bioconductor.org/biocLite.R')
## biocLite(c('affy', 'limma')) help(package='limma') installed.packages() %>% .[, c('Package', 'LibPath')] to find all libary path test
## 'C:\\Users\\chenl\\Dropbox\\FBNNet\\ChildhoodLeukeamiaDataFile\\GSE2677_RAW'
##'HG-U133_Plus_2cdf'

#'A function to read Affy row data
#' @rdname "FBNBioinformatics"
#'@export
getRelatedAffyRawData <- function(cellDirectory, 
                                  cdfname = "HG-U133_Plus_2", 
                                  outPutEvalue = FALSE) {
  futile.logger::flog.info(sprintf("Enter getRelatedAffyRawData zone: cellDirectory=%s, cdfname=%s, outPutEvalue=%s",
                                   cellDirectory,
                                   cdfname,
                                   outPutEvalue))
  dir.create(file.path(getwd(), "output"), showWarnings = FALSE)
  # 1 Read in probe level data The affy package will automatically download the appropriate array annotation when you require it. However, if you wish you may
  # download and install the cdf environment you need from http://www.bioconductor.org/packages/release/data/annotation/ manually. If there is no cdf
  # environment currently built for your particular chip and you have access to the CDF file then you may use the makecdfenv package to create one yourself.
  # To make the cdf packaes, Microsoft Windows users will need to use the tools described here: http://cran.r-project.org/bin/windows/rw-FAQ.html
  
  # FIRST solution mycdf <- read.cdffile(cdfFileName) source('http://www.bioconductor.org/biocLite.R') biocLite('affy') biocLite(cdfname) require('affy')
  # library(affy) biocLite('HG-U133_Plus_2cdf')
  listCellFiles <- dir(path = cellDirectory, pattern = "*\\.CEL", full.names = TRUE)
  affydata <- affy::ReadAffy(filenames = listCellFiles)
  # indicate you want to use the custom cdf If you don't specify the cdfname, BioConductor will use the default Affymetrix cdf.
  affydata@cdfName = cdfname
  
  # raw expression data
  expdata <- affy::exprs(affydata)
  if (outPutEvalue) {
    write.csv(expdata, file = paste("output/EValue", ".csv"))
  }
  
  
  samp <- affy::sampleNames(affydata)
  probes <- affy::featureNames(affydata)
  
  res <- list()
  res[[1]] <- affydata
  names(res)[[1]] <- "AffyBatchObject"
  
  res[[2]] <- expdata
  names(res)[[2]] <- "ExpressionData"
  
  res[[3]] <- samp
  names(res)[[3]] <- "SampleName"
  
  res[[4]] <- probes
  names(res)[[4]] <- "Probes"
  futile.logger::flog.info("Leave getRelatedAffyRawData zone.")
  res
}

#' A benchmark type function to test a complete process of FBN model
#' Step 2, normalizing data
#' @rdname "FBNBioinformatics"
#' @export
normalizeTimesereisRawData <- function(rawData, 
                                       method = c("RMA", "GCRMA", "MAS5")) {
  futile.logger::flog.info(sprintf("Enter normalizeTimesereisRawData zone: length of rawData=%s, method=%s",
                                   length(rawData),
                                   method))
  # source('http://www.bioconductor.org/biocLite.R') biocLite('affyPLM') require(affyPLM) Normalizing Data The Affy package has implementations of a number of
  # normalization methods for single-channel arrays. This includes (among others):
  #- mas5() - Affymetrix's own normalization program
  #- rma() - 'Robust Multi-Chip' average
  #- gcrma() - A bias-corrected RMA
  # GCRMA is good but takes significantly longer than RMA, so RMA is the most commonly used cat('test') rawData <- as.numeric(as.character(rawData[[1]]))
  # nvals <- rma(rawData) bydefault all normalized data has been log2 processed
  dir.create(file.path(getwd(), "output"), showWarnings = FALSE)
  if (method == "GCRMA") {
    nvals <- gcrma::gcrma(rawData$AffyBatchObject)
    ned <- Biobase::exprs(nvals)
  } else if (method == "MAS5") {
    nvals <- affy::mas5(rawData$AffyBatchObject)
    ned <- log2(Biobase::exprs(nvals))  #scale to log2 based
    
  } else {
    nvals <- affy::rma(rawData$AffyBatchObject)  #, normalize = TRUE, background = TRUE, bgversion = 2, destructive = TRUE)
    ned <- Biobase::exprs(nvals)
  }
  
  # normalised expression data
  
  nsamp <- rawData$SampleName
  nprobes <- rownames(ned)
  
  # check the feature name
  if (!identical(rawData$Probes, nprobes)) {
    write.csv(rawData$Probes, file = paste("output/rawData_probes", ".csv"))
    write.csv(nprobes, file = paste("output/normalized_probes", ".csv"))
    stop("The features of the data is not identical to the features of normalized data")
  }
  
  res <- list()
  res[[1]] <- nvals
  names(res)[[1]] <- "NormalizedObject"
  
  res[[2]] <- ned
  names(res)[[2]] <- "NormalizedExpData"
  
  res[[3]] <- nsamp
  names(res)[[3]] <- "SampleName"
  
  res[[4]] <- nprobes
  names(res)[[4]] <- "FeatureNames"
  futile.logger::flog.info("Leave normalizeTimesereisRawData zone.")
  res
}


#' Download GenomeData Exprs
#' 
#' @rdname "FBNBioinformatics"
#' @export
downloadGenomeDataExprs <- function(ges_assess_no) {
  futile.logger::flog.info(sprintf("Enter downloadGenomeDataExprs zone: ges_assess_no=%s",
                                   ges_assess_no))
    # use gse2553 <- getGEO('GSE2677',GSEMatrix=TRUE) to download GSE2677 data for example. To detect time of day dependent gene expression in human epidermis
    # suction blister samples from 20 healthy subjects were obtained at three different time points throughout the day. RNA from 20 subjects were used to
    # perform whole genome microarray analysis. Microarrays from 19 subjects showed sufficient quality to perform analysis for differential gene expression. We
    # detected significant differential expression levels for several canonical clock genes such as Per1, Per2, Per3, Bmal1 and Rev-Erb_alpha throughout the
    # day. In total we identified 294 genes that showed significant circadian gene expression including several transcription factors and rate limiting enzymes.
    # To our knowledge this is the first study to investigate genome wide circadian gene expression in human epidermis.
    # 1) gathering data from gse35635 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35635
    ges_assess <- GEOquery::getGEO(ges_assess_no, GSEMatrix = TRUE)
    eset_gse <- ges_assess[[1]]
    futile.logger::flog.info("Leave downloadGenomeDataExprs zone.")
    Biobase::assayDataElement(eset_gse, "exprs")
}

#'This method is used to map probeset names with gene names
#'@param probesets A vector of probeset names
#'@param names_mapped the names that have been mapped/processed
#' @export
mapProbesetNames <- function(probesets, 
                             names_mapped = NULL) {
    data("DAVID_Gene_List")
    if (is.null(probesets)) 
      return(c())
    
    if (length(probesets) == 0) 
      return(c())
    
    if (is.null(names_mapped)) 
      probesetMapping <- mapToGeneNameWithhgu133plus2Db(probesets) 
    else 
      probesetMapping <- names_mapped
  
    names(DAVID_Gene_List)[[1]] = "Probeset"
    probesetGeneNameMappings <- lapply(probesets, function(name) {
      geneNames <- convert_probeset_to_gene(probesetMapping, name)$SYMBOL
      geneNames_DAVID <- DAVID_Gene_List[DAVID_Gene_List$Probeset %in% c(name), ]$Symbol
      if (length(geneNames) > 0) {
        return(geneNames)
      } else if (!is.null(geneNames_DAVID)) {
        return(as.character(geneNames_DAVID))
      } else {
        return("N/A")
      }
    })
    names(probesetGeneNameMappings) <- probesets
    probesetGeneNameMappings
}

#'A benchmark type function to test a complete process of FBN model
#'All sub time series must contain the same number of timepoints
#'Step 3, identify significantly expressed genes, which are strongly related with the samples / study purposes
#'@param orderSampleTimeSeries A sorted time series data, which is the output of the method reorderSampleTimeSeries
#'@param cutOffInduction a threshold that identify genes as folds. If the cutOffInduction is 2, the differential genes are identified based on 2 folds
#'@param cutOffRepression a threshold that identify genes as folds. If the cutOffRepression is 2, the differential genes are identified based on 2 folds
#'@param majority A criteria that make a gene as differential
#'@param needLog2scale If it is true, then all gene values will be processed using log2
#'@param probesetGeneNameMappings gene mapping file
#' @export
identifyDifferentiallyExpressedGenes <- function(orderSampleTimeSeries, 
                                                 cutOffInduction = 1, 
                                                 cutOffRepression = 1, 
                                                 majority = 7, 
                                                 needLog2scale = FALSE, 
                                                 probesetGeneNameMappings = NULL, 
                                                 nameTab = "RMA") {
  futile.logger::flog.info(sprintf("Enter identifyDifferentiallyExpressedGenes zone:cutOffInduction=%s, cutOffRepression=%s, majority=%s, needLog2scale=%s, nameTab=%s",
                                   cutOffInduction,
                                   cutOffRepression,
                                   majority,
                                   needLog2scale,
                                   nameTab))
  
  originalData <- orderSampleTimeSeries
  cutOffInduction <- cutOffInduction
  cutOffRepression <- cutOffRepression
  print(paste("cutOffInduction=", cutOffInduction, "; cutOffRepression=", cutOffRepression, "; majority=", majority, sep = ""))
  
  # The expression values from RMA are log2 transformed, so to calculate the log ratio you simply subtract one from the other. log2(x) - log2(y) = log2(x/y)
  # Note here the the log ratio gives you the fold change directly. A log ratio of 1 = 2-fold up regulated and a log ratio of -1 = 2-fold down regulated (when
  # comparing x vs y).
  folddata <- lapply(originalData, function(input) {
    cols <- colnames(input)
    if (length(cols) < 2) {
      stop("The input must have columns more than one")
    }
    
    # get initial data
    if (needLog2scale) {
      res <- round(log2(input[, 2]/input[, 1]), 5)
    } else {
      # already in log 2 scare when normalized
      res <- (input[, 2] - input[, 1])
    }
    
    namesOfDiff <- paste("2", "To", "1", sep = "", collapse = "")
    
    # loop through other combinations
    for (i in seq_along(cols)) {
      if (i > 2) {
        k <- 1
        while (k <= (i - 1)) {
          
          if (needLog2scale) {
            res <- cbind(res, round(log2(input[, i]/input[, k]), 5))
          } else {
            # already in log 2 scare when normalized
            res <- cbind(res, (input[, i] - input[, k]))
          }
          
          namesOfDiff <- c(namesOfDiff, paste(i, "To", k, sep = "", collapse = ""))
          k <- k + 1
        }
      }
    }
    
    if (is.vector(res)) {
      res <- matrix(res, length(res), byrow = TRUE)
    }
    colnames(res) <- namesOfDiff
    res
  })
  
  diffExpressed <- list()
  # anovation
  first_mat <- folddata[[1]]
  colnamelist <- colnames(first_mat)
  rownamelist <- rownames(first_mat)
  
  # map probeID with gene name
  if (is.null(probesetGeneNameMappings)) {
    probesetGeneNameMappings <- mapProbesetNames(rownamelist)
  }
  dir.create(file.path(getwd(), "output"), showWarnings = FALSE)
  dir.create(file.path(getwd(), "output", "differential"), showWarnings = FALSE)
  ## mats <- list()
  for (i in seq_along(colnamelist)) {
    tempMat <- do.call(cbind, lapply(folddata, function(subList) subList[, i]))
    ## mats[i] <- tempMat
    write.csv(tempMat, file = paste("output/differential/MValue_", colnamelist[i], "_", nameTab, ".csv"))
    rownamelist <- rownames(tempMat)
    totalSamples <- ncol(tempMat)
    # print(each column for comparision)
    lapply(1:ncol(tempMat), function(colIndex, logRatio) {
      colvalues <- tempMat[, colIndex]
      names(colvalues) <- rownames(tempMat)
      colNames <- colnames(tempMat)[colIndex]
      colvalues <- colvalues[which(colvalues > logRatio)]
      filteredNames <- names(colvalues)
      convertedGeneNames <- probesetGeneNameMappings[filteredNames]
      res <- cbind(colvalues, convertedGeneNames)
      colnames(res) <- c("M Value", "Gene Name")
      rownames(res) <- filteredNames
      # write.csv(res, file = paste("output/differential/Induction_", colNames, colnamelist[i], ".csv"))
    }, cutOffInduction)
    
    lapply(1:ncol(tempMat), function(colIndex, logRatio) {
      colvalues <- tempMat[, colIndex]
      names(colvalues) <- rownames(tempMat)
      colNames <- colnames(tempMat)[colIndex]
      colvalues <- colvalues[which(colvalues < (-1) * logRatio)]
      filteredNames <- names(colvalues)
      convertedGeneNames <- probesetGeneNameMappings[filteredNames]
      res <- cbind(colvalues, convertedGeneNames)
      colnames(res) <- c("M Value", "Gene Name")
      rownames(res) <- filteredNames
      # write.csv(res, file = paste("output/differential/Repression_", colNames, colnamelist[i], ".csv"))
    }, cutOffRepression)
    
    # names(difflist)<-colnames(tempMat) Induced
    newRownames <- c()
    newMatInduced <- do.call(rbind, lapply(1:nrow(tempMat), function(rowindex, logRatio, totalSamples, majority) {
      # get row vector
      rowvector <- round(as.numeric(tempMat[rowindex, ]), 6)
      probesetName <- rownames(tempMat)[rowindex]
      rowvector2 <- rowvector[which(rowvector >= logRatio)]
      
      numOfMatchCriteria <- length(rowvector2)
      
      if (numOfMatchCriteria >= majority) {
        return(c(probesetName, tempMat[rowindex, ], majority = majority))
      }else{
        return(NULL)
      }
    }, cutOffInduction, totalSamples, majority))
    
    # rownames(newMatInduced)<-newRownames
    
    # repressed
    newRownames <- c()
    newMatRepressed <- do.call(rbind, lapply(1:nrow(tempMat), function(rowindex, logRatio, totalSamples, majority) {
      # get row vector
      rowvector <- round(as.numeric(tempMat[rowindex, ]), 6)
      probesetName <- rownames(tempMat)[rowindex]
      rowvector2 <- rowvector[which(rowvector <= (-1 * logRatio))]
      numOfMatchCriteria <- length(rowvector2)
      
      if (numOfMatchCriteria >= majority) {
        return(c(probesetName, tempMat[rowindex, ], majority = majority))
      } else {
        return(NULL)
      }
    }, cutOffRepression, totalSamples, majority))
    # rownames(newMatRepressed)<-newRownames
    
    diffindex <- length(diffExpressed) + 1
    diffExpressed[[diffindex]] <- list()
    diffExpressed[[diffindex]][[1]] <- newMatInduced[, -1]
    diffExpressed[[diffindex]][[2]] <- newMatRepressed[, -1]
    diffExpressed[[diffindex]][[3]] <- as.vector(newMatInduced[, 1])
    diffExpressed[[diffindex]][[4]] <- as.vector(newMatRepressed[, 1])
    diffExpressed[[diffindex]][[5]] <- probesetGeneNameMappings[newMatInduced[, 1]]
    diffExpressed[[diffindex]][[6]] <- probesetGeneNameMappings[newMatRepressed[, 1]]
    
    names(diffExpressed[[diffindex]])[[1]] <- "Induced_M_Value"
    names(diffExpressed[[diffindex]])[[2]] <- "Repressed_M_Value"
    names(diffExpressed[[diffindex]])[[3]] <- "Induced_ProbeID"
    names(diffExpressed[[diffindex]])[[4]] <- "Repressed_ProbeID"
    names(diffExpressed[[diffindex]])[[5]] <- "Induced_Genes"
    names(diffExpressed[[diffindex]])[[6]] <- "Repressed_Genes"
    
    names(diffExpressed)[[diffindex]] <- colnamelist[i]
  }
  
  # find all common genes in the result, the result should be
  commonNames <- unique(unlist(lapply(diffExpressed, function(subdata) c(unlist(subdata[["Induced_ProbeID"]]), unlist(subdata[["Repressed_ProbeID"]])))))
  
  filteredData <- list()
  filteredData[[1]] <- diffExpressed
  filteredData[[2]] <- probesetGeneNameMappings[commonNames]
  names(filteredData)[[1]] <- "DifferentialExpression"
  names(filteredData)[[2]] <- "DifferentialMappings"
  futile.logger::flog.info("Leave identifyDifferentiallyExpressedGenes zone.")
  filteredData
}


#' getSpecificExpressedGenes
#' 
#' @param orderSampleTimeSeries the original timeseries data that contain continual values
#' @param genelist the target genes
#' @export
getSpecificExpressedGenes <- function(orderSampleTimeSeries, genelist = c()) {
  futile.logger::flog.info(sprintf("Enter getSpecificExpressedGenes zone: orderSampleTimeSeries=%s, genelist=%s",
                                   head(orderSampleTimeSeries),
                                   head(genelist)))
  alldifferExpressednames <- unique(rownames(orderSampleTimeSeries[[1]]))  #find the maximum common gene set
  probesetMapping <- mapToGeneNameWithhgu133plus2Db(alldifferExpressednames)
  filteredDataKnown <- lapply(orderSampleTimeSeries, function(subdata) subdata[probesetMapping$MappedProbeset[, 1], ])
  # filteredDataUnKnown <- lapply(orderSampleTimeSeries, function(subdata) subdata[probesetMapping$UnMappedProbeset[, 1], ])
  filteredDataUnKnown <- list()
  # update known probset names to gene names
  filteredDataKnown <- lapply(filteredDataKnown, function(mtx) {
    geneNames <- probesetMapping$MappedProbeset[which(probesetMapping$MappedProbeset$PROBEID %in% rownames(mtx)), ]$SYMBOL
    rownames(mtx) <- geneNames
    return(mtx)
  })
  
  if (length(genelist) > 0) {
    mappedtargetGenes <- probesetMapping$MappedProbeset[which(probesetMapping$MappedProbeset$PROBEID %in% genelist), ]$SYMBOL
    unmappedtargetGenes <- probesetMapping$UnMappedProbeset[which(probesetMapping$UnMappedProbeset$PROBEID %in% genelist), ]$PROBEID
    genelist <- unique(c(mappedtargetGenes, unmappedtargetGenes))
  }
  commonNames <- c()
  # find all common genes in the result, the result should be
  for (i in seq_along(filteredDataKnown)) {
    subset <- filteredDataKnown[[i]]
    if (i > 1) {
      subset <- subset[rownames(subset) %in% commonNames, ]
      commonNames <- rownames(subset)
    } else {
      commonNames <- rownames(subset)
    }
  }
  # allCommonGeneNames<-unique(unlist(sapply(filteredDataKnown,function(x)rownames(x)))) #find the maximum common gene set
  filteredDataKnown <- lapply(filteredDataKnown, function(subdata) subdata[commonNames, ])
  if (length(genelist) > 0) {
    filteredDataKnown <- lapply(filteredDataKnown, function(subdata) subdata[rownames(subdata) %in% genelist, ])
    filteredDataUnKnown <- lapply(filteredDataUnKnown, function(subdata) subdata[rownames(subdata) %in% genelist, ])
  }
  
  sampleNames <- names(filteredDataKnown)
  allDiffData <- lapply(sampleNames, function(sampleName) rbind(filteredDataKnown[[sampleName]], filteredDataUnKnown[[sampleName]]))
  names(allDiffData) <- sampleNames
  filteredData <- list()
  filteredData[[1]] <- filteredDataKnown
  filteredData[[2]] <- filteredDataUnKnown
  filteredData[[3]] <- allDiffData
  filteredData[[4]] <- probesetMapping
  names(filteredData)[[1]] <- "KnownDiffExpressedTimeSeriesData"
  names(filteredData)[[2]] <- "UnKnownDiffExpressedTimeSeriesData"
  names(filteredData)[[3]] <- "AllDiffExpressedTimeSeriesData"
  names(filteredData)[[4]] <- "ProbesetMapping"
  futile.logger::flog.info("Leave getSpecificExpressedGenes zone.")
  filteredData
}

#'A method that map the row names (probesets) of a specific timeseries data into gene names
#'Old function name mapPropersetsWithhgu133plus2Db
#'@param orderSampleTimeSeries A sorted time series data, which is the output of the method reorderSampleTimeSeries
#' @export
convertTimeseriesProbsetNameToGeneName <- function(orderSampleTimeSeries) {
  futile.logger::flog.info(sprintf("Enter convertTimeseriesProbsetNameToGeneName zone: orderSampleTimeSeries=%s",
                                   head(orderSampleTimeSeries)))
  alldifferExpressednames <- c()
  for (i in seq_along(orderSampleTimeSeries)) {
    row_names <- rownames(orderSampleTimeSeries[[i]])
    alldifferExpressednames <- c(alldifferExpressednames, row_names)
    alldifferExpressednames <- unique(alldifferExpressednames)
  }
  probesetMapping <- mapToGeneNameWithhgu133plus2Db(alldifferExpressednames)
  mappingsheet <- probesetMapping$DeDuplicatedMapping
  
  mapto <- lapply(alldifferExpressednames, function(name) as.character(mappingsheet[which(mappingsheet$original == name), ][, 2]))
  cond <- sapply(mapto, function(x) length(x) > 0)
  mapto <- mapto[cond]
  mapto <- unique(unlist(mapto))
  
  convert_data <- lapply(orderSampleTimeSeries, function(subdata, mappingsheet, probesetMapping) {
    
    subdata <- subdata[which(rownames(subdata) %in% mappingsheet$original), ]
    gene_names <- mappingsheet[mappingsheet$original %in% rownames(subdata), "mapTo"]
    gene_names <- mapProbesetNames(rownames(subdata), probesetMapping)
    col_names <- colnames(subdata)
    res <- cbind(subdata, gene_names)
    colnames(res)[length(col_names) + 1] <- "geneName"
    #remove duplicate by mean
    res <- as.data.frame(res) 
    res[,(length(col_names) + 1)] <- as.character(res[,(length(col_names) + 1)])
    for(i in seq_len(length(col_names))) {
      res[, i] <- as.numeric(res[, i])
    }
    res <- dplyr::group_by(res, geneName)
    res <- dplyr::summarise_all(res, mean) 
    newrowNames <-dplyr::pull(res, geneName)
    res <- as.matrix(res[,-1])
    rownames(res) <- newrowNames
    res
    ##subdata[!duplicated(rownames(subdata)), ]
  }, mappingsheet, probesetMapping)

  filteredData <- list()
  filteredData[[1]] <- convert_data
  filteredData[[2]] <- probesetMapping
  filteredData[[3]] <- mappingsheet
  filteredData[[4]] <- alldifferExpressednames
  
  names(filteredData)[[1]] <- "convert_data"
  names(filteredData)[[2]] <- "ProbesetMapping"
  names(filteredData)[[3]] <- "mappingsheet"
  names(filteredData)[[4]] <- "alldifferExpressednames"
  futile.logger::flog.info("Leave convertTimeseriesProbsetNameToGeneName zone.")
  filteredData
}

#' Annotation method
#' 
#' The method maps the probesets with gene names
#' 
#' @param probesets the target probeset names
#' 
#' @export
mapToGeneNameWithhgu133plus2Db <- function(probesets) {
  futile.logger::flog.info(sprintf("Enter mapToGeneNameWithhgu133plus2Db zone: probesets=%s",
                                   head(probesets)))
  
  # if the following package hasen't been installed, install them source('http://bioconductor.org/biocLite.R') biocLite('annotate') biocLite('hgu133plus2.db')
  # require("AnnotationDbi")
  # require("hgu133plus2.db")
  # find all types keyTypes <- keytypes(hgu133plus2.db)
  res <- list()
  ## 70542
  getGeneNames <- AnnotationDbi::select(hgu133plus2.db::hgu133plus2.db, 
                                        keys = probesets, 
                                        columns = c("UNIGENE", "SYMBOL", "ENTREZID", "GENENAME"), 
                                        keytype = "PROBEID")
  ## remove RNA type of non genes' probeset
  ## could be none
  ## 59211 after removed NA 
  getGeneNames <- getGeneNames[!is.na(getGeneNames$SYMBOL) & 
                                 !is.na(getGeneNames$ENTREZID) & 
                                 !is.na(getGeneNames$GENENAME), ] 

  ## remove duplicate records based on probeid 44109/59211, duplicates 15867
  mapped <- getGeneNames[!duplicated(getGeneNames[, "PROBEID"]), ]  #de-duplicated by PROBEID
  mapped <- mapped [order(mapped$ENTREZID),]
  ## remove duplicated records based on entrezid, 20764
  deduplicatedmapped <- mapped[!duplicated(mapped[, "ENTREZID"]), ]  #de-duplicated by ENTREZID
  
  mappingCol1 <- c()
  mappingCol2 <- c()
  num_rows <- nrow(deduplicatedmapped)
  mappings <- lapply(seq_len(num_rows), function(i, mapped, deduplicatedmapped) {
    entrezId <- deduplicatedmapped[i, "ENTREZID"]
    probesetIds <- mapped[mapped$ENTREZID == entrezId, ]$PROBEID
    matchedduplicatedId <- rep(deduplicatedmapped[i, "PROBEID"], length(probesetIds))
    cbind(probesetIds, matchedduplicatedId)
  }, mapped, deduplicatedmapped)
  mappingsheet <- do.call(rbind, mappings)
  colnames(mappingsheet) <- c("original", "mapTo")
  
  res[["MappedProbeset"]] <- mapped
  res[["DeDuplicatedMappedProbeset"]] <- deduplicatedmapped
  res[["DeDuplicatedMapping"]] <- as.data.frame(mappingsheet)
  class(res) <- "AnnotationMappedObject"
  futile.logger::flog.info("Leave mapToGeneNameWithhgu133plus2Db zone.")
  res
}

#' Innternal method
#' @noRd
convert_probeset_to_gene <- function(annotationMappedObject, probeset) {
  if (!inherits(annotationMappedObject, "AnnotationMappedObject")) 
    stop("the value of the parameter: annotationMappedObject must be inherited from the class of AnnotationMappedObject")
  mapping <- as.data.frame(annotationMappedObject$DeDuplicatedMapping)
  get_map_to <- mapping[which(mapping$orginal %in% c(probeset)), ]$mapto
  annotationMappedObject$DeDuplicatedMappedProbeset[which(annotationMappedObject$DeDuplicatedMappedProbeset$PROBEID %in% c(get_map_to)), ]
}


#'Retrieving genome sequence data using SeqinR. The code is a sample from the book 'a little book of r for bioinformatic'
#'
#'@param accession NBCI accession number
#'@return genome sequence
#'@export
getncbiseq <- function(accession) {
  # first find which ACNUC database the accession is stored in:
  dbs <- c("genbank", "refseq", "refseqViruses", "bacterial")
  numdbs <- length(dbs)
  for (i in 1:numdbs) {
    db <- dbs[i]
    seqinr::choosebank(db)
    # check if the sequence is in ACNUC database ?\200\231db?\200\231:
    resquery <- try(query(".tmpquery", paste("AC=", accession)), silent = TRUE)
    if (!(inherits(resquery, "try-error"))) {
      queryname <- "query2"
      thequery <- paste("AC=", accession, sep = "")
      seqinr::query("queryname", "thequery")
      # see if a sequence was retrieved:
      seq <- seqinr::getSequence(seqinr::query2$req[[1]])
      seqinr::closebank()
      return(seq)
    }
    seqinr::closebank()
  }
  print(paste("ERROR: accession", accession, "was not found"))
}

#'@export
outputFasta <- function(names, sequenceData) {
  newfile <- paste(names, ".fasta", sep = "")
  print(paste("Output sequence to", newfile))
  seqinr::write.fasta(names = names, sequences = seqinr::dengueseq, file.out = newfile)
}

#'@export
readFasta <- function(fastaFileName) {
  print(paste("read sequence from", fastaFileName))
  seqinr::read.fasta(file = fastaFileName)
}

#'@export
slidingwindowplot <- function(windowsize, inputseq) {
  starts <- seq(1, length(inputseq) - windowsize, by = windowsize)
  n <- length(starts)  # Find the length of the vector 'starts'
  chunkGCs <- numeric(n)  # Make a vector of the same length as vector 'starts', but just containing
  for (i in 1:n) {
    chunk <- inputseq[starts[i]:(starts[i] + windowsize - 1)]
    chunkGC <- seqinr::GC(chunk)
    print(chunkGC)
    chunkGCs[i] <- chunkGC
  }
  plot(starts, chunkGCs, type = "b", xlab = "Nucleotide start position", ylab = "GC content")
}

########## Poly fix polyfit(y,x,maxDegree) fits all polynomials up to degree max degree; y is vector for response varaible, x for predictor; creates an object of
########## class 'PolyRegression', consisting of outputs from the various regression models, plus the orginal data
polyfit <- function(y, x, maxDegree) {
  pwrs <- powers(x, maxDegree)  #form powers or predictor variable
  lmout <- list()  #start to build class
  class(lmout) <- "PolyRegression"  # create a new class
  for (i in 1:maxDegree) {
    lmo <- lm(y ~ pwrs[, 1:i])
    # extend the lm class here, with the cross-validated predictions
    lmo$fitted.xvvalues <- lvoneout(y, pwrs[, 1:i, drop = F])
    lmout[[i]] <- lmo
  }
  lmout$x <- x
  lmout$y <- y
  lmout
}

# generic print() for an object fits of class 'polyreg': print cross-validated mean-squarted prediction errors
print.PolyRegression <- function(fits) {
  maxdeg <- length(fits) - 2  #count lm() outputs only, not $x and $y
  n <- length(fits$y)
  tbl <- matrix(nrow = maxdeg, ncol = 1)
  cat("mean squared prediction errors, by degree \n")
  colnames(tbl) <- "MSPE"
  for (i in 1:maxdeg) {
    fi <- fits[[i]]
    errs <- fits$y - fi$fitted.xvvalues
    spe <- sum(errs^2)
    tbl[i, 1] <- spe/n
  }
  print(tbl)
}

# generic plot(); plots fits against raw data
plot.PolyRegression <- function(fits) {
  plot(fits$x, fits$y, xlab = "X", ylab = "Y")  #plot data points as background
  maxdg <- length(fits) - 2
  cols <- c("red", "green", "blue")
  dg <- curvecount <- 1
  while (dg <- maxdg) {
    prompt <- paste("RETURN for XV fit for degree", dg, " or type degree", "or q for quit")
    rl <- readline(prompt)
    dg <- if (rl == "") 
      dg else if (rl != "q") 
        as.integer(rl) else break
    lines(fits$x, fits[[dg]]$fitted.values, col = cols[curvecount%%3 + 1])
    dg <- dg + 1
    curvecount <- curvecount + 1
  }
}

# forms matrix of powers of the vector x, through degree dg
powers <- function(x, dg) {
  pw <- matrix(x, nrow = length(x))  #transpose X
  prod <- x
  for (i in 2:dg) {
    prod <- prod * x
    pw <- cbind(pw, prod)
  }
  pw
}

#' finds cross-validated predicted values; could be made much fasgter via matrix-update methods
#' @noRd
lvoneout <- function(y, xmat) {
  n <- length(y)
  predy <- vector(length = n)
  for (i in 1:n) {
    # repress, leavingout ith observation
    lmo <- lm(y[-i] ~ xmat[-i, ])
    betahat <- as.vector(lmo$coef)
    # the 1 accommodates the constant term
    predy[i] <- betahat %*% c(1, xmat[i, ])
  }
  predy
}

#' polynomial function of x, coefficients cfs
#'
#' @noRd
calculate_polynomial <- function(x, cfs) {
  val <- cfs[1]
  prod <- 1
  dg <- length(cfs) - 1
  for (i in 1:dg) {
    prod <- prod * x
    val <- val + cfs[i + 1] * prod
  }
}

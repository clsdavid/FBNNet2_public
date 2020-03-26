
#' get Statistic Measures
#' 
#' A method to calculate the statistic measurements based on feature names
#'
#' @param featurenames a list of feature names
#' @param timeseriesdata the original timeseries data that contains continual values
#' @return statstic measures
#' 
getStatisticMeasures <- function(featurenames, timeseriesdata) {
  
  # generate data
  
  res <- lapply(featurenames, function(feature) {
    sampledata <- unlist(lapply(timeseriesdata, function(sample) sample[feature, ]))
    std_sampledata <- stats::sd(sampledata)  #standard deviation
    mean_sampledata <- mean(sampledata)
    length_sampledata <- length(sampledata)
    sem_sampledata <- std_sampledata/sqrt(length(sampledata))  #Standard error of the mean
    max_sampledata <- max(sampledata)
    min_sampledata <- min(sampledata)
    median_sampledata <- stats::median(sampledata)
    
    result <- c(id = feature, std = std_sampledata, mean = mean_sampledata, length = length_sampledata, sem = sem_sampledata, max = max_sampledata, 
      min = min_sampledata, median = median_sampledata)
    return(result)
  })
  
  measures <- do.call(rbind, res)
  rownames(measures) <- measures[, 1]
  measures <- measures[order(rownames(measures)), ]
  measures
}

#' This method is used to converts normalized
#' timeseries data into a list of samples
#' @param normalizedData An output of the method 
#' normalizeTimesereisRawData that normalized timeseries data
#' @param func A function that specified how to split the column names
#' @param splitor Seperator.
#'
#'@export
convertIntoSampleTimeSeries <- function(normalizedData, func = function(x) paste(x[1], x[2], x[3], sep = "-"), splitor = as.character("-")) {
  
  mat <- normalizedData
  # Obtain the last part of each column names
  groups <- sapply(strsplit(x = colnames(mat), split = splitor), func)
  print(groups)
  # Go through each unique column name and extract the corresponding columns
  res <- lapply(unique(groups), function(x) mat[, which(groups == x)])
  names(res) <- unique(groups)
  
  # need to convert column name into numbers str_extract_all(names(df), '[0-9]+')
  res
}

#'This method is used to sort time series order based on columns
#'@param convertedTimeSeries An output of the method 
#'convertIntoSampleTimeSeries that converts normalized timeseries data
#' into a list of samples
#'@param func A function that specified how to split the column names
#'@param splitor Separator.
#'@export
reorderSampleTimeSeries <- function(convertedTimeSeries, func = function(x) x[length(x)], splitor = as.character("-")) {
  len <- length(convertedTimeSeries)
  res <- lapply(seq_len(len), function(index) {
    result <- convertedTimeSeries[[index]]
    ocolnames <- colnames(result)
    ocolnames <- gsub("\\..*", "", ocolnames)
    ocolnames <- sapply(strsplit(x = ocolnames, split = splitor), func)
    ocolnames <- as.integer(stringr::str_extract_all(ocolnames, "[0-9]+"))
    colnames(result) <- ocolnames
    result2 <- result[, order(as.numeric(colnames(result)))]
    result2 <- result2[order(rownames(result2)), ]
    return(result2)
  })
  names(res) <- names(convertedTimeSeries)
  res
}


#' An function to divide large data into small groups
#' @param  discretedTimeSeriesdata discreted timeseries data
#' @param maxElements the max elements to divide the 
#' discreted timeseries data into
#' @param maxK The maximum level that can be drilled into
dividedDataIntoSubgroups <- function(discretedTimeSeriesdata, maxElements = 20, maxK = 4) {
  futile.logger::flog.info(sprintf("Enter dividedDataIntoSubgroups zone:
                   maxElements=%s", maxElements))
  genes <- unique(unlist(lapply(discretedTimeSeriesdata, rownames)))
  
  len_genes <- length(genes)
  n_subgroup <- ceiling(len_genes/maxElements)
  
  sub_genes <- split(genes, rep(1:n_subgroup, each = maxElements)[1:len_genes])
  sub_matrix <- lapply(sub_genes, function(sub, matx_data) {
    lapply(matx_data, function(matx, sub) {
      matx[rownames(matx) %in% sub, ]
    }, sub)
  }, discretedTimeSeriesdata)
  
  if (maxK > length(sub_genes)) {
    maxK <- length(sub_genes)
  }
  
  combined_groups <- list()
  processed <- c()
  for (i in seq_along(sub_genes)) {
    processed <- c(processed, i)
    sub_genes2 <- sub_genes[-processed]
    for (j in seq_along(sub_genes2)) {
      if (identical(sub_genes[[i]], sub_genes2[[j]])) 
        (next)()
      
      newset <- unique(c(sub_genes[[i]], sub_genes2[[j]]))
      newset <- genes[which(genes %in% newset)]
      combined_groups[[length(combined_groups) + 1]] <- newset
    }
  }
  res_combined_groups <- list()
  if (maxK > 1) {
    for (m in seq_len(maxK - 1)) {
      for (i in seq_along(combined_groups)) {
        group <- combined_groups[[i]]
        conds <- vapply(sub_genes, function(sub) !all(sub %in% group), logical(1))
        filered <- sub_genes[conds]
        for (j in seq_along(filered)) {
          newset <- unique(c(group, filered[[j]]))
          newset <- genes[genes %in% newset]
          conds2 <- sapply(res_combined_groups, function(sub, check) {
          all(check %in% sub)
          }, newset)
          if (!any(conds2)) {
          res_combined_groups[[length(res_combined_groups) + 1]] <- newset
          futile.logger::flog.info(sprintf("dividedDataIntoSubgroups zone: group=%s", paste(newset, sep = ",", collapse = ",")))
          }
        }
      }
    }
  } else {
    res_combined_groups <- combined_groups
  }
  
  futile.logger::flog.info(sprintf("(dividedDataIntoSubgroups) generates: %s combined groups", length(res_combined_groups)))
  
  res <- list(clusters = sub_genes, clusters_timeseries = sub_matrix, combined_clusters = res_combined_groups, discretedTimeSeriesdata = discretedTimeSeriesdata, 
    target_genes = genes)
  class(res) <- "ClusteredTimeseriesData"
  res
}

#' An internal function that divides large data into small groups.
#' @param discretedTimeSeriesdata discreted timeseries data
#' @param maxElements the max elements to divide the discreted timeseries data into
#' @param maxK The max levels.
#' @return A list of data
#' @export
dividedDataIntoSmallgroups <- function(discretedTimeSeriesdata, maxElements = 20, maxK = 4) {
  futile.logger::flog.info(sprintf("Enter dividedDataIntoSmallgroups zone: maxElements=%s", maxElements))
  genes <- unique(unlist(lapply(discretedTimeSeriesdata, rownames)))
  
  len_genes <- length(genes)
  n_subgroup <- ceiling(len_genes/maxElements)
  
  sub_genes <- split(genes, rep(1:n_subgroup, each = maxElements)[1:len_genes])
  
  futile.logger::flog.info(sprintf("(dividedDataIntoSmallgroups) generates: %sgroups", length(sub_genes)))
  
  res <- list(clusters = sub_genes, discretedTimeSeriesdata = discretedTimeSeriesdata, conditional_genes = genes)
  class(res) <- "ClusteredTimeseriesData"
  res
}

#'A methiod to reduce the timeseries data based on the gene list
#'@param timeseries a list of samples that contains genes on 
#'row and time steps on columns
#'@param genelist An vector of genes
#' @export
getRelatedGeneTimeseries <- function(timeseries, genelist = c()) {
  lapply(timeseries, function(sheet) sheet[rownames(sheet) %in% genelist, ])
}

#'generate all combaination binary data based on an vector of genes
#'@param genelist An vector of genes
#'@param begin return begin index
#'@param last  return last index
#'@export
generateAllCombinationBinary <- function(genelist = c(), begin = 1, last = 0) {
  if (length(genelist) == 0) {
    stop("The genelist is empty")
  }
  
  if (last == 0) {
    last = 2^(length(genelist))
  }
  res <- list()
  for (gene in genelist) {
    newindex <- length(res) + 1
    res[[newindex]] <- c(0, 1)
    names(res)[[newindex]] <- gene
  }
  result <- expand.grid(res)[begin:last, ]
  result <- t(result)
  ncols <- ncol(result)
  split(result, rep(seq_len(ncols), each = nrow(result)))
}

#'A method to generate binary data randomly
#'@param genelist An vector of genes
#'@param maxState that should be less tha 2^length(genelist)
#'@examples
#'individualgenes<-c('a','b','c')
#'testdata2<-randomGenerateBinary(individualgenes,maxState = 4)
#' @export
randomGenerateBinary <- function(genelist = c(), maxState = 0) {
  if (length(genelist) == 0) {
    stop("The genelist is empty")
  }
  
  if (maxState == 0) {
    maxState <- 2^length(genelist)
  }
  
  result <- list()
  index <- 1
  while (index <= maxState) {
    res <- c()
    for (gene in genelist) {
      res <- c(res, as.numeric(randomSelection(0.5)))
    }
    index <- index + 1
    names(res) <- genelist
    result[[length(result) + 1]] <- res
  }
  
  result
}

#'A method to generate BoolNet type of timeseries data
#'
#' @param network A traditional Boolea type of network that can be used for BoolNet
#' @param initialStates A list of initial states
#' @param numMeasurements the number of timepoints that need to reconstruct
#' @param type the type of the network in traditional Boolean modelling
#' @param geneProbabilities optional if type is probabilistic, and then 
#' it is a need to specify the probabilities for each gene
#' @export
genereateBoolNetTimeseries <- function(network, initialStates, numMeasurements, type = c("synchronous", "asynchronous", "probabilistic"), geneProbabilities = NULL) {
  # initialStates<-initialStates[match(network$genes,rownames(initialStates)),]
  lapply(initialStates, function(state) {
    res <- state
    startState <- state
    for (j in 2:numMeasurements) {
      startState <- BoolNet::stateTransition(network, startState, type = type, geneProbabilities = geneProbabilities)
      res <- cbind(res, startState)
    }
    rownames(res) <- network$genes
    colnames(res) <- seq_len(ncol(res))
    return(res)
  })
}

#'This method compare two timeseries data and generate similar report
#'@param timeseriesdata1 the source time series data
#'@param timeseriesdata2 tge target time series data
#' @export
generateSimilaryReport <- function(timeseriesdata1, timeseriesdata2) {
  # validate network types
  similar <- checkSimilarity(timeseriesdata1, timeseriesdata2)
  
  cond1 <- vapply(similar, function(entry) entry[[1]] == "similar", logical(1))
  cond2 <- vapply(similar, function(entry) entry[[1]] == "likely", logical(1))
  cond3 <- vapply(similar, function(entry) entry[[1]] == "verysimilar", logical(1))
  cond4 <- vapply(similar, function(entry) entry[[1]] == "unlikely", logical(1))
  cond5 <- vapply(similar, function(entry) entry[[1]] == "veryunlikely", logical(1))
  
  res <- list()
  res[[1]] <- similar
  
  
  res[[2]] <- similar[cond1]
  res[[3]] <- similar[cond2]
  res[[4]] <- similar[cond3]
  res[[5]] <- similar[cond4]
  res[[6]] <- similar[cond5]
  
  names(res)[[1]] <- "SimilarityReport"
  names(res)[[2]] <- "Similar"
  names(res)[[3]] <- "Likely"
  names(res)[[4]] <- "verysimilar"
  names(res)[[5]] <- "unlikely"
  names(res)[[6]] <- "veryunlikely"
  
  # benchmark result
  pmcond <- vapply(similar, function(entry) as.numeric(entry[[2]]) == 1, logical(1))
  mmcond <- vapply(similar, function(entry) as.numeric(entry[[2]]) < 1, logical(1))
  
  pm <- similar[pmcond]
  mm <- similar[mmcond]
  
  mmvalues <- unlist(lapply(mm, function(cond) as.numeric(cond[[2]])))
  avgMMvalues <- ifelse(length(mm) > 0, sum(mmvalues)/length(mm), 0)
  
  errorRate <- ((1 - avgMMvalues) * length(mm))/(length(pm) + length(mm))
  perfectMatchedRate <- length(pm)/(length(pm) + length(mm))
  missMatchedRate <- length(mm)/(length(pm) + length(mm))
  accurateRate <- (length(pm) + avgMMvalues * length(mm))/(length(pm) + length(mm))
  
  res[[7]] <- errorRate
  res[[8]] <- accurateRate
  res[[9]] <- missMatchedRate
  res[[10]] <- perfectMatchedRate
  names(res)[[7]] <- "ErrorRate"
  names(res)[[8]] <- "AccurateRate"
  names(res)[[9]] <- "MissMatchedRate"
  names(res)[[10]] <- "PerfectMatchedRate"
  res
}

#' A method to duplicate a row
#' @param x the row value that need to be repeated
#' @param n the number of repeats
#' @export
rep_row <- function(x, n) {
  matrix(rep(x, each = n), nrow = n)
}

#' A method to duplicate a column
#' @param x the column value that need to be repeated
#' @param n the number of repeats
#' @export
rep_col <- function(x, n) {
  matrix(rep(x, each = n), ncol = n, byrow = TRUE)
}

#' A method to bind two matrix by row
#' @param x matrix one
#' @param y matrix two
#' @export
rbind_all_columns <- function(x, y) {
  
  x.diff <- setdiff(colnames(x), colnames(y))
  y.diff <- setdiff(colnames(y), colnames(x))
  
  x[, c(as.character(y.diff))] <- NA
  
  y[, c(as.character(x.diff))] <- NA
  
  return(rbind(x, y))
}

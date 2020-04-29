
#' An internal function that divides large data into small groups.
#' @param vector A vector that need to be splited
#' @param maxElements the max elements to divide the discreted timeseries data into
#' @return A list of data
#' @export
dividedVectorIntoSmallgroups <- function(vector, maxElements = 20) {
  if (!is.vector(vector)) {
    stop("The parameter 'vector' must be a vector")
  }
  futile.logger::flog.info(sprintf("Enter dividedVectorIntoSmallgroups zone: maxElements=%s", maxElements))

  len_vector <- length(vector)
  n_subgroup <- ceiling(len_vector/maxElements)
  
  sub_vector <- split(vector, rep(1:n_subgroup, each = maxElements)[1:len_vector])
  
  futile.logger::flog.info(sprintf("(dividedVectorIntoSmallgroups) generates: %sgroups", 
                                   length(sub_vector)))
  
  list(clusters = sub_vector, original = vector)
}

#'A function to reduce the timeseries data based on the gene list
#'@param timeseries a list of samples that contains genes on 
#'row and time steps on columns
#'@param genelist An vector of genes
#' @export
getRelatedGeneTimeseries <- function(timeseries, genelist = c()) {
  lapply(timeseries, function(sheet) sheet[rownames(sheet) %in% genelist, ])
}

#' Generate all combination of binary data based on an vector of genes
#' @param genelist An vector of genes
#' @param begin The begin index
#' @param last  The last index, the default is 0 means 2^length(genelist)
#' @examples
#' individualgenes<-c('a','b','c')
#' testdata2 <- generateAllCombinationBinary(individualgenes)
#' testdata2
#'@export
generateAllCombinationBinary <- function(genelist = c(), begin = 1, last = 0) {
  if (length(genelist) == 0) {
    stop("The genelist is empty")
  }
  
  if (last == 0) {
    last <- 2^(length(genelist))
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
#'@param maxState that should be less that 2^length(genelist)
#'@examples
#'individualgenes<-c('a','b','c')
#'testdata2 <- randomGenerateBinary(individualgenes,maxState = 4)
#'testdata2
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

#'A method to generate BoolNet type of timeseries data for generating BoolNet type of 
#'testing timeseries data
#'
#' @param network A traditional Boolean type of network that can be used for BoolNet
#' @param initialStates A list of initial states
#' @param numMeasurements the number of timepoints that need to reconstruct
#' @param type the type of the network in traditional Boolean modelling
#' @param geneProbabilities optional if type is probabilistic, and then 
#' it is a need to specify the probabilities for each gene
#' @export
genereateBoolNetTimeseries <- function(network, 
                                       initialStates, 
                                       numMeasurements, 
                                       type = c("synchronous", "asynchronous", "probabilistic"), 
                                       geneProbabilities = NULL) {
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
#'@param timeseriesdata1 The source time series data
#'@param timeseriesdata2 The target time series data
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

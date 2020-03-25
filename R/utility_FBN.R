#'Remove duplidate timeseries sample data
#'
#'@param timeseriescube A list of time series data
#'@return The return a list of time series data
#'@examples
#' mat1<-matrix(c('1','2','3','4','5','6','7','8','9'),3,3)
#' mat2<-matrix(c('1','2','3','4','5','6','7','8','9'),3,3)
#' listtest<-list(mat1,mat2)
#' FBNDataReduction(listtest)
#'@export
FBNDataReduction <- function(timeseriescube) {
    duplicateIndexes <- duplicated(timeseriescube)
    timeseriescube[!duplicateIndexes]
}


#'Check the similarity between time series
#'
#'@param originalTimeseriesCube The original data set that
#' contains samples and each sample contains genes and time points
#'@param reconstructedTimeSeriesCube The reconstructed data set
#' that contains samples and each sample contains genes and time points
#'@return similarity report
#'@examples
#' ##coming later
#'@export
checkSimilarity <- function(originalTimeseriesCube, reconstructedTimeSeriesCube) {
    res <- list()
    if (!identical(length(originalTimeseriesCube), length(reconstructedTimeSeriesCube))) {
        stop("The length of each timeseries data must be identical")
    }
    
    for (i in seq_along(originalTimeseriesCube)) {
        if (!identical(dim(originalTimeseriesCube[[i]]), dim(reconstructedTimeSeriesCube[[i]]))) {
            stop("The dimension of each timeseries data must be identical")
        }
        res[[i]] <- similarityBetweenMatrix(originalTimeseriesCube[[i]],
                                            reconstructedTimeSeriesCube[[i]], i)
    }
    res
}


#'Generate similarity report
#'
#'@param similarityreport The raw similarity report which was created by the function checkSimilarity
#'@return An organized similarity report
#'@examples
#' ##coming later
#'@export
generateSimilarReport <- function(similarityreport) {
    cond1 <- sapply(similarityreport, function(entry) as.numeric(entry[[2]]) >= 0.9)
    cond2 <- sapply(similarityreport, function(entry) as.numeric(entry[[2]]) >= 0.8 & as.numeric(entry[[2]]) < 0.9)
    cond3 <- sapply(similarityreport, function(entry) as.numeric(entry[[2]]) >= 0.7 & as.numeric(entry[[2]]) < 0.8)
    cond4 <- sapply(similarityreport, function(entry) as.numeric(entry[[2]]) >= 0.6 & as.numeric(entry[[2]]) < 0.7)
    cond5 <- sapply(similarityreport, function(entry) as.numeric(entry[[2]]) >= 0.5 & as.numeric(entry[[2]]) < 0.6)
    cond6 <- sapply(similarityreport, function(entry) as.numeric(entry[[2]]) >= 0.4 & as.numeric(entry[[2]]) < 0.5)
    cond7 <- sapply(similarityreport, function(entry) as.numeric(entry[[2]]) >= 0.3 & as.numeric(entry[[2]]) < 0.4)
    cond8 <- sapply(similarityreport, function(entry) as.numeric(entry[[2]]) >= 0.2 & as.numeric(entry[[2]]) < 0.3)
    cond9 <- sapply(similarityreport, function(entry) as.numeric(entry[[2]]) >= 0.1 & as.numeric(entry[[2]]) < 0.2)
    cond10 <- sapply(similarityreport, function(entry) as.numeric(entry[[2]]) < 0.1)
    
    res <- list()
    res[[1]] <- similarityreport[cond1]
    res[[2]] <- similarityreport[cond2]
    res[[3]] <- similarityreport[cond3]
    res[[4]] <- similarityreport[cond4]
    res[[5]] <- similarityreport[cond5]
    res[[6]] <- similarityreport[cond6]
    res[[7]] <- similarityreport[cond7]
    res[[8]] <- similarityreport[cond8]
    res[[9]] <- similarityreport[cond9]
    res[[10]] <- similarityreport[cond10]
    names(res)[[1]] <- "A"
    names(res)[[2]] <- "B"
    names(res)[[3]] <- "C"
    names(res)[[4]] <- "D"
    names(res)[[5]] <- "E"
    names(res)[[6]] <- "F"
    names(res)[[7]] <- "G"
    names(res)[[8]] <- "H"
    names(res)[[9]] <- "I"
    names(res)[[10]] <- "J"
    res
}

#' A function moves the sub list items to its parrent list
#' 
#' This function is used to simplify the complex list
#' 
#' @param x A complex list
#' 
#' @return A simplified list
dissolve <- function(x) {
    combi = list()
    operator <- function(x, name = NULL) {
        if (is.list(x)) {
            for (i in seq(x)) {
                operator(x[[i]], names(x)[[i]])
            }
        } else {
            combi[[length(combi) + 1]] <<- x
            names(combi)[[length(combi)]] <<- name
        }
    }
    operator(x)
    combi
}


#' This function is used to compare the similarity of two matrixes 
#' 
#' @param timeseries1 A matrix
#' @param timeseries2 The second matrix to compare with
#' @param index A label to distinguish the result.
#' 
#' @return A vector to show the similarity information about the two matrixes.
similarityBetweenMatrix <- function(timeseries1, 
                                    timeseries2, 
                                    index) {
    if (!identical(dim(timeseries1), dim(timeseries2))) {
        stop("The two matrixes must have the same dimensions")
    }
    
    differ <- abs(timeseries1 - timeseries2)
    
    # correlation<-cor(c(timeseries1),c(timeseries2))
    zerosum <- length(differ[differ == 0])
    correlation <- zerosum/length(differ)
    
    if (is.na(correlation) | is.null(correlation)) {
        correlation <- 0
    }
    
    if (correlation <= 0.2) 
        return(c("veryunlikely", correlation, index))
    
    if (correlation > 0.2 & correlation <= 0.4) 
        return((c("unlikely", correlation, index)))
    
    if (correlation > 0.4 & correlation <= 0.6) 
        return((c("likely", correlation, index)))
    
    if (correlation > 0.6 & correlation <= 0.8) 
        return((c("similar", correlation, index)))
    
    if (correlation > 0.8) 
        return((c("verysimilar", correlation, index)))
}



#' This function is used to check whether or not the data 
#' is the right type for \code{FBNNet}.
#' 
#' @param timeseries_data The timeseries data
#' @return Error or NULL
CheckRightTypeTimeseriesData <- function(timeseries_data) {
    if (!is.list(timeseries_data)) 
        stop("The type of timeseries_data must be LIST")
    
    check <- sapply(timeseries_data, is.matrix)
    if (any(check) == FALSE) 
        stop("The element of the data must be matrix or data.frame")
    
    NULL
}

#' A simple function to check a value is numeric or not
#' @param x A value that need to be checked.
#' @return Error or NULL
checkNumeric <- function(x) {
    if (!is.numeric(x)) 
        stop("The input is not a type of numeric")
    
    NULL
}
#' A simple function to check a value is probability type
#' data.
#' @param x A value that need to be checked.
#' @return Error or NULL
checkProbabilityTypeData <- function(x) {
    if (!is.numeric(x) && (x > 1 || x < 0)) {
        stop("The input is not a type of probability or a value between 0 and 1")
    }
    NULL
}
#' A simple function to check a value is boolean type
#' data.
#' @param x A value that need to be checked.
#' @return TRUE or FALSE
isBooleanTypeTimeseriesData <- function(x) {
    conds <- sapply(x, function(mat) {
        f_mat <- factor(mat)
        if (all(unique(levels(f_mat)) %in% c(0, 1)) || unique(all(levels(f_mat)) %in% c(FALSE, TRUE))) {
            TRUE
        } else {
            FALSE
        }
    })
    
    all(conds)
}

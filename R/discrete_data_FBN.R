#'A method to discrete timeseries data based on the threshod of mean value
#'@param timeSeriesData A list of time series data containing original value
#'@export
discreteTimeSeries <- function(timeSeriesData, method = c("average", "distance")) {
    res <- timeSeriesData
    
    if (is.matrix(timeSeriesData)) 
        fullData <- timeSeriesData else if (is.list(timeSeriesData)) {
        # in case of list, paste all matrices before clustering
        fullData <- do.call(cbind, timeSeriesData)
    } else {
        stop("The type of time seriesdata is not supported")
    }
    
    # std_sampledata<-lapply(1:nrow(fullData),function(rowIndex)sd(fullData[rowIndex,])) #standard deviation names(std_sampledata)<-rownames(fullData)
    
    mean_sampledata <- lapply(1:nrow(fullData), function(rowIndex) mean(fullData[rowIndex, ]))  #average
    names(mean_sampledata) <- rownames(fullData)
    
    # length_sampledata<-length(fullData)
    
    ## sem_sampledata<-std_sampledata/sqrt(length_sampledata) #Standard error of the mean
    
    # max_sampledata<-lapply(1:nrow(fullData),function(rowIndex)max(fullData[rowIndex,])) names(max_sampledata)<-rownames(fullData)
    
    # min_sampledata<-lapply(1:nrow(fullData),function(rowIndex)min(fullData[rowIndex,])) names(min_sampledata)<-rownames(fullData)
    
    
    
    # use average as a thread to discrete time series
    switch(match.arg(method), average = {
        print("Discretizating data using average")
        mean_sampledata <- lapply(1:nrow(fullData), function(rowIndex) mean(fullData[rowIndex, ]))  #average
        names(mean_sampledata) <- rownames(fullData)
        mean_sampledata_max <- do.call(cbind, mean_sampledata)
        
        colmean <- t(mean_sampledata_max)
        if (is.matrix(res)) {
            numofcol <- ncol(fullData)
            matrixMean <- rep.col(colmean, numofcol)
            colNames <- colnames(fullData)
            rowNames <- rownames(fullData)
            binarizedTimeSeries <- sapply(as.data.frame(fullData - (fullData%%matrixMean)), function(vr) as.numeric(vr > 0))
            rownames(binarizedTimeSeries) <- rowNames
            colnames(binarizedTimeSeries) <- colNames
        } else if (is.list(res)) {
            # in case of list, paste all matrices before clustering all matrixs in the list must have the same dimensions
            binarizedTimeSeries <- lapply(res, function(m) {
                numofcol <- ncol(m)
                matrixMean <- rep.col(colmean, numofcol)
                fullData <- abs(m) + 1  #add 1 to avoid 0
                colNames <- colnames(m)
                rowNames <- rownames(m)
                sub_res <- sapply(as.data.frame(fullData - (fullData%%matrixMean)), function(vr) as.numeric(vr > 0))
                rownames(sub_res) <- rowNames
                colnames(sub_res) <- colNames
                return(sub_res)
            })
        }
    }, distance = {
        print("Discretizating data using distance weight")
        mean_sampledata <- lapply(1:nrow(fullData), function(rowIndex) mean(fullData[rowIndex, ]))  #average
        names(mean_sampledata) <- rownames(fullData)
        mean_sampledata_max <- do.call(cbind, mean_sampledata)
        
        colmean <- t(mean_sampledata_max)
        
        weigthMedian <- lapply(1:nrow(fullData), function(rowIndex) median(sort(sqrt((fullData[rowIndex, ] - mean(fullData[rowIndex, ]))^2))))
        names(weigthMedian) <- rownames(fullData)
        weigthMedian_max <- do.call(cbind, weigthMedian)
        
        colmedian <- t(weigthMedian_max)
        
        if (is.matrix(res)) {
            numofcol <- ncol(fullData)
            matrixMean <- rep.col(colmean, numofcol)
            matrixMedian <- rep.col(colmedian, numofcol)
            colNames <- colnames(fullData)
            rowNames <- rownames(fullData)
            binarizedTimeSeries <- sapply(as.data.frame(sqrt((fullData - matrixMean)^2) - matrixMedian), function(vr) as.numeric(vr > 0))
            rownames(binarizedTimeSeries) <- rowNames
            colnames(binarizedTimeSeries) <- colNames
        } else if (is.list(res)) {
            # in case of list, paste all matrices before clustering all matrixs in the list must have the same dimensions
            binarizedTimeSeries <- lapply(res, function(m) {
                numofcol <- ncol(m)
                matrixMean <- rep.col(colmean, numofcol)
                matrixMedian <- rep.col(colmedian, numofcol)
                fullData <- abs(m) + 1  #add 1 to avoid 0
                colNames <- colnames(m)
                rowNames <- rownames(m)
                sub_res <- sapply(as.data.frame(sqrt((fullData - matrixMean)^2) - matrixMedian), function(vr) as.numeric(vr > 0))
                rownames(sub_res) <- rowNames
                colnames(sub_res) <- colNames
                return(sub_res)
            })
        }
    }, stop("'method' must be one of \"average\",\"kmeans\""))
}


## down graded
#' @export
# dividedDiscreteDataintosmallgroups <- function(originalTimeSeriesData, discretedTimeSeriesdata, minElements = 10, maxElements = 30, membexp = 2) {
# findAllLeaves <- function(treegoups) { res1 <- list() len <- length(treegoups) nm <- names(treegoups)[1] if (nm == 'leaf') { res1[length(res1) + 1] <-
# treegoups[1] } else { newlen <- length(res1) + 1 res1[[newlen]] <- list() res1[[newlen]] <- findAllLeaves(treegoups[[1]]) } if (len == 2) { nm <-
# names(treegoups)[2] if (nm == 'leaf') { res1[length(res1) + 1] <- treegoups[2] } else { newlen2 <- length(res1) + 1 res1[[newlen2]] <- list()
# res1[[newlen2]] <- dissolve(findAllLeaves(treegoups[[2]])) } } return(res1) } fuzzygroups <- dividedintosmallgroups(originalTimeSeriesData, minElements,
# maxElements, membexp) groups <- dissolve(findAllLeaves(fuzzygroups)) res <- lapply(groups, function(subgroup, discretedTimeSeriesdata) { subnames <-
# rownames(subgroup) res <- lapply(discretedTimeSeriesdata, function(mtx, subnames) { res1 <- mtx[rownames(mtx) %in% subnames, ] return(res1) }, subnames)
# return(res) }, discretedTimeSeriesdata) names(res) <- c(1:length(res)) class(res) <- 'ClusteredTimeseriesData' return(res) }



DiscretedDataReduction <- function(discretedTimeSeries) {
    numberOfSamples <- length(discretedTimeSeries)
    getAllFalseStatus <- lapply(discretedTimeSeries, function(sample) {
        testM <- as.data.frame(sample)
        stest <- sapply(1:nrow(testM), function(index) all(test[index, ] == 0))
        return(stest)
    })
    matrixFalse <- do.call(cbind, getAllFalseStatus)
    getAllTrueStatus <- lapply(discretedTimeSeries, function(sample) {
        testM <- as.data.frame(sample)
        stest <- sapply(1:nrow(testM), function(index) all(test[index, ] == 1))
        return(stest)
    })
    matrixTrue <- do.call(cbind, getAllTrueStatus)
}

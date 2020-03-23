#'A method to discrete timeseries data based on the threshod of mean value
#'@param timeSeriesData A list of time series data containing original value
#'@export
discreteTimeSeries <- function(timeSeriesData, method = c("average", "distance")) {
    res <- timeSeriesData
    
    if (is.matrix(timeSeriesData)) {
        fullData <- timeSeriesData 
    } 
    else if (is.list(timeSeriesData)) {
        # in case of list, paste all matrices before clustering
        fullData <- do.call(cbind, timeSeriesData)
    } else {
        stop("The type of time seriesdata is not supported")
    }

    nrows <- nrow(fullData)
    mean_sampledata <- lapply(seq_len(nrows), 
                              function(rowIndex) mean(fullData[rowIndex, ]))
    names(mean_sampledata) <- rownames(fullData)
    # use average as a thread to discrete time series
    switch(match.arg(method), average = {
        print("Discretizating data using average")
        nrows <- nrow(fullData)
        mean_sampledata <- lapply(seq_len(nrows), function(rowIndex) mean(fullData[rowIndex, ]))
        names(mean_sampledata) <- rownames(fullData)
        mean_sampledata_max <- do.call(cbind, mean_sampledata)
        
        colmean <- t(mean_sampledata_max)
        if (is.matrix(res)) {
            numofcol <- ncol(fullData)
            matrixMean <- rep.col(colmean, numofcol)
            colNames <- colnames(fullData)
            rowNames <- rownames(fullData)
            binarizedTimeSeries <- sapply(as.data.frame(fullData - (fullData %% matrixMean)), function(vr) as.numeric(vr > 0))
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
                sub_res <- sapply(as.data.frame(fullData - (fullData %% matrixMean)), function(vr) as.numeric(vr > 0))
                rownames(sub_res) <- rowNames
                colnames(sub_res) <- colNames
                return(sub_res)
            })
        }
    }, distance = {
        print("Discretizating data using distance weight")
        nrows <- nrow(fullData)
        mean_sampledata <- lapply(seq_len(nrows), function(rowIndex) mean(fullData[rowIndex, ]))  #average
        names(mean_sampledata) <- rownames(fullData)
        mean_sampledata_max <- do.call(cbind, mean_sampledata)
        
        colmean <- t(mean_sampledata_max)
        nrows <- nrow(fullData)
        weigthMedian <- lapply(seq_len(nrows), function(rowIndex) {
            stats::median(sort(sqrt((fullData[rowIndex, ] - mean(fullData[rowIndex, ]))^2)))
            })
        names(weigthMedian) <- rownames(fullData)
        weigthMedian_max <- do.call(cbind, weigthMedian)
        
        colmedian <- t(weigthMedian_max)
        
        if (is.matrix(res)) {
            numofcol <- ncol(fullData)
            matrixMean <- rep.col(colmean, numofcol)
            matrixMedian <- rep.col(colmedian, numofcol)
            colNames <- colnames(fullData)
            rowNames <- rownames(fullData)
            binarizedTimeSeries <- sapply(as.data.frame(sqrt((fullData - matrixMean)^2) - matrixMedian), function(vr) as.numeric(vr > 
                0))
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
                sub_res <- sapply(as.data.frame(sqrt((fullData - matrixMean)^2) - matrixMedian), function(vr) as.numeric(vr > 
                  0))
                rownames(sub_res) <- rowNames
                colnames(sub_res) <- colNames
                return(sub_res)
            })
        }
    }, stop("'method' must be one of \"average\",\"kmeans\""))
}
#' get Statistic Measures
#' 
#' A method to calculate the statistic measurements based on feature names
#'
#' @param featurenames a list of feature names, i.e., the gene name
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

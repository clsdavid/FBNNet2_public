
#' This method is used to converts normalized
#' timeseries data into a list of samples
#' @param normalizedData An output of the method 
#' normalizeTimesereisRawData that normalized timeseries data
#' @param func A function that specified how to split the column names
#' @param splitor A separator.
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

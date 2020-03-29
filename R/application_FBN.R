#' Generate Fundamental Boolean Model type of Network
#'
#' This is the main entry of the package FBNNet that can be used to mine the
#'  gene regulatory network.
#' 
#' @param timeseries_data A list of timeseries data of samples.
#' @param method Specify a method to discrete the data in the range of
#'  ('kmeans', 'edgeDetector', 'scanStatistic') for the function 
#'  \code{BoolNet::binarizeTimeSeries} to convert the data from numeric value
#'   to boolean value.
#' @param maxK The maximum deep of the Orchard Cube can mine into.
#' @param useParallel optional, by default it is TRUE to run the network
#'  inference algorithm in parallel. FALSE without parallel
#' @param max_deep_temporal, a setting for Temporal Fundamental Boolean model
#'  that specifies the maximum temporal space
#' @param threshold_confidence A threshold of confidence (between 0 and 1) 
#' that used to filter the Fundamental Boolean functions
#' @param threshold_error A threshold of error rate (between 0 and 1) that used
#'  to filter the Fundamental Boolean functions
#' @param threshold_support A threshold of support (between 0 and 1) that used
#'  to filter the Fundamental Boolean functions
#' @param maxFBNRules The maximum rules per type (Activation and Inhibition)
#'  per gene can be mined, the rest will be discarded
#' @param network_only optional for Debug purpose, if TRUE, only output the
#'  networks only, otherwise, output the Orchard cube as well. Warning, 
#' turn off this may cause memory leaking if the number of nodes is too large.
#'  divide them into subgroups.
#' @param verbose Optional, if it is TRUE, then output the logger information to
#' the console.
#' @return An object of a list contains Fundamental Boolean Network, Orchard
#'  cube (optional) and discreted timeseries data
#' @author Leshi Chen, leshi, chen@lincolnuni.ac.nz, chenleshi@hotmail.com
#' @keywords Fundamental Boolean Network
#'
#' @references Chen et al.(2018), Front. Physiol., 25 September 2018, 
#' (\href{https://doi.org/10.3389/fphys.2018.01328}{Front. Physiol.})
#' @references Mussel, Hopfensitz et al. 2010, BoolNet - an R package
#'  for generation, reconstruction and analysis of Boolean networks
#' 
#' @examples  
#' data('yeastTimeSeries')
#' network <- generateFBMNetwork(yeastTimeSeries)
#' network
#' ## draw the general graph
#' FBNNetwork.Graph(network)
#' 
#' network <- generateFBMNetwork(yeastTimeSeries, verbose = TRUE)
#' network
#' network <- generateFBMNetwork(yeastTimeSeries,
#'                               method = "kmeans")
#' network
#' network <- generateFBMNetwork(yeastTimeSeries,
#'                               method = "edgeDetector")
#' network
#' network <- generateFBMNetwork(yeastTimeSeries,
#'                               method = "scanStatistic")
#' network
#' 
#' res <- generateFBMNetwork(yeastTimeSeries, network_only = FALSE)
#' res
#' @export
generateFBMNetwork <- function(
  timeseries_data, 
  method = c("kmeans", "edgeDetector", "scanStatistic"), 
  maxK = 4, 
  useParallel = FALSE, 
  max_deep_temporal = 1, 
  threshold_confidence = 1, 
  threshold_error = 0, 
  threshold_support = 1e-05,
  maxFBNRules = 5, 
  network_only = TRUE,
  verbose = FALSE) {
  if (is.matrix(timeseries_data)) {
    timeseries_data <- list(timeseries_data)
  }
  
  method <- match.arg(method)
  CheckRightTypeTimeseriesData(timeseries_data)
  checkProbabilityTypeData(threshold_confidence)
  checkProbabilityTypeData(threshold_error)
  checkProbabilityTypeData(threshold_support)
  checkNumeric(maxFBNRules)
  if (!is.logical(network_only)) {
    network_only <- TRUE
  }
  if (!is.logical(verbose)) {
    verbose <- TRUE
  }
  if (verbose) {
    # The following code now sets the flog threshold to the log_level defined in the 2nd argument of the function above
    
    futile.logger::flog.threshold(9)
    
    # This could be changed so that every detail is printed (not only errors and warnings). Do this by setting log_level argument to futile.logger::TRACE for full info.
    
    # Now, we create an option of printing and storing results from the checks in somewhere other than the console. Below we are saying
    # that if the log_appender argument is not set to console (and set to a file name) then R will create a file to store the information and outputs from checks. 
    
    # if(log_appender != "console")
    # {
    #   futile.logger::flog.appender(futile.logger::appender.file(log_appender))
    # }
  } else {
    futile.logger::flog.threshold(1)
  }
  
  futile.logger::flog.info(sprintf("Enter generateFBMNetwork zone: 
          method=%s,
          maxK = %s, 
          useParallel = %s, 
          max_deep_temporal = %s,
          threshold_confidence = %s,
          threshold_error = %s,
          threshold_support = %s,
          maxFBNRules = %s,
          network_only = %s,
          verbose = %s", 
    method, 
    maxK, 
    useParallel, 
    max_deep_temporal, 
    threshold_confidence, 
    threshold_error, 
    threshold_support, 
    maxFBNRules,
    network_only,
    verbose))
  if (!isBooleanTypeTimeseriesData(timeseries_data)) {
    timeseries_data <- BoolNet::binarizeTimeSeries(timeseries_data, 
                                                   method = method)$binarizedMeasurements
  }
  genes <- rownames(timeseries_data[[1]])

  futile.logger::flog.info(sprintf("Run generateFBMNetwork with a single cube"))
  time1 <- as.numeric(Sys.time())
  cube <- constructFBNCube(target_genes = genes,
                           conditional_genes = genes,
                           timeseriesCube = timeseries_data,
                           maxK = maxK,
                           temporal = max_deep_temporal,
                           useParallel = useParallel)
  time2 <- as.numeric(Sys.time())
  print(paste("Total cost ", time2 - time1, " seconds to construct a FBN Cube.", sep = "", collapse = ""))
  time1 <- as.numeric(Sys.time())
  network <- mineFBNNetwork(fbnGeneCube = cube, 
                 threshold_confidence  = threshold_confidence,
                 threshold_error = threshold_error,
                 threshold_support = threshold_support,
                 maxFBNRules = maxFBNRules, 
                 useParallel = useParallel)
  time2 <- as.numeric(Sys.time())
  print(paste("Total cost ", time2 - time1, " seconds to mine Fundamental Boolean Functions ", sep = "", collapse = ""))
  
  if (network_only) {
    network
  } else {
    list(cube = cube, network = network)
  }
}


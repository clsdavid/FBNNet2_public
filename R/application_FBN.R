# ~~~~~~~~1~~~~~~~~~2~~~~~~~~~3~~~~~~~~~4~~~~~~~~~5~~~~~~~~~6~~~~~~~~~7~~~~~~~~~8
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
#' @param parallel_on_group optional, if TRUE will run the parallel on the
#'  sub gene groups, otherwise on the genes level
#' @param max_deep_temporal, a setting for Temporal Fundamental Boolean model
#'  that specifies the maximum temporal space
#' @param threshold_confidence A threshod of confidence (between 0 and 1) 
#' that used to filter the Fundamental Boolean functions
#' @param threshold_error A threshod of error rate (between 0 and 1) that used
#'  to filter the Fundamental Boolean functions
#' @param threshold_support A threshod of support (between 0 and 1) that used
#'  to filter the Fundamental Boolean functions
#' @param maxFBNRules The maximum rules per type (Activation and Inhibition)
#'  per gene can be mined, the rest will be discarded
#' @param maxGenesForSingleCube The maximum number of genes for a single
#'  cube. If there are more than the maximum value, the system will 
#' @param runtype The type of return object.
#' @param network_only optional for Debug purpose, if TRUE, only output the
#'  networks only, otherwise, output the Orchard cube as well. Warning, 
#' turn off this may cause memory leaking if the number of nodes is too large.
#'  divide them into subgroups.
#' @return An object of a list contains Fundamental Boolean Network, Orchard
#'  cube (optional) and discreted timeseries data
#' @author Leshi Chen, leshi, chen@lincolnuni.ac.nz, chenleshi@hotmail.com
#' @keywords Fundamental Boolean Network, Boolean Network,
#'  Genetic Regulatory Network
#'
#' @references Chen et al.(2018), Front. Physiol., 25 September 2018, 
#' (\href{https://doi.org/10.3389/fphys.2018.01328}{Front. Physiol.})
#' @references Mussel, Hopfensitz et al. 2010, BoolNet - an R package
#'  for generation, reconstruction and analysis of Boolean networks
#' 
#' @examples  
#' library(BoolNet)
#' data('yeastTimeSeries')
#' network <- generateFBMNetwork(yeastTimeSeries)
#' network
#' @export
generateFBMNetwork <- function(
  timeseries_data, 
  method = c("kmeans", "edgeDetector", "scanStatistic"), 
  maxK = 4, 
  useParallel = FALSE, 
  parallel_on_group = FALSE, 
  max_deep_temporal = 7, 
  threshold_confidence = 1, 
  threshold_error = 0, 
  threshold_support = 1e-05,
  maxFBNRules = 5, 
  maxGenesForSingleCube = 20, 
  runtype = 0, 
  network_only = TRUE) {
  if (is.matrix(timeseries_data)) {
    timeseries_data <- list(timeseries_data)
  }
  
  CheckRightTypeTimeseriesData(timeseries_data)
  checkProbabilityTypeData(threshold_confidence)
  checkProbabilityTypeData(threshold_error)
  checkProbabilityTypeData(threshold_support)
  checkNumeric(maxFBNRules)
  
  futile.logger::flog.info(sprintf("Enter generateFBMNetwork zone: method=%s,
          maxK=%s, 
          useParallel=%s, 
          parallel_on_group=%s,
          max_deep_temporal=%s,
          threshold_confidence=%s,
          threshold_error=%s,
          threshold_support=%s,
          maxFBNRules=%s,
          maxGenesForSingleCube=%s", 
    method, 
    maxK, 
    useParallel, 
    parallel_on_group, 
    max_deep_temporal, 
    threshold_confidence, 
    threshold_error, 
    threshold_support, 
    maxFBNRules, 
    maxGenesForSingleCube))
  if (!isBooleanTypeTimeseriesData(timeseries_data)) {
    timeseries_data <- BoolNet::binarizeTimeSeries(timeseries_data, method = method)$binarizedMeasurements
  }
  genes <- rownames(timeseries_data[[1]])

  futile.logger::flog.info(sprintf("Run generateFBMNetwork with a single cube"))
  cube <- constructFBNCube(genes, genes, timeseries_data, maxK = maxK, temporal = max_deep_temporal, useParallel = useParallel)
  mineFBNNetwork(cube, maxFBNRules = maxFBNRules, useParallel = useParallel)
}

#' A method to convert a vector of gene names to annotated gene details
#' 
#' @param genes A vector of genes
#' @param  filename The name of the output file such as xx.csv
#' @export
output_annotated_genes <- function(genes, filename) {
  ## DAVID_gene_list <- NULL
  utils::data("DAVID_gene_list", overwrite = TRUE)
  mapped_genes <- with(DAVID_gene_list, {
    DAVID_gene_list[DAVID_gene_list$Symbol %in% genes, ]
  })
  distic_mapped_genes <- with(mapped_genes, {
    dplyr::distinct(mapped_genes, Symbol, .keep_all = TRUE)
  })
  utils::write.csv(distic_mapped_genes, file = paste0("temp/", filename))
}

#' Generate Fundamental Boolean Model type of Network
#'
#' This is the main entry of the package FBNNet that can be used to mine the gene regulatory network.
#' 
#' @param timeseries_data A list of timeseries data of samples.
#' @param method an optonal parameter to discrete the data in the range of ("average", "kmeans", "edgeDetector", "scanStatistic"). 
#'   The system applies the function provided by the BoolNet package to convert the data from numeric continually value to boolean value.
#' @param maxK The maximum deep of the Orchard Cube can mine into.
#' @param useParallel optional, by default it is TRUE to run the network inference algorithm in parallel. FALSE without parallel
#' @param parallel_on_group optional, if TRUE will run the parallel on the sub gene groups, otherwise on the genes level
#' @param max_deep_temporal, a setting for Temporal Fundamental Boolean model that specifies the maximum temporal space
#' @param threshold_confidence A threshod of confidence (between 0 and 1) that used to filter the Fundamental Boolean functions
#' @param threshold_error A threshod of error rate (between 0 and 1) that used to filter the Fundamental Boolean functions
#' @param threshold_support A threshod of support (between 0 and 1) that used to filter the Fundamental Boolean functions
#' @param maxFBNRules The maximum rules per type (Activation and Inhibition) per gene can be mined, the rest will be discarded
#' @param maxGenesForSingleCube The maximum number of genes for a single cube. If there are more than the maximum value, the system will 
#' @param network_only optional for Debug purpose, if TRUE, only output the networks only, otherwise, output the Orchard cube as well. Warning, 
#' turn off this may cause memory leaking if the number of nodes is too large.
#'  divide them into subgroups.
#' @param usingFunnayCluster optional to use Fanny cluster algorithm to divide large nodes into groups.
#' @return An object of a list contains Fundamental Boolean Network, Orchard cube (optional) and discreted timeseries data
#' @author Leshi Chen, leshi, chen@lincolnuni.ac.nz, chenleshi@hotmail.com
#' @keywords Fundamental Boolean Network, Boolean Network, Genetic Regulatory Network
#'
#' @references Chen et al.(2018), Front. Physiol., 25 September 2018, 
#' (\href{https://doi.org/10.3389/fphys.2018.01328}{Front. Physiol.})
#' @references Mussel, Hopfensitz et al. 2010, BoolNet - an R package for generation, reconstruction and analysis of Boolean networks
#' 
#' @examples  
#' library(BoolNet)
#' data("yeastTimeSeries")
#' network <- generateFBMNetwork(yeastTimeSeries)
#' network
#' @export
generateFBMNetwork <- function(timeseries_data,
                               method = "kmeans", 
                               maxK = 4,
                               useParallel = FALSE,
                               parallel_on_group = FALSE,
                               max_deep_temporal = 7,
                               threshold_confidence = 1, 
                               threshold_error = 0, 
                               threshold_support = 0.00001, 
                               maxFBNRules = 5,
                               maxGenesForSingleCube = 20,
                               runtype = 0,#0, 1, 2
                               network_only = TRUE) {
  if (is.matrix(timeseries_data)) {
    timeseries_data <- list(timeseries_data)
  }
  
  CheckRightTypeTimeseriesData(timeseries_data)
  checkProbabilityTypeData(threshold_confidence)
  checkProbabilityTypeData(threshold_error)
  checkProbabilityTypeData(threshold_support)
  checkNumeric(maxFBNRules)
  
  futile.logger::flog.info(sprintf("Enter generateFBMNetwork zone: method=%s, maxK=%s, useParallel=%s, 
                                   parallel_on_group=%s, max_deep_temporal=%s, threshold_confidence=%s, threshold_error=%s, threshold_support=%s,
                                   maxFBNRules=%s, maxGenesForSingleCube=%s",
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
  raw_timeseries_data <- timeseries_data
  if (!isBooleanTypeTimeseriesData(timeseries_data))
  {
    if (method == "average") {
      timeseries_data <- discreteTimeSeries(timeseries_data, method = "average")
    } else {
      require(BoolNet)
      timeseries_data <- BoolNet::binarizeTimeSeries(timeseries_data,
                                                     method = method)$binarizedMeasurements
    }
  }
  genes <- rownames(timeseries_data[[1]])
  if (runtype == 0) {
    futile.logger::flog.info(sprintf("Run generateFBMNetwork with a single cube"))
    cube <- constructFBNCube(genes, 
                                genes, 
                                timeseries_data, 
                                maxK = maxK, 
                                temporal = max_deep_temporal, 
                                useParallel = useParallel)
    network <- mineFBNNetwork(cube, 
                              maxFBNRules = maxFBNRules,
                              useParallel = useParallel)
    return(network)
  }
  
  if (runtype == 1) {
    futile.logger::flog.info(sprintf("Run generateFBMNetwork with cimbined clustering"))
    cube_cluster <- clusterdDiscreteData(raw_timeseries_data,
                                              timeseries_data,
                                              minElements = 2,
                                              maxElements = maxGenesForSingleCube)
    this_cube_cluster <- constructFBNCubeAndNetworkInClusters(cube_cluster,
                                                              maxK = maxK,
                                                              temporal = max_deep_temporal,
                                                              useParallel = useParallel,
                                                              parallel_on_group = parallel_on_group,
                                                              output_cube = FALSE)
  } else if (runtype == 2) {
    futile.logger::flog.info(sprintf("Run generateFBMNetwork with combined sub groups"))
    cube_cluster <-  dividedDataIntoSubgroups(timeseries_data, maxGenesForSingleCube, maxK = maxK)
    this_cube_cluster <-  constructFBNCubeAndNetworkInClusters(cube_cluster,
                                                               maxK = maxK,
                                                               temporal = max_deep_temporal,
                                                               useParallel = useParallel,
                                                               parallel_on_group = parallel_on_group,
                                                               output_cube = FALSE)
  } else if (runtype == 3) {
    futile.logger::flog.info(sprintf("Run generateFBMNetwork with small target gene groups"))
    cube_cluster <-  dividedDataIntoSmallgroups(timeseries_data, maxGenesForSingleCube, maxK = maxK)
    this_cube_cluster <-  constructFBNCubeAndNetworkInSmallGroups(cube_cluster,
                                                                  maxK = maxK,
                                                                  temporal = max_deep_temporal,
                                                                  useParallel = useParallel,
                                                                  parallel_on_group = parallel_on_group,
                                                                  output_cube = FALSE)
  } else {
    stop("The runtype is not supported")
  }


    # merge all networks? seperate by clusters?
    network <- mergeClusterNetworks(clusteredFBNCube = this_cube_cluster, 
                                    threshold_error = threshold_error, 
                                    maxFBNRules = maxFBNRules)
    futile.logger::flog.info(sprintf("Leave generateFBMNetwork"))
    final_network <- filterNetworkConnections(network)
    sink(file = "temp/final_network.txt", type = "output")
    final_network
    sink()
    data("DAVID_gene_list")
    mapped_genes <- DAVID_gene_list[DAVID_gene_list$Symbol %in% final_network$genes,]
    distic_mapped_genes <- dplyr::distinct(mapped_genes, Symbol, .keep_all= TRUE)
    write.csv(distic_mapped_genes, file="temp/annotated_differencial_genes.csv")
    final_network
}

#' A method to convert a vector of gene names to annotated gene details
#' 
#' @param genes A vector of genes
#' @param  filename The name of the output file such as xx.csv
#' @export
output_annotated_genes <- function(genes, filename) {
  data("DAVID_gene_list")
  mapped_genes <- DAVID_gene_list[DAVID_gene_list$Symbol %in% genes,]
  distic_mapped_genes <- dplyr::distinct(mapped_genes, Symbol, .keep_all= TRUE)
  write.csv(distic_mapped_genes, file= paste0("temp/", filename))
}
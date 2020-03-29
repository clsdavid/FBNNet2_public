#' Create an Orchard cube
#'
#' This is the main function(s) to genereate a single Orchard Cube or a group
#' of cubes
#'
#' @param target_genes A vector of genes that will be treated as target genes
#' @param conditional_genes All genes that are available for building up the
#'  cube
#' @param timeseriesCube A list of samples in which a sample is a matrix that
#'  contains gene states 
#'  where genes in rows and time points in columns.
#' @param maxK The maximum level the cube can dig in
#' @param temporal A value that used to be 1 indicates the previous steps the
#'  current one can depend on
#' @param useParallel If it is TRUE, the constructing will run it in parallel,
#'  otherwise in a singl thread
#' @return An Orchard cube that contains all precomputed measures
#' @examples
#' require(BoolNet)
#' data('ExampleNetwork')
#' initialStates <- generateAllCombinationBinary(ExampleNetwork$genes)
#' trainingseries <- genereateBoolNetTimeseries(ExampleNetwork,
#'                                            initialStates,
#'                                            43,
#'                                            type='synchronous')
#' cube<-constructFBNCube(target_genes = ExampleNetwork$genes,
#'                        conditional_genes = ExampleNetwork$genes,
#'                        timeseriesCube = trainingseries,
#'                        maxK = 4,
#'                        temporal = 1,
#'                        useParallel = FALSE)
#' network <- mineFBNNetwork(cube,ExampleNetwork$genes)
#' network
#' ## draw the general graph
#' FBNNetwork.Graph(network)
#' @export
constructFBNCube <- function(target_genes, 
                             conditional_genes, 
                             timeseriesCube, 
                             maxK = 5, 
                             temporal = 1, 
                             useParallel = FALSE) {
  CheckRightTypeTimeseriesData(timeseriesCube)
  checkNumeric(maxK)
  checkNumeric(temporal)
  if (!is.logical(useParallel)) {
    useParallel <- FALSE
  }
  if (is.null(target_genes) || 
     !is.character(target_genes) || 
     length(target_genes) == 0 ||
     !is.vector(target_genes)) {
    stop("The 'target_genes' cannot be NULL or EMPTY or other type rather than Character")
  }
  if (is.null(conditional_genes) || 
      !is.character(conditional_genes) || 
      length(conditional_genes) == 0 ||
      !is.vector(conditional_genes)) {
    stop("The 'conditional_genes' cannot be NULL or EMPTY or other type rather than Character")
  }
  futile.logger::flog.info(sprintf("Enter constructFBNCube zone: 
        target_genes=%s genes and they are %s,
        conditional_genes=%s genes and they are %s,
        data_length=%s,
        maxK=%s, 
        temporal=%s,
        useParallel=%s", 
    length(target_genes), 
    paste(target_genes, sep = ", ", collapse = ", "), 
    length(conditional_genes),
    paste(conditional_genes, sep = ", ", collapse = ", "), 
    length(timeseriesCube), 
    maxK, 
    temporal, 
    useParallel))
  
  ## construct gene tree by timeseries * samples * timepoints(columns)
  ## divid into sub groups
  
  internalloopByWhole <- function(i, target_genes, conditional_genes, maxK, temporal, mainParameters) {
    target_gene <- target_genes[[i]]
    process_cube_algorithm(target_gene, conditional_genes, maxK, temporal, mainParameters)
  }
  
  # try to improve the performance using vector
  res <- vector("list", length = length(target_genes))
  # need to verify this that may use removeduplicateCol
  reducedCube <- FBNDataReduction(timeseriesCube)
  
  getCurrentStates <- vector("list", length = temporal)
  getpreviousStates <- vector("list", length = temporal)
  getCurrentStates_c <- vector("list", length = temporal)
  getpreviousStates_c <- vector("list", length = temporal)
  # set up data sheet for each temporal
  index <- temporal
  while (index > 0) {
    getCurrentStates[[index]] <- extractGeneStateFromTimeSeriesCube(reducedCube, index)
    getpreviousStates[[index]] <- getCurrentStates[[index]]
    getCurrentStates_c[[index]] <- getCurrentStates[[index]]
    getpreviousStates_c[[index]] <- getCurrentStates[[index]]
    index <- index - 1
  }
  
  total_timepoints <- sum(vapply(reducedCube, function(timeshet) ncol(timeshet), integer(1)))
  total_samples <- length(reducedCube)
  all_gene_names <- rownames(reducedCube[[1]])
  
  mainParameters <- new.env(hash = TRUE, parent = globalenv())
  mainParameters$currentStates <- getCurrentStates
  mainParameters$previousStates <- getpreviousStates
  mainParameters$currentStates_c <- getCurrentStates_c
  mainParameters$previousStates_c <- getpreviousStates_c
  mainParameters$total_samples <- total_samples
  mainParameters$all_gene_names <- all_gene_names
  mainParameters$total_timepoints <- total_timepoints
  if (useParallel) {
    res <- doParallelWork(internalloopByWhole, target_genes, conditional_genes, maxK, temporal, mainParameters)
    
  } else {
    res <- doNonParallelWork(internalloopByWhole, target_genes, conditional_genes, maxK, temporal, mainParameters)
  }
  
  if (!is.null(res) && length(res) > 0) {
    cond1 <- vapply(res, function(entry) !is.null(entry), logical(1))
    res <- (res[cond1][unlist(lapply(res[cond1], length) != 0)])
    class(res) <- c("FBNCube")
  }

  if (useParallel) {
    closeAllConnections()
  }
  
  futile.logger::flog.info("Leave constructFBNCube zone.")
  res
}

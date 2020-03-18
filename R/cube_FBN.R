#' Create an Orchard cube
#'
#' This is the main function(s) to genereate a sigle Orchard Cube or a group
#' of cubes
#'
#' @param target_genes A vector of genes that will be treated as target genes
#' @param conditional_genes All genes that are available for building up the
#'  cube
#' @param timeseriesCube A list of samples in which a sample is a matrix that
#'  contains gene states 
#'  where genes in rows and time points in columms
#' @param maxK The maximum level the cube can dig in
#' @param temporal A value that used to be 1 indicates the previous steps the
#'  current one can depend on
#' @param useParallel If it is TRUE, the constructing will run it in parallel,
#'  otherwise in a singl thread
#' @param clusteredTimeseriesData clustered timeseries data
#' @return A Orcahrd cube that contains all precomputed measures
#' @examples
#' mat1<-matrix(c('1','0','0','1','0','0','0','1','1'),
#'              3,3, dimnames=list(c('gene1','gene2','gene3'),c('1','2','3')))
#' mat2<-matrix(c('1','1','0','1','0','1','1','1','0'),
#'              3,3, dimnames=list(c('gene1','gene2','gene3'),c('1','2','3')))
#' listtest<-list(mat1,mat2)
#' constructFBNCube(c('gene1','gene2'),
#'                  c('gene1','gene2','gene3'),
#'                  listtest, 4, 1, FALSE)
#' 
#' @export
constructFBNCube <- function(target_genes,
                             conditional_genes,
                             timeseriesCube,
                             maxK = 5,
                             temporal = 1,
                             useParallel = FALSE) {
    futile.logger::flog.info(
        sprintf("Enter constructFBNCube zone: 
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
    
    # construct gene tree by timeseries * samples * timepoints(columns)
    internalloopByWhole <- function(i,
                                    target_genes,
                                    conditional_genes,
                                    maxK,
                                    temporal,
                                    mainParameters) {
        target_gene <- target_genes[[i]]
        process_cube_algorithm(target_gene,
                               conditional_genes,
                               maxK,
                               temporal,
                               mainParameters)
    }
    
    time1 <- as.numeric(Sys.time())
    #try to improve the performance using vector
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
        getCurrentStates[[index]] <- 
            extractGeneStateFromTimeSeriesCube(reducedCube, index)
        getpreviousStates[[index]] <- getCurrentStates[[index]]
        getCurrentStates_c[[index]] <- getCurrentStates[[index]]
        getpreviousStates_c[[index]] <- getCurrentStates[[index]]
        index <- index - 1
    }
    
    total_timepoints <- sum(vapply(reducedCube, function(timeshet) ncol(timeshet),
                                   integer(1)));
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
        res <- doParallelWork(internalloopByWhole,
                              target_genes,
                              conditional_genes,
                              maxK,
                              temporal,
                              mainParameters)
        
    } else {
        res <- doNonParallelWork(internalloopByWhole, 
                                 target_genes,
                                 conditional_genes,
                                 maxK, temporal,
                                 mainParameters)
    }
    
    if (!is.null(res) & length(res) > 0) {
        cond1 <- sapply(res, function(entry) !is.null(entry))
        res <- (res[cond1][unlist(lapply(res[cond1], length) != 0)])
        class(res) <- c("FBNCube")
    }
    time2 <- as.numeric(Sys.time())
    print(paste("Total cost ",
                time2 - time1,
                " seconds to construct gene cube with size of ",
                ((object.size(res)/1024)/1024), 
        sep = "", collapse = ""))
    
    if (useParallel) {
        closeAllConnections()
    }
    futile.logger::flog.info("Leave constructFBNCube zone.")
    res
}



#' @rdname constructFBNCube
#' @export
constructFBNCubeAndNetworkInSmallGroups <- function(
    groupedTimeseriesData,
    maxK = 5,
    temporal = 1,
    useParallel = FALSE,
    threshold_confidence = 1, 
    threshold_error = 0,
    threshold_support = 1e-05,
    maxFBNRules = 20, 
    output_cube = FALSE, 
    parallel_on_group = FALSE) {
    futile.logger::flog.info(
        sprintf("Enter constructFBNCubeAndNetworkInSmallGroups zone:
                clusteredTimeseriesData=%s,
                maxK=%s,
                useParallel=%s,
                threshold_confidence=%s,
                threshold_error=%s,
                threshold_support=%s,
                maxFBNRules=%s,
                parallel_on_group=%s", 
        length(groupedTimeseriesData), 
        maxK,
        useParallel,
        threshold_confidence,
        threshold_error, 
        threshold_support, 
        maxFBNRules, 
        parallel_on_group))
    
    if (parallel_on_group == TRUE) {
        ## ToDo need to find out how to do parallel inside parallel
        useParallel = FALSE  
    }
    
    ## get clusters
    combined_clusters <- groupedTimeseriesData$clusters
    timeseries_data <- groupedTimeseriesData$discretedTimeSeriesdata
    conditional_genes <- groupedTimeseriesData$conditional_genes
    ## build combination clusters
    print("Starting.....")
    if (parallel_on_group) {
        clusteredFBNCube <- doParallelWork(
            internal_constructCubeWithGroups, 
            combined_clusters,
            conditional_genes, 
            timeseries_data, 
            maxK, 
            temporal,
            useParallel,
            threshold_confidence,
            threshold_error,
            threshold_support,
            maxFBNRules, 
            output_cube)
        
    } else {
        clusteredFBNCube <- doNonParallelWork(
            internal_constructCubeWithGroups,
            combined_clusters, 
            conditional_genes,
            timeseries_data, 
            maxK,
            temporal,
            useParallel,
            threshold_confidence,
            threshold_error, 
            threshold_support, 
            maxFBNRules, 
            output_cube)
    }
    
    
    class(clusteredFBNCube) <- "ClusteredFBNCube"
    print("End.....")
    futile.logger::flog.info(
        sprintf("Leave constructFBNCubeAndNetworkInSmallGroups zone."))
    clusteredFBNCube
}


#' @rdname 'constructFBNCube'
#' @export
constructFBNCubeAndNetworkInClusters <- function(clusteredTimeseriesData, maxK = 5, temporal = 1, useParallel = FALSE, threshold_confidence = 1, 
    threshold_error = 0, threshold_support = 1e-05, maxFBNRules = 20, output_cube = FALSE, parallel_on_group = FALSE) {
    futile.logger::flog.info(sprintf("Enter constructFBNCubeAndNetworkInClusters zone: clusteredTimeseriesData=%s, maxK=%s, useParallel=%s, threshold_confidence=%s, threshold_error=%s, threshold_support=%s, maxFBNRules=%s, parallel_on_group=%s", 
        length(clusteredTimeseriesData), maxK, useParallel, threshold_confidence, threshold_error, threshold_support, maxFBNRules, 
        parallel_on_group))
    
    if (parallel_on_group == TRUE) {
        useParallel = FALSE  ## ToDo need to find out how to do parallel inside parallel
    }
    
    ## get clusters
    combined_clusters <- clusteredTimeseriesData$combined_clusters
    timeseries_data <- clusteredTimeseriesData$discretedTimeSeriesdata
    target_genes <- clusteredTimeseriesData$target_genes
    ## build combination clusters
    print("Starting.....")
    if (parallel_on_group) {
        clusteredFBNCube <- doParallelWork(internal_constructCubeWithClusters, combined_clusters, target_genes, timeseries_data, 
            maxK, temporal, useParallel, threshold_confidence, threshold_error, threshold_support, maxFBNRules, output_cube)
        
    } else {
        clusteredFBNCube <- doNonParallelWork(internal_constructCubeWithClusters, combined_clusters, target_genes, timeseries_data, 
            maxK, temporal, useParallel, threshold_confidence, threshold_error, threshold_support, maxFBNRules, output_cube)
    }
    
    
    class(clusteredFBNCube) <- "ClusteredFBNCube"
    print("End.....")
    futile.logger::flog.info(sprintf("Leave constructFBNCubeAndNetworkInClusters zone."))
    clusteredFBNCube
}

#' @noRd
internal_constructCubeWithClusters <- function(index, combined_clusters, target_genes, timeseries_data, maxK, temporal, useParallel, 
    threshold_confidence, threshold_error, threshold_support, maxFBNRules, output_cube) {
    res <- list()
    # get row values i.e. genes
    genes <- combined_clusters[[index]]
    futile.logger::flog.info(sprintf("enter internal_constructCubeWithClusters: genes=%s", paste(genes, sep = ", ", collapse = ", ")))
    
    if (output_cube) {
        res[[1]] <- constructFBNCube(target_genes, genes, timeseries_data, maxK, temporal, useParallel)
        names(res)[[1]] <- "Cube"
        
        res[[2]] <- searchFBNNetworkCore(res[[1]], genes = target_genes, threshold_confidence = threshold_confidence, threshold_error = threshold_error, 
            threshold_support = threshold_support, maxFBNRules = maxFBNRules, useParallel = useParallel)
        names(res)[[2]] <- "NetworkCores"
        
        res[[3]] <- genes
        names(res)[[3]] <- "Genes"
        
    } else {
        this_cube <- constructFBNCube(target_genes, genes, timeseries_data, maxK, temporal, useParallel)
        res[[1]] <- searchFBNNetworkCore(this_cube, genes = target_genes, threshold_confidence = threshold_confidence, threshold_error = threshold_error, 
            threshold_support = threshold_support, maxFBNRules = maxFBNRules, useParallel = useParallel)
        names(res)[[1]] <- "NetworkCores"
        
        res[[2]] <- genes
        names(res)[[2]] <- "Genes"
        
    }
    
    futile.logger::flog.info(sprintf("leave internal_constructCubeWithClusters."))
    list(index = res)
}

#' @noRd
internal_constructCubeWithGroups <- function(index, combined_clusters, conditional_genes, timeseries_data, maxK, temporal, 
    useParallel, threshold_confidence, threshold_error, threshold_support, maxFBNRules, output_cube) {
    res <- list()
    # get row values i.e. genes
    genes <- combined_clusters[[index]]
    futile.logger::flog.info(sprintf("enter internal_constructCubeWithClusters: genes=%s", paste(genes, sep = ", ", collapse = ", ")))
    
    if (output_cube) {
        res[[1]] <- constructFBNCube(genes, conditional_genes, timeseries_data, maxK, temporal, useParallel)
        names(res)[[1]] <- "Cube"
        
        res[[2]] <- searchFBNNetworkCore(res[[1]], genes = genes, threshold_confidence = threshold_confidence, threshold_error = threshold_error, 
            threshold_support = threshold_support, maxFBNRules = maxFBNRules, useParallel = useParallel)
        names(res)[[2]] <- "NetworkCores"
        
        res[[3]] <- genes
        names(res)[[3]] <- "Genes"
        
    } else {
        this_cube <- constructFBNCube(genes, conditional_genes, timeseries_data, maxK, temporal, useParallel)
        res[[1]] <- searchFBNNetworkCore(this_cube, genes = genes, threshold_confidence = threshold_confidence, threshold_error = threshold_error, 
            threshold_support = threshold_support, maxFBNRules = maxFBNRules, useParallel = useParallel)
        names(res)[[1]] <- "NetworkCores"
        
        res[[2]] <- genes
        names(res)[[2]] <- "Genes"
        
    }
    
    futile.logger::flog.info(sprintf("leave internal_constructCubeWithClusters."))
    list(index = res)
}

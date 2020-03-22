
#~~~~~~~~1~~~~~~~~~2~~~~~~~~~3~~~~~~~~~4~~~~~~~~~5~~~~~~~~~6~~~~~~~~~7~~~~~~~~~8
#'ClusterTimeSeries is a method to cluster timeseries data into small sets of
#' timeseries data
#' 
#' @param timeSeriesData A list of timeseries matrix (genes on row and time
#'  points on columns)
#' @param originalTimeSeriesData The original timeseries data that contains
#'  continual values
#' @param discretedTimeSeriesdata Discreted timeseries data that contains
#'  boolean values
#' @param method Clustering method
#' @param modelParameters An object of ModelParameter
#' @param requireddata The expression data that are going to be clustered.
#'  The data must have time steps on colums and genes on rows
#' @param minElements Each cluster must have elements more or equal than the
#'  value specified
#' @param maxElements Each cluster must have elements less or equal than the
#'  value specified
#' @param membexp A value for fussy algorithm to specify the degree of fussy.
#'  Big number, big fussy
#' @param useParallel Turned on parallel
#' 
#'@return Clustered timeseries data
#'@examples
#' ## kmeansParameters<-list(type='kmeans',
#' ## numOfClusters=10,
#' ## nstart=100,
#' ## iter.max=1000)
#' ## class(kmeansParameters)='ModelParameter'
#' ## hierarchicalParameters<-list(type='hierarchical',
#' ## distmethod='euclidean',
#' ## hclustmethod='ward.D2')
#' ## class(hierarchicalParameters)='ModelParameter'
#' ## nbclustParameters<-list(type='nbclust',
#' ## distmethod='euclidean',
#' ## min.nc=2,
#' ## max.nc=10,
#' ## nbmethod='complete',
#' ## nbindex='all')
#' ## class(nbclustParameters)='ModelParameter'
#' ## dianaParameters<-list(type='diana',
#' ## numOfClusters=10,palette='jco')
#' ## class(dianaParameters)='ModelParameter'
#' ## fuzzyParameters<-list(type='fanny',
#' ## metric = 'euclidean', stand = FALSE)
#' ## class(fuzzyParameters)='ModelParameter'
#'
#' ## sortedtimeseries<-constructTestGSE2677Data('Your Cel Data folder',
#' ## useGCRMA=FALSE)
#' ## clustered_kmeans<-clusterTimeSeries(sortedtimeseries,
#' ## method='kmeans',
#' ## kmeansParameters)
#' ## clustered_hierarchical<-clusterTimeSeries(sortedtimeseries,
#' ## method='hierarchical',
#' ## hierarchicalParameters)
#' ## clustered_nbclust<-clusterTimeSeries(sortedtimeseries,
#' ## method='nbclust',
#' ## nbclustParameters)
#' ## clustered_diana<-clusterTimeSeries(sortedtimeseries,
#' ## method='diana',
#' ## dianaParameters)
#' ## clusteredFanny<-clusterTimeSeries(sortedtimeseries,
#' ## method='fanny',
#' ## fuzzyParameters)

#' @export
clusterTimeSeries <- function(timeSeriesData,
                              method = c("kmeans",
                                         "hierarchical",
                                         "diana",
                                         "fanny",
                                         "nbclust"),
                              modelParameters) {
    # require('cluster') require('factoextra') require('magrittr')
    if (is.matrix(timeSeriesData)) 
        clusterData <- timeSeriesData 
    else 
        stop("The type of timeSeriesData must be matrix")
    
    # switch between the different methods
    switch(match.arg(method), kmeans = {
        set.seed(123)
        if (!inherits(modelParameters, "ModelParameter")) 
            stop("modelParameters must be inherited from ModelParameter")
        
        if (modelParameters$type != "kmeans") 
            stop("The type of the ModelParameter cannot be applied to kmeans")
        
        numOfCluster <- modelParameters$numOfClusters
        if (is.na(numOfCluster) | is.null(numOfCluster)) 
            stop("The parameter numOfCluster is missing")
        
        nstart <- modelParameters$nstart
        iter.max <- modelParameters$iter.max
        
        thisCluster <- na.omit(clusterData)
        thisCluster <- scale(thisCluster)
        thisCluster <- kmeans(thisCluster,
                              numOfCluster,
                              nstart = nstart,
                              iter.max = iter.max)
        
        # visualize the cluster
        factoextra::fviz_cluster(thisCluster,
                                 clusterData,
                                 ellipse.type = "convex",
                                 palette = "jco",
                                 ggtheme = ggplot2::theme_minimal())
        
        return(list(type = "kmeans", cluster = thisCluster))
    }, hierarchical = {
        # It does not require to pre-specify the number of clusters to be generated.
        if (!inherits(modelParameters, "ModelParameter")) 
            stop("modelParameters must be inherited from ModelParameter")
        
        if (modelParameters$type != "hierarchical") 
            stop("The type of the ModelParameter cannot be applied to hierarchical")
        
        distmethod <- modelParameters$distmethod
        hclustmethod <- modelParameters$hclustmethod
        numOfCluster <- modelParameters$numOfClusters
        if (is.na(numOfCluster) | is.null(numOfCluster)) 
            stop("The parameter numOfCluster is missing")
        
        thisCluster <- na.omit(clusterData)
        thisCluster <- scale(thisCluster)
        thisCluster <- dist(thisCluster, method = distmethod)
        # Compute hierachical clustering
        thisCluster <- hclust(thisCluster, 
                              method = hclustmethod)  
        
        clusterCut <- cutree(thisCluster, numOfCluster)
        
        # Visualize using factoextra Cut in 4 groups and color by groups
        library(factoextra)
        # Add rectangle around groups
        fviz_dend(thisCluster,
                  k = numOfCluster,
                  cex = 0.5, 
                  k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
                  color_labels_by_k = TRUE, 
                  rect = TRUE  )
        
        return(list(type = match.arg(method),
                    cluster = thisCluster, 
                    clusterCut = clusterCut))
        
    }, diana = {
        if (!inherits(modelParameters, "ModelParameter")) 
            stop("modelParameters must be inherited from ModelParameter")
        
        if (modelParameters$type != "diana") 
            stop("The type of the ModelParameter cannot be applied to diana")
        
        numOfCluster <- modelParameters$numOfClusters
        if (is.na(numOfCluster) | is.null(numOfCluster)) 
            stop("The parameter numOfCluster is missing")
        
        palette <- modelParameters$palette
        if (is.na(numOfCluster) || is.null(numOfCluster)) {
            palette <- "jco"
        }
        # Compute diana()
        library(cluster)
        thisCluster <- na.omit(clusterData)
        thisCluster <- scale(thisCluster) 
        thisCluster <- diana(thisCluster,
                             stand = TRUE)
        
        clusterCut <- cutree(thisCluster, numOfCluster)
        
        # Plot the dendrogram Cut in four groups Color palette
        library(factoextra)
        fviz_dend(thisCluster,
                  cex = 0.5,
                  k = numOfCluster,
                  palette = palette)
        
        return(list(type = match.arg(method), 
                    cluster = thisCluster,
                    clusterCut = clusterCut))
    }, fanny = {
        # fuzzy cluster
        library(cluster)
        if (!inherits(modelParameters, "ModelParameter")) 
            stop("modelParameters must be inherited from ModelParameter")
        
        if (modelParameters$type != "fanny") 
            stop("The type of the ModelParameter cannot be applied to fanny")
        
        metric <- modelParameters$metric
        numOfCluster <- modelParameters$numOfClusters
        stand <- modelParameters$stand  #True or false
        membexp <- modelParameters$membexp
        # Remove missing values (NA) Scale the data # standardize variables
        clusterData <-  na.omit(clusterData)
        clusterData <- scale(clusterData)
        # Compute fuzzy clustering with k = 2
        thisCluster <- fanny(clusterData,
                             k = numOfCluster,
                             metric = metric,
                             stand = stand,
                             memb.exp = membexp)  
        
        factoextra::fviz_cluster(thisCluster,
                                 ellipse.type = "norm",
                                 repel = TRUE, 
                                 palette = "jco",
                                 ggtheme = ggplot2::theme_minimal(), 
            legend = "right")
        return(list(type = match.arg(method),
                    cluster = thisCluster))
        
    }, nbclust = {
        # Determining the optimal number of clusters
        if (!inherits(modelParameters, "ModelParameter")) 
            stop("modelParameters must be inherited from ModelParameter")
        
        if (modelParameters$type != "nbclust") 
            stop("The type of the ModelParameter cannot be applied to nbclust")
        
        #distance method like euclidean
        distmethod <- modelParameters$distmethod 
        min.nc <- modelParameters$min.nc
        max.nc <- modelParameters$max.nc
        nbmethod <- modelParameters$nbmethod
        nbindex <- modelParameters$nbindex
        # Remove missing values (NA) standardize variables
        thisCluster <- na.omit(clusterData)
        thisCluster <- scale(thisCluster)
        thisCluster <- NbClust::NbClust(thisCluster, 
                               distance = distmethod,
                               min.nc = min.nc,
                               max.nc = max.nc, 
                               method = nbmethod,
                               index = nbindex)
        
        factoextra::fviz_nbclust(thisCluster, 
                                 ggtheme = ggplot2::theme_minimal())
        
        return(list(type = match.arg(method), cluster = thisCluster))
        
    }, mclust = {
        # model based cluster
        
    }, stop("'method' must be one of \"kmeans\",\"edgeDetector\",\"scanStatistic\""))
}

fuzzyTreeCluster <- function(requireddata, 
                             minElements,
                             maxElements,
                             membexp = 3) {
    
    internalloop <- function(i,
                             groupsdata,
                             minElements,
                             maxElements,
                             fuzzyParameters) {
        
        
        subgroupData <- groupsdata[[i]]
        res <- list()
        rowValues <- rownames(subgroupData)
        lenOfsub <- length(rowValues)
        
        if (lenOfsub <= maxElements) {
            res[[i]] <- subgroupData
            names(res)[[i]] <- i
            return(res)
        }
        
        tryCatch({
            clusterobject <- clusterTimeSeries(subgroupData,
                                               method = "fanny",
                                               fuzzyParameters)
            clusters <- clusterobject$cluster$clustering
            
            
            fac <- factor(clusters)
            dat <- data.frame(cluster = clusters)
            groups <- split.data.frame(dat, dat$cluster)
            
            # ToDo: calculate group data
            groupsdata <- lapply(groups, function(g, subgroupData) {
                gnames <- rownames(g)
                return(subgroupData[rownames(subgroupData) %in% gnames, ])
            }, subgroupData)
            
            result <- doNonParallelWork(internalloop,
                                        groupsdata,
                                        minElements,
                                        maxElements,
                                        fuzzyParameters)
            conds <- sapply(result, function(x) !is.null(x))
            result <- result[conds]
            res[[i]] <- result
            names(res)[[i]] <- "sub"
            return(res)
        }, error = function(e) {
            emas <- e
            stop(sprintf("Error executing FBN Cluster: %s", e$message))
        })
    }
    
    if (is.matrix(requireddata)) {
        fullData <- requireddata
    } else if (is.list(requireddata)) {
        # in case of list, paste all matrices before clustering
        fullData <- do.call(cbind, requireddata)
    } else {
        stop("The type of time seriesdata is not supported")
    }
    
    rowValues <- rownames(fullData)
    
    # cluster parameter defined the tree always divided by two
    fuzzyParameters <- list(type = "fanny", 
                            metric = "euclidean",
                            stand = FALSE, 
                            numOfClusters = 2, 
                            membexp = membexp)
    class(fuzzyParameters) = "ModelParameter"
    
    clusterobject <- clusterTimeSeries(fullData, 
                                       method = "fanny",
                                       fuzzyParameters)
    
    clusters <- clusterobject$cluster$clustering
    
    fac <- factor(clusters)
    dat <- data.frame(cluster = clusters)
    groups <- split.data.frame(dat, dat$cluster)
    
    groupsdata <- list()
    # ToDo: calculate initial group data
    groupsdata <- lapply(groups, function(g, fullData) {
        gnames <- rownames(g)
        return(fullData[rownames(fullData) %in% gnames, ])
    }, fullData)
    
    
    
    time1 <- as.numeric(Sys.time())
    
    res <- list()
    
    if (length(rowValues) == 0) {
        return(list())
    }
    
    res <- doNonParallelWork(internalloop,
                             groupsdata, 
                             minElements, 
                             maxElements, 
                             fuzzyParameters)
    
    
    cond1 <- vapply(res, function(entry) !is.null(entry) && !is.na(entry),
                    logical(1))
    res <- (res[cond1][unlist(lapply(res[cond1], length) != 0)])
    time2 <- as.numeric(Sys.time())
    print(paste("Total cost ",
                time2 - time1, 
                " seconds to run FuzzyTree algorithm ",
                sep = "", 
                collapse = ""))
    return(res)
}


#' internal for calling from the function clusterdDiscreteData
#' @noRd
dividedintosmallgroups <- function(timeSeriesData, 
                                   minElements = 10,
                                   maxElements = 30, 
                                   membexp = 2) {
    if (is.matrix(timeSeriesData)) {
        fullData <- timeSeriesData
    } else if (is.list(timeSeriesData)) {
        fullData <- do.call(cbind, timeSeriesData)
    } else {
        stop("The type of time seriesdata is not supported")
    }
    fuzzgroups <- fuzzyTreeCluster(fullData,
                                   minElements,
                                   maxElements,
                                   membexp)
    return(fuzzgroups)
}


#' @export
clusterdDiscreteData <- function(originalTimeSeriesData,
                                 discretedTimeSeriesdata,
                                 minElements = 10,
                                 maxElements = 30,
                                 membexp = 2) {
    findAllLeaves <- function(treegoups) {
        res1 <- list()
        len <- length(treegoups)
        
        nm <- names(treegoups)[1]
        if (nm == "sub") {
            res1[length(res1) + 1] <- treegoups[1]
        } else {
            newlen <- length(res1) + 1
            res1[[newlen]] <- list()
            res1[[newlen]] <- findAllLeaves(treegoups[[1]])
        }
        
        if (len == 2) {
            nm <- names(treegoups)[2]
            if (nm == "sub") {
                res1[length(res1) + 1] <- treegoups[2]
            } else {
                newlen2 <- length(res1) + 1
                res1[[newlen2]] <- list()
                res1[[newlen2]] <- dissolve(findAllLeaves(treegoups[[2]]))
            }
        }
        conds <- sapply(res1, function(x) !is.null(x))
        return(res1[conds])
    }
    
    conds <- sapply(originalTimeSeriesData, 
                    function(mat) any(duplicated(rownames(mat))))
    genes <- unique(unlist(lapply(discretedTimeSeriesdata, rownames)))
    if (any(conds)) 
        stop("duplicated row names found")
    
    fuzzygroups <- dividedintosmallgroups(originalTimeSeriesData,
                                          minElements, 
                                          maxElements,
                                          membexp)
    clusters <- dissolve(findAllLeaves(fuzzygroups))
    
    clusters <- lapply(clusters, function(cluster) unique(rownames(cluster)))
    
    combined_clusters <- list()
    processed <- c()
    for (i in seq_along(clusters)) {
        processed <- c(processed, i)
        clusters2 <- clusters[-processed]
        for (j in seq_along(clusters2)) {
            if (identical(clusters[[i]], clusters2[[j]])) 
                (next)()
            
            newset <- unique(c(clusters[[i]], clusters2[[j]]))
            newset <- genes[which(genes %in% newset)]
            combined_clusters[[length(combined_clusters) + 1]] <- newset
        }
    }
    futile.logger::flog.info(
        sprintf("(clusterdDiscreteData) generates: %s combined groups", 
                length(combined_clusters)))
    res <- list(combined_clusters = combined_clusters,
                discretedTimeSeriesdata = discretedTimeSeriesdata,
                target_genes = genes)
    class(res) <- "ClusteredTimeseriesData"
    res
}

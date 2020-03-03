# The follow functions should be used internally and none of them should be exposed to outside user.
#'Executing a list of processes in parallel
#'
#'@param parallelFuc A process that will be run in parallel
#'@param listitems The main list for each process
#'@param ... other parameters
#'@return No return
#'@examples
#'coming later
#' @export
doParallelWork <- function(parallelFuc, listitems, ...) {
    if (!is.function(match.fun(parallelFuc))) {
        stop("The parameter is not a function")
    }
    if (is.null(listitems)) {
        stop("The type of listitems must be a list")
    }

    cl <- parallel::makeCluster(parallel::detectCores()[1] - 1, type = "SOCK")
    # * Do the work .....**
    res <- parallel::clusterApply(cl, 1:length(listitems), parallelFuc, listitems, ...)
    parallel::stopCluster(cl)
    res <- unlist(res, recursive = FALSE)
    # closeAllConnections()
    cond1 <- sapply(res, function(entry) !is.null(entry))
    res[cond1][unlist(lapply(res[cond1], length) != 0)]
    res
}

#'Executing a list of processes not in parallel
#'
#'@param parallelFuc A process that will be run in parallel
#'@param listitems The main list for each process
#'@param ... other parameters
#'@return No return
#'@examples
#'coming later
#' @export
doNonParallelWork <- function(parallelFuc, listitems, ...) {
    # foreach environment, the option stringsAsFactors is required to be set to FALSE
    options(stringsAsFactors = FALSE)
    if (!is.function(match.fun(parallelFuc))) {
        stop("The parameter is not a function")
    }
    if (is.null(listitems)) {
        stop("The type of listitems must be a list")
    }

    res <- lapply(1:length(listitems), parallelFuc, listitems, ...)
    res <- unlist(res, recursive = FALSE)
    cond1 <- sapply(res, function(entry) !is.null(entry))
    # remove the unwant outer list
    res[cond1][unlist(lapply(res[cond1], length) != 0)]
    res
}

#'Executing a list of processes in parallel with decreasing list
#'
#'@param parallelFuc A process that will be run in parallel
#'@param listitems The main list for each process
#'@param unprocessedListitems The remain list items that haven't been processed yet
#'@param ... other parameters
#'@return No return
#'@examples
#' coming later
#' @export
doNonParallelWorkDecrease <- function(parallelFuc, listitems, unprocessedListitems, ...) {
    # foreach environment, the option stringsAsFactors is required to be set to FALSE
    options(stringsAsFactors = FALSE)
    if (!is.function(match.fun(parallelFuc))) {
        stop("The parameter is not a function")
    }
    if (is.null(listitems)) {
        stop("The type of listitems must be a list")
    }
    
    res <- lapply(1:length(listitems), function(k,listitems, ...) {
        if (length(unprocessedListitems) > 0) {
            reindex <- which(listitems[k] %in% unprocessedListitems)
            unprocessedListitems2 <- unprocessedListitems[-reindex]
            unprocessedListitems <- unprocessedListitems2
            parallelFuc(k, listitems, unprocessedListitems2, ...)
        }
    }, listitems, ...)
    
    res <- unlist(res, recursive = FALSE)

    cond1 <- sapply(res, function(entry) !is.null(entry))
    # remove the unwant outer list
    res[cond1][unlist(lapply(res[cond1], length) != 0)]
}
unregister <- function() {
    env <- foreach:::.foreachGlobals
    rm(list = ls(name = env), pos = env)
}

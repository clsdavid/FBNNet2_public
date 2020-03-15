#~~~~~~~~1~~~~~~~~~2~~~~~~~~~3~~~~~~~~~4~~~~~~~~~5~~~~~~~~~6~~~~~~~~~7~~~~~~~~~8
#'A function to find all possible FBM (Fundamental Boolean models) attractors
#'
#'@param fbnnetwork An object of FBNNetwork
#'@param startStates A list of initial states, the row names of each state must
#' be matched with the genes
#'@param genes a list of genes which index order must match with the current
#' state
#'@param type A type of Boolean network update schema choosen from synchronous,
#' asynchronous and time step based
#'@param genesOn It is a vector of genes that are marked as On
#'@param genesOff It is a vector of genes that are marked as Off
#'@param maxSearch The maximum timesteps that the system will try to search.
#'@return Attractor objects
#' @author Leshi Chen, leshi, chen@lincolnuni.ac.nz, chenleshi@hotmail.com
#' @keywords Fundamental Boolean Network, Boolean Network, Genetic Regulatory
#'  Network
#'
#' @references Chen et al.(2018), Front. Physiol., 25 September 2018, 
#' (\href{https://doi.org/10.3389/fphys.2018.01328}{Front. Physiol.})
#' @references Mussel, Hopfensitz et al. 2010, 
#'  BoolNet - an R package for generation, reconstruction and analysis of
#'   Boolean networks
#' 
#'@examples
#' requires(BoolNet)
#' network<-loadNetwork('testthat/others/cellcycle.txt')
#' initialStates<-generateAllCombinationBinary(network$genes)
#' trainingseries<-genereateBoolNetTimeseries(network,
#'  initialStates,
#'  43,
#'   type='synchronous')
#' cube<-constructFBNCube(network$genes,
#' network$genes,
#' trainingseries,
#' 4,
#' 1,
#' TRUE)
#' NETWORK2<-mineFBNNetwork(cube,network$genes)
#' attractor<-searchForAttractors(NETWORK2,
#' initialStates,
#' network$genes)
#' print(attractor)
#' FBNNetwork.Graph.DrawAttractor(NETWORK2,attractor,2)
#' @export
searchForAttractors <- function(fbnnetwork,
                                startStates = list(), 
                                genes,
                                type = c("synchronous", "asynchronous"),
                                genesOn = c(), 
    genesOff = c(), maxSearch = 1000) {
    if (!(inherits(fbnnetwork, "FundamentalBooleanNetwork"))) 
        stop("Network must be inherited from FundamentalBooleanNetwork")
    
    type <- match.arg(type, c("synchronous", "asynchronous"))
    
    # random generated test data
    if (length(startStates) == 0) {
        startStates <- generateAllCombinationBinary(genes)
    }
    
    if (length(genesOn) > 0 & is.vector(genesOn)) {
        genesOnIndexes <- which(genes %in% genesOn)
        lapply(startStates, function(state) state[genesOnIndexes] <- 1)
    }
    
    if (length(genesOff) > 0 & is.vector(genesOff)) {
        genesOffIndexes <- which(genes %in% genesOff)
        lapply(startStates, function(state) state[genesOffIndexes] <- 0)
    }
    
    resultList <- list()
    basinStates <- list()
    stateAssighed <- list()
    tempbasinStates <- list()
    
    tryCatch({
        for (s in seq_along(startStates)) {
            iniState <- startStates[[s]]
            searchedStates <- list()
            names(iniState) <- genes
            
            numrow <- length(iniState)
            
            mat <- matrix(0,
                          nrow = numrow,
                          ncol = maxSearch,
                          byrow = FALSE,
                          dimnames = list(genes, sequence(maxSearch)))
            # get initial state
            mat[, 1] <- iniState
            premat <- mat[, 1]
            
            
            if (length(stateAssighed) > 0) {
                if (Position(function(x) identical(x,
                                                   iniState),
                             stateAssighed,
                             nomatch = 0) > 0) {
                  next
                }
            }
            
            found <- FALSE
            maxS <- 1
            searchedStates[[length(searchedStates) + 1]] <- iniState
            tempbasinStates[[length(tempbasinStates) + 1]] <- iniState
            decayIndex <- c()
            
            while (!found & maxS <= maxSearch) {
                k <- maxS + 1
                result <- getFBMSuccessor(fbnnetwork,
                                          premat,
                                          k,
                                          genes, 
                                          type = type,
                                          decayIndex)
                nextState <- result$nextState
                decayIndex <- result$decayIndex
                
                names(nextState) <- genes
                # check if an attractor is found
                if (length(stateAssighed) > 0) {
                  if (Position(function(x) identical(x,
                                                     nextState),
                               stateAssighed, nomatch = 0) > 0) {
                    break
                  }
                }
                # correction
                if (length(basinStates) > 0) {
                  foundI <- 0
                  lapply(seq_along(basinStates), function(statesIndex) {
                    states <- basinStates[[statesIndex]]
                    if (Position(function(x) identical(x, nextState),
                                 states, nomatch = 0) > 0) {
                      foundI <- statesIndex
                    }
                  })
                  if (length(tempbasinStates) > 0 & foundI > 0) {
                    basinStates[[foundI]] <- unique(
                        dissolve(list(dissolve(basinStates[[foundI]]),
                                 tempbasinStates)))
                    break
                  }
                }
                
                # find vector in list of vectors
                searchedIndex <- Position(function(x) identical(x, nextState),
                                          searchedStates, nomatch = 0)
                
                if (searchedIndex > 0) {
                  found <- TRUE
                  # check if the attractor has been found
                  attractors <- lapply(resultList, function(result) result[[1]])
                  if (!Position(function(x) identical(x, nextState), 
                                attractors, nomatch = 0) > 0) {
                    currentlength <- length(resultList) + 1
                    resultList[[currentlength]] <- list()
                    resultList[[currentlength]] <- dissolve(
                        list(dissolve(searchedStates[searchedIndex:length(searchedStates)]), 
                        nextState))
                    basinStates[[currentlength]] <- list()
                    temp <- dissolve(tempbasinStates)
                    basinStates[[currentlength]] <- temp[!temp %in% resultList[[currentlength]]]
                    names(basinStates)[[currentlength]] <- currentlength
                    stateAssighed <- unique(dissolve(list(stateAssighed, resultList[[currentlength]])))
                  }
                } else {
                  found <- FALSE
                  searchedStates[[length(searchedStates) + 1]] <- nextState
                  tempbasinStates[[length(tempbasinStates) + 1]] <- nextState
                }
                iniState <- nextState
                mat[, k] <- nextState
                premat <- mat[, 1:k]
                maxS <- maxS + 1
            }
        }
        res <- list()
        res[[1]] <- resultList
        names(res)[[1]] <- "Attractors"
        res[[2]] <- genes
        names(res)[[2]] <- "Genes"
        res[[3]] <- basinStates
        names(res)[[3]] <- "BasinOfAttractor"
        class(res) <- c("FBMAttractors")
        return(res)
    }, error = function(e) {
        mes <- e$message
        stop(sprintf("Error executing FBN model: %s", e$message))
    })
}

#'A method to fix a specific genes in the FBM networks
#'@param network A fundamentalBooleanNetwork type of Network
#'@param fixIndices, a vector of gene indexes that are intended
#' to be fixed
#'@param values, a vector of values for genes. 0 means fixing
#' the gene to inhibition 1 means fixing the gene to activation
#'  and -1 means no fixing
#'@noRd
networkFixUpdate <- function(network, fixIndices, values) {
    if (!(inherits(network, "FundamentalBooleanNetwork"))) 
        stop("Network must be inherited from FundamentalBooleanNetwork")
    
    if (length(fixIndices) != length(values) && length(values) != 1) 
        stop("fixIndices and values must have the same number of elements!")
    
    if (any(is.na(network$fixed[fixIndices]))) 
        stop("fixIndices contains invalid indices!")
    
    if (any(values != 0 & values != 1 & values != -1)) 
        stop("Please supply only 0, 1, or -1 in values!")
    
    network$fixed[fixIndices] <- as.integer(values)
    
    network
}

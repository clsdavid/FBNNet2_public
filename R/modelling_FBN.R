
#'An utility function to verify whether or not the input expression and the required gene input are marched
#'
#'@param geneState Pre gene input state
#'@param expression The expression of the target regulatory function
#'@return TRUE or FALSE
#'@examples
#' ##coming later
isSatisfied <- function(geneState, expression) {
    res <- TRUE
    # inistate<-geneState[[1]][[1]]
    if (identical(expression, "1")) {
        # over expressed?
        return(TRUE)
    }
    
    if (identical(expression, "0")) {
        # over inhibited
        return(FALSE)
    }
    
    for (i in seq_along(geneState)) {
        index <- which(expression == names(geneState)[[i]])[1]
        if (length(index) == 0 | length(index) > 1) {
            res <- FALSE
        }
        
        if (index > 1) {
            # if contain negation
            if (identical(expression[index - 1], "!")) {
                res <- res & !(as.numeric(geneState[[i]]) == 1L)
            } else {
                if (!identical(expression[index], "!") & !identical(expression[index], "&") & !identical(expression[index], 
                  "(") & !identical(expression[index], ")")) {
                  res <- res & (as.numeric(geneState[[i]]) == 1L)
                }
            }
        } else {
            if (!identical(expression[index], "!") & !identical(expression[index], "&") & !identical(expression[index], "(") & 
                !identical(expression[index], ")")) {
                res <- res & (as.numeric(geneState[[i]]) == 1L)
            }
        }
    }
    return(res)
}


#'A random selection function denoted as P[]
#'
#'@param probability A probability of an event
#'@return TRUE or FALSE indicating whether or not the event has been selected
#'@examples
#' ##coming later
#' @export
randomSelection <- function(probability) {
    sample(c(FALSE, TRUE), size = 1, replace = TRUE, prob = c(1 - probability, probability))
}

#'This method is used to reconstruct time series data
#'
#'@param fbnnetwork An object of FBNNetwork
#'@param initialStates A list of initial states
#'@param type Specify the type of the Fundamental Boolean model (synchronous or asynchronous)
#'@param maxTimepoints The max time points that are going to be constructed
#'@return A list object that contains reconstructed time series and FBN network
#'@examples
#' mat1<-matrix(c('1','0','0','1','0','0','0','1','1'),3,3, dimnames=list(c('gene1','gene2','gene3'),c('1','2','3')))
#' mat2<-matrix(c('1','1','0','1','0','1','1','1','0'),3,3, dimnames=list(c('gene1','gene2','gene3'),c('1','2','3')))
#' listtestOriginal<-list(mat1,mat2)
#' cube<-constructFBNCube(c('gene1','gene2'),c('gene1','gene2','gene3'),listtestOriginal,4,1,FALSE)
#' network<-mineFBNNetwork(cube,c('gene1','gene2'))
#' mat1<-matrix(c('1','1','0'),3,1, dimnames=list(c('gene1','gene2','gene3'),c('1')))
#' mat2<-matrix(c('1','1','1'),3,1, dimnames=list(c('gene1','gene2','gene3'),c('1')))
#' listtestInitial<-list(mat1,mat2)
#' result<-reconstructTimeseries(network,listtestInitial,1,'synchronous')
#' @export
reconstructTimeseries <- function(fbnnetwork, initialStates, type = c("synchronous", "asynchronous"), maxTimepoints = 100, 
    useParallel = FALSE) {
    
    
    if (!is.numeric(maxTimepoints) | maxTimepoints <= 0L | !all.equal(maxTimepoints, as.integer(maxTimepoints))) 
        stop("maxTimepoints must be integer")
    
    
    if (!(inherits(fbnnetwork, "FundamentalBooleanNetwork"))) 
        stop("Network must be inherited from FundamentalBooleanNetwork")
    
    type <- type[1]
    # code
    
    genes <- fbnnetwork$genes
    
    internal_fn <- function(i, p_initialStates, p_fbnnetwork, p_genes, p_type, p_maxTimepoints) {
        initialState <- p_initialStates[[i]]
        res <- list()
        res[[i]] <- transitionStates(initialState, p_fbnnetwork, p_genes, p_type, p_maxTimepoints)
        return(res)
    }
    
    if (useParallel) {
        reconstructed <- doParallelWork(internal_fn, initialStates, fbnnetwork, genes, type, maxTimepoints)
        
    } else {
        reconstructed <- doNonParallelWork(internal_fn, initialStates, fbnnetwork, genes, type, maxTimepoints)
    }
    
    # doNonParallelWork
    if (!is.null(reconstructed) & length(reconstructed) > 0) {
        # remove null entry
        cond1 <- sapply(reconstructed, function(entry) !is.null(entry))
        reconstructed <- (reconstructed[cond1][unlist(lapply(reconstructed[cond1], length) != 0)])
        class(reconstructed) <- c("FBNTimeSeries")
    }
    
    if (useParallel) {
        closeAllConnections()
    }
    reconstructed
}


# private method

transitionStates <- function(initialState, fbnNetwork, genes, type = c("synchronous", "asynchronous"), maxTimepoints) {
    
    numrow <- length(genes)
    rowNames <- genes
    colNames <- c(1:maxTimepoints)
    mat <- matrix(0, nrow = numrow, ncol = length(colNames), byrow = FALSE, dimnames = list(rowNames, colNames))
    
    vector <- unlist(initialState)
    # get initial state
    mat[, 1] <- vector
    k <- 2
    premat <- mat[, 1:k - 1]
    decayIndex <- c()
    timestepTrack <- list()
    while (k <= length(colNames)) {
        nextState <- getFBMSuccessor(fbnNetwork = fbnNetwork, previous_states = premat, current_step = k, genes = rowNames, 
            type = type, decayIndex = decayIndex)
        mat[, k] <- nextState$nextState
        premat <- mat[, 1:k]
        decayIndex <- nextState$decayIndex
        k <- k + 1
    }
    
    mat
}

#'This method is used to calculate the next state
#'
#'@param fbnnetwork An object of FBNNetwork
#'@param previous_states A vector of current gene state
#'@param genes a list of genes which index order must match with the current state
#'@param type A type of Boolean network update schema choosen from synchronous, asynchronous based. Asynchronous will randomly pick up a gene to process at time.
#'@param decayIndex An value indicates the period of time when to degrade an activated gene if no activators presented. It is usually one time step
#'@return A list object that contains reconstructed time series and FBN network
#'@examples
#' require(BoolNet)
#' require(FBNNet)
#' network<-loadNetwork('testthat/others/cellcycle.txt')
#' trainingseries<-FBNDataReduction(generateTimeSeries(network,2000,43))
#' cube<-constructFBNCube(network$genes,network$genes,trainingseries,4,1,TRUE)
#' NETWORK2<-mineFBNNetwork(cube,network$genes)
#' state<-c('0','1','1','0','1','1','1','1','0','0')
#' names(state)<-c('CycD','Rb','E2F','CycE','CycA','p27','Cdc20','Cdh1','UbcH10','CycB')
#' getFBMSuccessor(NETWORK2,state,names(state),'synchronous')
#' @export
getFBMSuccessor <- function(fbnNetwork, previous_states, current_step, genes, type = c("synchronous", "asynchronous"), decayIndex = c()) {
    internalFun <- function(gene, interactions, geneState, previous_states, current_step, timedecay, decayIndex = 1) {
        # gene, a target gene
        genefunctions <- interactions[[gene]]
        ini <- geneState
        decay <- timedecay
        if (decay < 1) 
            decay <- 1L
        if (length(genefunctions) > 0) {
            # find all activators' probabilities
            condOfActivation <- sapply(genefunctions, function(activator) activator$type == 1L)
            funcOfActivators <- genefunctions[condOfActivation]
            # find all inhibitors' probabilities
            condOfInhibitors <- sapply(genefunctions, function(inhibitor) inhibitor$type == 0L)
            funcOfInhibitors <- genefunctions[condOfInhibitors]
        } else {
            funcOfActivators <- list()
            funcOfInhibitors <- list()
        }
        
        prFA <- FALSE
        probabilityFA <- 0
        timeporal_applied <- FALSE
        if (length(funcOfActivators) > 0) {
            for (fa in seq_along(funcOfActivators)) {
                fbnName <- names(funcOfActivators)[[fa]]
                
                fa_timestep <- as.numeric(funcOfActivators[[fbnName]]$timestep)
                adapted_timestep <- current_step - fa_timestep
                if (ncol(previous_states) < adapted_timestep || adapted_timestep < 1) {
                  timeporal_applied <- TRUE
                  (next)()  #skip not enough timestep
                }
                
                pregeneInput <- dissolve(lapply(funcOfActivators[[fbnName]]$input, function(geneindex) {
                  res <- list()
                  res[[1]] <- previous_states[geneindex, adapted_timestep]
                  names(res)[[1]] <- genes[[geneindex]]
                  return(res)
                }))
                
                probability <- getProbabilityFromFunctionInput(1, funcOfActivators[[fbnName]]$expression, funcOfActivators[[fbnName]]$probability, 
                  pregeneInput)
                
                if (!length(probability) == 0) {
                  prFA <- prFA | randomSelection(probability)
                }
            }
        }
        
        prFD <- FALSE
        probabilityFD <- 0
        if (length(funcOfInhibitors) > 0) {
            for (fd in seq_along(funcOfInhibitors)) {
                fbnName <- names(funcOfInhibitors)[[fd]]
                fd_timestep <- as.numeric(funcOfInhibitors[[fbnName]]$timestep)
                adapted_timestep <- current_step - fd_timestep
                if (ncol(previous_states) < adapted_timestep || adapted_timestep < 1) {
                  timeporal_applied <- TRUE
                  (next)()  #skip not enough timestep
                }
                
                pregeneInput <- dissolve(lapply(funcOfInhibitors[[fbnName]]$input, function(geneindex) {
                  res <- list()
                  res[[1]] <- previous_states[geneindex, adapted_timestep]
                  names(res)[[1]] <- genes[[geneindex]]
                  return(res)
                }))
                
                probability <- getProbabilityFromFunctionInput(0, funcOfInhibitors[[fbnName]]$expression, funcOfInhibitors[[fbnName]]$probability, 
                  pregeneInput)
                
                if (!length(probability) == 0) {
                  prFD <- prFD | randomSelection(probability)
                }
            }
        }
        
        # nothing happen then deactivate when the decay time is due
        if (decay > 0 && !timeporal_applied) {
            if (prFA | prFD) {
                decayIndex <- 1L
            } else {
                # check decay
                if (as.numeric(decayIndex) >= as.numeric(decay)) {
                  ini <- FALSE
                  decayIndex <- 1L
                } else {
                  decayIndex <- as.numeric(decayIndex) + 1L
                }
            }
        }
        
        result <- (ini | prFA) & (!prFD)
        res <- list()
        res[[1]] <- result
        names(res)[[1]] <- "nextState"
        res[[2]] <- decayIndex
        names(res)[[2]] <- "decayIndex"
        # result<-(prFA)&(!prFD) nextState<-c(nextState,as.numeric(result))
        return(res)
    }
    
    
    
    # main body
    interactions <- fbnNetwork$interactions
    names(previous_states) <- genes
    if (!is.matrix(previous_states)) {
        previous_states <- t(t(previous_states))
        rownames(previous_states) <- genes
        colnames(previous_states) <- 1
    }
    
    
    # decayIndex<-0L
    nextState <- c()
    
    # fixed when the value fixed != -1 and fixed ==0
    fixedgeneIndex <- which(fbnNetwork$fixed != -1 || fbnNetwork$fixed == 0)
    fixedgenes <- fbnNetwork$genes[fixedgeneIndex]
    
    type <- match.arg(type, c("synchronous", "asynchronous"))
    
    if (length(decayIndex) == 0) {
        decayIndex <- rep(1L, length(genes))
        names(decayIndex) <- genes
    }
    
    ############################################################################################################# need to revise decayIndex[[gene]],currentState,timestepTrack to ensure it should work as expect randomly pick up a gene
    ############################################################################################################# to process at time
    if (type == "asynchronous") {
        nonfixedgenes <- genes[!genes %in% fixedgenes]  #if gene is fixed then get the current gene stat
        # randomly pickup one
        set.seed(100)
        gene <- sample(nonfixedgenes, 1)
        timedecay <- fbnNetwork$timedecay[[gene]]
        nextState <- previous_states[, current_step - 1]
        ini <- previous_states[gene, current_step - 1] == 1L
        result <- internalFun(gene, interactions, ini, previous_states, current_step, timedecay, decayIndex[[gene]])
        names(nextState) <- genes
        names(decayIndex) <- genes
        
        decayIndex[[gene]] <- result[[2]]
        # update next state for the random choosen gene
        nextState[[gene]] <- as.numeric(result[[1]])
        
    } else {
        s_res <- lapply(seq_along(genes), function(j) {
            gene <- genes[[j]]
            ini <- previous_states[gene, current_step - 1] == 1L
            timedecay <- fbnNetwork$timedecay[[gene]]
            # if gene is fixed then get the current gene state
            if (gene %in% fixedgenes) {
                return(ini, 0)
            }
            result <- internalFun(gene, interactions, ini, previous_states, current_step, timedecay, decayIndex[[gene]])
            return(c(as.numeric(result[[1]]), result[[2]]))
        })
        nextState <- sapply(s_res, function(entry) entry[[1]])
        decayIndex <- sapply(s_res, function(entry) entry[[2]])
        names(nextState) <- genes
        names(decayIndex) <- genes
    }
    res <- list()
    res[[1]] <- nextState
    names(res)[[1]] <- "nextState"
    res[[2]] <- decayIndex
    names(res)[[2]] <- "decayIndex"
    return(res)
}


#'An internal method to get the probability of regulatory function from the cube
#'
#'@param FBNCube An object of orchard cube
#'@param targetgene The target gene
#'@param funcType type of function
#'@param FBNExpression expression of function
#'@param preGeneInputs pre input genes' state
#'@return A probablity of the target regulatory function
#'@examples
#' ##coming later
getProbabilityFromFunctionInput <- function(funcType, FBNExpression, FBNProbability, preGeneInputs) {
    # internal functions
    isInputStateMatchedFBNFunction <- function(inputstate, expression) {
        splitedexpression <- splitExpression(expression, 1, FALSE)
        isSatisfied(preGeneInputs, splitedexpression)
    }
    
    if (!funcType %in% c(1, 0)) {
        stop("The function type value must be in the range of \"0\" and \"1\"")
    }
    
    if (isInputStateMatchedFBNFunction(preGeneInputs, FBNExpression)) {
        probability <- as.numeric(FBNProbability)
    } else {
        probability <- 0
    }
    
    as.numeric(probability)
}

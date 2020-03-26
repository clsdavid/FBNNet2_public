
#'Mine FBN Networks from an Orchard cube
#'
#'@param fbnGeneCube A pre constructed Orchard cube
#'@param genes The target genes in the output
#'@param useParallel An option turns on parallel
#'@param threshold_confidence A threshod of confidence (between 0 and 1)
#' that used to filter the Fundamental Boolean functions
#'@param threshold_error A threshod of error rate (between 0 and 1) that
#' used to filter the Fundamental Boolean functions
#'@param threshold_support A threshod of support (between 0 and 1) that
#' used to filter the Fundamental Boolean functions
#'@param maxFBNRules The maximum rules per type (Activation and Inhibition)
#' per gene can be mined or filtered, the rest will be discarded
#'@return A object of FBN network
#'@author Leshi Chen, leshi, chen@lincolnuni.ac.nz
#'@keywords Fundamental Boolean Network, Boolean Network, Genetic Regulatory Network
#'
#'@references Chen et al.(2018), Front. Physiol., 25 September 2018, 
#'(\href{https://doi.org/10.3389/fphys.2018.01328}{Front. Physiol.})
#'@examples
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
#' NETWORK <- mineFBNNetwork(cube,ExampleNetwork$genes)
#' NETWORK
#' @export
mineFBNNetwork <- function(fbnGeneCube, genes = NULL, useParallel = FALSE, threshold_confidence = 1, threshold_error = 0, threshold_support = 1e-05, 
  maxFBNRules = 5) {
  futile.logger::flog.info(sprintf("Enter mineFBNNetwork zone: genes=%s, useParallel=%s", paste(genes, sep = ", ", collapse = ", "), useParallel))
  checkProbabilityTypeData(threshold_confidence)
  checkProbabilityTypeData(threshold_error)
  checkProbabilityTypeData(threshold_support)
  checkNumeric(maxFBNRules)
  
  if (length(fbnGeneCube) == 0) 
    return(NULL)
  
  if (is.null(genes)) {
    genes <- names(fbnGeneCube)
  }
  time1 <- as.numeric(Sys.time())
  
  genes <- genes[genes %in% names(fbnGeneCube)]
  
  midle_result <- searchFBNNetworkCore(fbnGeneCube, genes, useParallel, threshold_confidence, threshold_error, threshold_support, maxFBNRules)
  finalresult <- mineFBNNetworkWithCores(midle_result, genes, threshold_error, maxFBNRules)
  time2 <- as.numeric(Sys.time())
  print(paste("Total cost ", time2 - time1, " seconds to mine Fundamental Boolean Functions ", sep = "", collapse = ""))
  futile.logger::flog.info(sprintf("Leave mineFBNNetwork zone"))
  finalresult
}

#' An internal function
#' 
#' @param searchFBNNetworkCore A result of \code{searchFBNNetworkCore}
#' @param genes Genes that involved.
#' @param threshold_error A threshod of error rate (between 0 and 1) 
#' that used to filter the Fundamental Boolean functions
#' @param maxFBNRules The maximum rules per type (Activation and Inhibition) 
#' per gene can be mined or filtered, the rest will be discarded
mineFBNNetworkWithCores <- function(searchFBNNetworkCore, genes = NULL, threshold_error, maxFBNRules) {
  if (is.null(genes)) {
    genes <- names(searchFBNNetworkCore)
  }
  
  futile.logger::flog.info(sprintf("Enter mineFBNNetworkWithCores zone"))
  midle_result <- mineFBNNetworkStage2(searchFBNNetworkCore, threshold_error, maxFBNRules)
  finalresult <- convertMinedResultToFBNNetwork(midle_result, genes)
  
  if (length(finalresult) > 0) {
    cond1 <- sapply(finalresult, function(entry) !is.null(entry))
    if (length(cond1) > 0) {
      finalresult <- (finalresult[cond1][unlist(lapply(finalresult[cond1], length) != 0)])
    } else {
      stop("No Network generated")
    }
  } else {
    stop("No Network generated")
  }
  
  finalresult <- filterNetworkConnections(finalresult)
  
  class(finalresult) <- c("FundamentalBooleanNetwork", "BooleanNetworkCollection")
  futile.logger::flog.info(sprintf("Leave mineFBNNetworkWithCores zone"))
  finalresult
}

#' Internal method
#' 
#' @param factors The value that needs to be duplicated.
#' 
removeDuplicates <- function(factors) {
  cond1 <- sapply(factors, function(x) x[[4]])
  return(factors[!duplicated(cond1)])
}


#' Internal method
#'@param fbnGeneCube A pre constructed Orchard cube
#'@param genes The target genes in the output
#'@param useParallel An option turns on parallel
#'@param threshold_confidence A threshod of confidence (between 0 and 1)
#' that used to filter the Fundamental Boolean functions
#'@param threshold_error A threshod of error rate (between 0 and 1) that
#' used to filter the Fundamental Boolean functions
#'@param threshold_support A threshod of support (between 0 and 1) that
#' used to filter the Fundamental Boolean functions
#'@param maxFBNRules The maximum rules per type (Activation and Inhibition)
#' per gene can be mined or filtered, the rest will be discarded
searchFBNNetworkCore <- function(fbnGeneCube, genes, useParallel = FALSE, threshold_confidence = 1, threshold_error = 0, threshold_support = 1e-04, 
  maxFBNRules = 5) {
  
  if (useParallel) {
    useParallel = FALSE
    futile.logger::flog.info(sprintf("The parallel for network is not support"))
  }
  futile.logger::flog.info(sprintf("Enter searchFBNNetworkCore zone: useParallel=%s", useParallel))
  
  
  network_configs <- new.env(hash = TRUE)
  network_configs$threshold_confidence <- threshold_confidence
  network_configs$threshold_error <- threshold_error
  network_configs$threshold_support <- threshold_support
  network_configs$maxFBNRules <- maxFBNRules
  ## main recursive function to search for rults (fruits)
  internalloop <- function(i, genes, env) {
    
    recursiveMiningFBNFunction <- function(fbnGeneCubeStem, targetgene, currentsubgene, groupbyGene, network_configs, findTrue, findFalse) {
      threshold_confidence <- network_configs$threshold_confidence
      threshold_error <- network_configs$threshold_error
      threshold_support <- network_configs$threshold_support
      
      if (is.null(fbnGeneCubeStem[[currentsubgene]])) {
        return(NULL)
      }
      
      if (is.null(groupbyGene)) {
        groupbyGene <- currentsubgene
      }
      # identity<-paste(input,'_',sep='',collapse = '')
      res <- list()
      resT <- list()
      resF <- list()
      errorActivator <- 0
      errorInhibitor <- 0
      pickT_causality_test <- 0
      pickF_causality_test <- 0
      pickTvalue <- 0
      pickFvalue <- 0
      pickTallconfidence <- 0
      pickFallconfidence <- 0
      # defauilt time step is 1 (synchronous update schema)
      timestepT <- 1
      timestepF <- 1
      bestFitP <- 0
      bestFitN <- 0
      
      # get measures at the current level
      input <- sort(fbnGeneCubeStem[[currentsubgene]][["Input"]])
      pickT <- fbnGeneCubeStem[[currentsubgene]][["ActivatorAndInhibitor"]][["Activator"]]
      pickF <- fbnGeneCubeStem[[currentsubgene]][["ActivatorAndInhibitor"]][["Inhibitor"]]
      
      
      if (!is.null(pickT)) {
        # pvalue<-pickT['chiSQ']
        pickTvalue <- as.numeric(pickT["Confidence"])  #avgSignal_T + sdSignal_T#
        pickTallconfidence <- as.numeric(pickT["all_confidence"])
        pickTmaxconfidence <- as.numeric(pickT["max_confidence"])
        ## pickTMutualInfo <- abs(as.numeric(pickT['MutualInfo']))
        pickTsupport <- as.numeric(pickT["support"])
        pickT_causality_test <- as.numeric(pickT["causality_test"])  #*******
        errorActivator <- as.numeric(pickT["Noise"])
        identityT <- pickT["Identity"]
        pickTType <- pickT["type"]
        timestepT <- as.numeric(pickT["timestep"])
        bestFitP <- as.numeric(pickT["bestFitP"])
      }
      
      if (!is.null(pickF)) {
        # pvalue<-pickF['chiSQ']
        pickFvalue <- as.numeric(pickF["Confidence"])  #avgSignal_F + sdSignal_F#
        pickFallconfidence <- as.numeric(pickF["all_confidence"])
        pickFmaxconfidence <- as.numeric(pickF["max_confidence"])
        ## pickFMutualInfo <- abs(as.numeric(pickF['MutualInfo']))
        pickFsupport <- as.numeric(pickF["support"])
        pickF_causality_test <- as.numeric(pickF["causality_test"])  #******
        errorInhibitor <- as.numeric(pickF["Noise"])
        identityF <- pickF["Identity"]
        pickFType <- pickF["type"]
        timestepF <- as.numeric(pickF["timestep"])
        bestFitN <- as.numeric(pickF["bestFitN"])
      }
      
      # pickTcosine
      if (!is.null(pickT) && pickTsupport >= threshold_support && pickT_causality_test >= 1 && pickTvalue >= threshold_confidence) {
        getT <- pickT
      } else {
        getT <- NULL
      }
      
      if (!is.null(pickF) && pickFsupport >= threshold_support && pickF_causality_test >= 1 && pickFvalue >= threshold_confidence) {
        getF <- pickF
      } else {
        getF <- NULL
      }
      
      
      newindex <- length(res) + 1
      if (!is.null(getT)) {
        
        res[[newindex]] <- c(targets = targetgene, factor = getT[1], type = 1, identity = identityT, error = round(errorActivator, 5), P = round(pickTvalue, 
          4), support = pickTsupport, timestep = timestepT, input = paste(input, collapse = ","), numOfInput = length(input), causality_test = round(pickT_causality_test, 
          5), GroupBy = groupbyGene, all_confidence = pickTallconfidence, dimensionType = pickTType, bestFitP = bestFitP)
        
        newindex <- length(res) + 1
        
      }
      
      if (!is.null(getF)) {
        # identityF<-paste(identity,0,sep='',collapse = '')
        res[[newindex]] <- c(targets = targetgene, factor = getF[1], type = 0, identity = identityF, error = round(errorInhibitor, 5), P = round(pickFvalue, 
          4), support = pickFsupport, timestep = timestepF, input = paste(input, collapse = ","), numOfInput = length(input), causality_test = round(pickF_causality_test, 
          5), GroupBy = groupbyGene, all_confidence = pickFallconfidence, dimensionType = pickFType, bestFitN = bestFitN)
        
        newindex <- length(res) + 1
      }
      
      
      # go through sub levels
      findTrue <- findTrue || !is.null(getT)
      findFalse <- findFalse || !is.null(getF)
      
      if (!(findTrue && findFalse)) {
        if (!is.null(fbnGeneCubeStem[[currentsubgene]]$SubGenesT)) {
          fbnGeneCubeStemT <- fbnGeneCubeStem[[currentsubgene]]$SubGenesT
          nextgenesT <- names(fbnGeneCubeStemT)
          indexT <- 1
          temp_resTU <- lapply(seq_along(nextgenesT), function(j) {
          nextgene <- nextgenesT[[j]]
          dissolve(recursiveMiningFBNFunction(fbnGeneCubeStem = fbnGeneCubeStemT, targetgene = targetgene, currentsubgene = nextgene, 
            groupbyGene = groupbyGene, network_configs = network_configs, findTrue = findTrue, findFalse = findFalse))
          })
          
          for (j in seq_along(temp_resTU)) {
          nextgene <- nextgenesT[[j]]
          resTU <- temp_resTU[[j]]
          if (length(resTU) > 0) {
            condT <- sapply(resTU, function(entry) !is.null(entry))
            resT[[indexT]] <- resTU[unlist(lapply(resTU[condT], length) != 0)]
            indexT <- indexT + 1
          }
          }
          if (length(resT) > 0) {
          resT <- dissolve(resT)
          }
        }
        if (!is.null(fbnGeneCubeStem[[currentsubgene]]$SubGenesF)) {
          fbnGeneCubeStemF <- fbnGeneCubeStem[[currentsubgene]]$SubGenesF
          nextgenesF <- names(fbnGeneCubeStemF)
          indexF <- 1
          temp_resFU <- lapply(seq_along(nextgenesF), function(j) {
          nextgene <- nextgenesF[[j]]
          dissolve(recursiveMiningFBNFunction(fbnGeneCubeStem = fbnGeneCubeStemF, targetgene = targetgene, currentsubgene = nextgene, 
            groupbyGene = groupbyGene, network_configs = network_configs, findTrue = findTrue, findFalse = findFalse))
          })
          for (j in seq_along(nextgenesF)) {
          nextgene <- nextgenesF[[j]]
          resFU <- temp_resFU[[j]]
          if (length(resFU) > 0) {
            condF <- sapply(resFU, function(entry) !is.null(entry))
            resF[[indexF]] <- resFU[unlist(lapply(resFU[condF], length) != 0)]
            indexF <- indexF + 1
          }
          }
          if (length(resF) > 0) {
          resF <- dissolve(resF)
          }
        }
      }
      
      if (length(resT) > 0) {
        res <- list(res, resT)
      }
      
      if (length(resF) > 0) {
        res <- list(res, resF)
      }
      
      res <- dissolve(res)
      
      if (length(res) > 0) {
        cond1 <- sapply(res, function(entry) !is.null(entry))
        res <- (res[cond1][unlist(lapply(res[cond1], length) != 0)])
      }
      
      
      res
    }
    
    ## entry part
    network_configs <- env$configs
    
    targetGene <- genes[[i]]
    res <- list()
    res[[i]] <- list()
    names(res)[[i]] <- targetGene
    
    if (length(env$cube) == 0) 
      return(NULL)
    
    if (is.null(env$cube[[targetGene]])) 
      return(NULL)
    
    if (is.null(env$cube[[targetGene]]$SubGenes)) 
      return(NULL)
    
    currentStem <- env$cube[[targetGene]]$SubGenes
    nextgenes <- names(currentStem)
    
    temp_res <- lapply(seq_along(nextgenes), function(k) {
      currentGene <- nextgenes[[k]]
      dissolve(recursiveMiningFBNFunction(fbnGeneCubeStem = currentStem, targetgene = targetGene, currentsubgene = currentGene, groupbyGene = NULL, 
        network_configs = network_configs, findTrue = FALSE, findFalse = FALSE))
    })
    rm(list = "currentStem")
    
    for (k in seq_along(temp_res)) {
      currentGene <- nextgenes[[k]]
      subres <- temp_res[[k]]
      if (length(subres) > 0) {
        cond1 <- sapply(subres, function(entry) !is.null(entry))
        subres <- (subres[cond1][unlist(lapply(subres[cond1], length) != 0)])
        if (!is.null(subres)) {
          if (length(subres) > 0) {
          index <- length(res[[i]]) + 1
          # get result and remove duplicates
          resultsub <- dissolve(subres)
          res[[i]][[index]] <- resultsub
          names(res[[i]])[[index]] <- currentGene
          }
        }
      }
    }
    
    preresponse <- removeDuplicates(dissolve(res[[i]]))
    
    if (length(preresponse) == 0) {
      return(list())
    }
    
    res[[i]] <- preresponse
    res
  }
  
  res <- list()
  
  if (length(genes) == 0) {
    return(list())
  }
  
  if (is.null(fbnGeneCube)) {
    return(list())
  }
  
  mainParameters <- new.env(hash = TRUE, parent = globalenv())
  mainParameters$cube <- fbnGeneCube
  mainParameters$configs <- network_configs
  
  if (useParallel) {
    res <- doParallelWork(internalloop, genes, mainParameters)
    
  } else {
    res <- doNonParallelWork(internalloop, genes, mainParameters)
  }
  
  cond1 <- sapply(res, function(entry) !is.null(entry))
  res <- (res[cond1][unlist(lapply(res[cond1], length) != 0)])
  
  
  
  if (useParallel) {
    closeAllConnections()
  }
  
  futile.logger::flog.info(sprintf("Leave searchFBNNetworkCore zone"))
  rm(network_configs)
  res
  ## wrap result
}


#' Internal method
#' @param res The stage 1 result
#' @param threshold_error The error threshold
#' @param maxFBNRules The maximum FBN rules.
mineFBNNetworkStage2 <- function(res, threshold_error = 0, maxFBNRules = 5) {
  futile.logger::flog.info(sprintf("Enter mineFBNNetworkStage2 zone"))
  
  if (is.null(res) | length(res) == 0) {
    return(list())
  }
  
  finalFilteredlist <- list()
  filteredres <- list()
  targetgenes <- names(res)
  for (i in seq_along(targetgenes)) {
    target <- targetgenes[i]
    finalFilteredlist[[i]] <- list()
    names(finalFilteredlist)[[i]] <- target
    processed <- c()
    ruleset <- res[[target]]
    
    for (j in seq_along(ruleset)) {
      rule <- ruleset[[j]]
      processed <- c(j)
      ruleset2 <- ruleset[-processed]
      for (k in seq_along(ruleset2)) {
        rule2 <- ruleset2[[k]]
        if (as.numeric(rule[["numOfInput"]]) < as.numeric(rule2[["numOfInput"]]) && as.numeric(rule[["type"]]) == as.numeric(rule2[["type"]]) && 
          as.numeric(rule[["timestep"]]) == as.numeric(rule2[["timestep"]]) && all(splitExpression(rule[["input"]], 2, FALSE) %in% splitExpression(rule2[["input"]], 
          2, FALSE) == TRUE)) {
          finalFilteredlist[[i]][[length(finalFilteredlist[[i]]) + 1]] <- rule2
        }
      }
    }
    filteredres[[target]] <- res[[target]][!res[[target]] %in% unique(finalFilteredlist[[i]])]
  }
  
  cond1 <- sapply(filteredres, function(entry) !is.null(entry))
  filteredres <- (filteredres[cond1][unlist(lapply(filteredres[cond1], length) != 0)])
  
  filteredres <- lapply(filteredres, function(entry) {
    if (is.null(entry) | length(entry) == 0) {
      return(NULL())
    }
    activators <- list()
    inhibitors <- list()
    
    for (e in seq_along(entry)) {
      if (as.numeric(entry[[e]]["error"]) <= threshold_error && as.numeric(entry[[e]]["type"]) == 1) {
        activators[[length(activators) + 1]] <- c(entry[[e]], threshold_error = threshold_error)
      }
      
      if (as.numeric(entry[[e]]["error"]) <= threshold_error && as.numeric(entry[[e]]["type"]) == 0) {
        inhibitors[[length(inhibitors) + 1]] <- c(entry[[e]], threshold_error = threshold_error)
      }
    }
    
    cond1 <- sapply(activators, function(entry) !is.null(entry))
    if (length(cond1) > 0) {
      filteredActivators <- lapply(activators[cond1], length)
      activators <- (activators[cond1][unlist(filteredActivators != 0)])
    }
    
    cond1 <- sapply(inhibitors, function(entry) !is.null(entry))
    if (length(cond1) > 0) {
      filteredInhibits <- lapply(inhibitors[cond1], length)
      inhibitors <- (inhibitors[cond1][unlist(filteredInhibits != 0)])
    }
    
    activators <- activators[order(as.numeric(sapply(activators, "[[", "timestep")), as.numeric(sapply(activators, "[[", "error")), as.numeric(sapply(activators, 
      "[[", "numOfInput")), -as.numeric(sapply(activators, "[[", "support")))]
    if (length(activators) > maxFBNRules) {
      activators <- activators[1:maxFBNRules]
    }
    
    inhibitors <- inhibitors[order(as.numeric(sapply(inhibitors, "[[", "timestep")), as.numeric(sapply(inhibitors, "[[", "error")), as.numeric(sapply(inhibitors, 
      "[[", "numOfInput")), -as.numeric(sapply(inhibitors, "[[", "support")))]
    if (length(inhibitors) > maxFBNRules) {
      inhibitors <- inhibitors[1:maxFBNRules]
    }
    
    finalRes <- list()
    finalRes[[1]] <- activators
    finalRes[[2]] <- inhibitors
    finalRes <- dissolve(finalRes)
    # futher remove based on maximum activators and maximum inhibitors
    return(finalRes)
  })
  cond1 <- sapply(filteredres, function(entry) !is.null(entry))
  filteredres <- (filteredres[cond1][unlist(lapply(filteredres[cond1], length) != 0)])
  
  class(filteredres) <- c("FBNTrueCubeMiner")
  futile.logger::flog.info(sprintf("Leave mineFBNNetworkStage2 zone"))
  filteredres
}


#' An internal function to check if an expression is atomic node
#' @param subExpression A sub Boolean expression
#' 
isAtomNode <- function(subExpression) {
  unwantList <- c("|", "&", "(", "[", "{", ")", "]", "}")
  res <- TRUE
  for (i in seq_along(subExpression)) {
    if (subExpression[[i]] %in% unwantList) {
      res <- FALSE
      break
    }
  }
  return(res)
}

#' An internal function to check if an expression is a sub expression
#' @param subExpression A sub Boolean expression
#' 
isSubExpression <- function(subExpression) {
  groupA <- c("(", "[", "{")
  groupB <- c(")", "]", "}")
  startBracket <- 0L
  endBracket <- 0L
  numOfBracket <- 0L
  res <- FALSE
  for (i in seq_along(subExpression)) {
    if (subExpression[[i]] %in% groupA) {
      if (numOfBracket == 0L) {
        startBracket <- i
      }
      numOfBracket <- numOfBracket + 1L
      next
    }
    
    if (subExpression[[i]] %in% groupB) {
      if (numOfBracket == 1L) {
        endBracket <- i
      }
      numOfBracket <- numOfBracket - 1L
    }
    
  }
  
  if (startBracket == 1L && endBracket == length(subExpression)) {
    res <- TRUE
  }
  return(res)
}

#' An internal function to remove the first brasket
#' @param splitedExpression A splited Boolean expression
#' 
removeFistLevelBracket <- function(splitedExpression) {
  res <- splitedExpression
  if (isSubExpression(splitedExpression)) {
    res <- res[-1]
    res <- res[-length(res)]
  }
  
  return(res)
}

#' An internal function to check a splited Boolean expression
#' is applied DeMorganLaw
#' 
#' @param splitedExpression A splited Boolean expression
#' 
isAppliedDeMorganLaw <- function(splitedExpression) {
  groupA <- c("(", "[", "{")
  groupB <- c(")", "]", "}")
  startBracket <- 0L
  endBracket <- 0L
  numOfBracket <- 0L
  res <- FALSE
  subExpression <- splitedExpression
  
  if (mode(splitedExpression) == "character" && identical(splitedExpression[[1]], "!")) {
    subExpression <- subExpression[-1]
  } else {
    return(FALSE)
  }
  
  for (i in seq_along(subExpression)) {
    if (subExpression[[i]] %in% groupA) {
      if (numOfBracket == 0L) {
        startBracket <- i
      }
      numOfBracket <- numOfBracket + 1L
      next
    }
    
    if (subExpression[[i]] %in% groupB) {
      if (numOfBracket == 1L) {
        endBracket <- i
      }
      numOfBracket <- numOfBracket - 1L
    }
    
  }
  
  if (startBracket == 1L && endBracket == length(subExpression)) {
    res <- TRUE
  }
  return(res)
}

#' An internal function to flat a splited Boolean expression
#' with DeMorganLaw
#' 
#' @param splitedExpression A splited Boolean expression
#'
flatDeMorganLaw <- function(splitedExpression) {
  if (!isAppliedDeMorganLaw(splitedExpression)) 
    return(splitedExpression)
  negation <- splitedExpression[1]
  partOfExpression <- splitedExpression[-1]
  splitexp <- removeFistLevelBracket(partOfExpression)
  res <- convertIntoExpressionTree(splitexp)
  if (is.element("&", splitexp) && !is.element("|", splitexp)) {
    for (i in seq_along(res)) {
      if (identical(res[[i]], "&")) {
        res[[i]] <- c("|")
      } else {
        res[[i]] <- c(negation, res[[i]])
      }
    }
    
  } else if (!is.element("&", splitexp) && is.element("|", splitexp)) {
    for (i in seq_along(res)) {
      if (identical(res[[i]], "|")) {
        res[[i]] <- c("&")
      } else {
        res[[i]] <- c(negation, res[[i]])
      }
    }
  }
  
  return(res)
}

#' An internal function to convert a splited Boolean expression
#' into an expression tree.
#' 
#' @param splitedExpression A splited Boolean expression
#'
convertIntoExpressionTree <- function(splitedExpression) {
  splitexp <- removeFistLevelBracket(splitedExpression)
  Oprators <- c("|", "&")
  
  groupA <- c("(", "[", "{")
  groupB <- c(")", "]", "}")
  
  containedBracket <- 0
  res <- list()
  leftCut <- 1L
  if (length(splitexp) == 1L) {
    res[[length(res) + 1]] <- splitexp[[1]]
    return(res)
  }
  
  if (!is.element("&", splitexp) && !is.element("|", splitexp)) {
    res <- splitexp
    return(res)
  }
  
  for (i in seq_along(splitexp)) {
    if (splitexp[[i]] %in% groupA) {
      containedBracket <- containedBracket + 1
      next
    }
    
    if (splitexp[[i]] %in% groupB) {
      containedBracket <- containedBracket - 1
    }
    
    if (containedBracket > 0) {
      next
    }
    
    if (splitexp[[i]] %in% Oprators) {
      
      expressionLeft <- splitexp[leftCut:(i - 1)]
      expressionRight <- splitexp[(i + 1):length(splitexp)]
      if (isSubExpression(expressionLeft)) {
        res[[length(res) + 1]] <- convertIntoExpressionTree(expressionLeft)
      } else {
        res[[length(res) + 1]] <- expressionLeft
      }
      
      res[[length(res) + 1]] <- splitexp[[i]]
      
      if (isSubExpression(expressionRight)) {
        res[[length(res) + 1]] <- convertIntoExpressionTree(expressionRight)
      } else {
        # need to split demorganlaw into !(A*B)=!A || !B or !(A||B)=!A * !B
        if (isAtomNode(expressionRight)) {
          res[[length(res) + 1]] <- expressionRight
          break
        }
        
        if (isAppliedDeMorganLaw(expressionRight)) {
          res[[length(res) + 1]] <- flatDeMorganLaw(expressionRight)
          break
        }
        leftCut <- i + 1
      }
    }
  }
  
  return(res)
}
#' An internal function to construct the FBN functions
#' with an expression tree.
#' 
#' @param expressionTree A expression tree
#'
constructFBNFunctions <- function(expressionTree) {
  stem <- expressionTree
  operators <- c("|", "&")
  res <- list()
  index <- 1L
  newLine <- FALSE
  for (i in seq_along(stem)) {
    if (!(length(stem[[i]]) == 1L && stem[[i]] %in% operators)) {
      if (mode(stem[[i]]) == "list") {
        subRes <- constructFBNFunctions(stem[[i]])
        
        if (length(res) == index) {
          pre <- res[[index]]
        } else {
          pre <- c(NA)
        }
        
        for (j in seq_along(subRes)) {
          if (newLine) {
          res[[index]] <- subRes[[j]]
          } else {
          if (!is.na(pre)[[1]]) {
            res[[index]] <- c(pre, subRes[[j]])
          } else {
            res[[index]] <- subRes[[j]]
          }
          }
          if (j < length(subRes)) {
          index <- index + 1
          }
        }
        
      } else {
        if (index == (length(res) + 1)) {
          res[[index]] <- stem[[i]]
        } else {
          if (length(res) > 1) {
          inx <- 1
          while (inx <= length(res)) {
            res[[inx]] <- c(res[[inx]], stem[[i]])
            inx <- inx + 1
          }
          
          } else {
          res[[index]] <- c(res[[index]], stem[[i]])
          }
        }
      }
    } else {
      if (length(stem[[i]]) == 1L && stem[[i]] == "|") {
        index <- index + 1L
        newLine <- TRUE
      } else if (length(stem[[i]]) == 1L && stem[[i]] == "&") {
        if (length(res) > 1) {
          inx <- 1
          while (inx <= length(res)) {
          res[[inx]] <- c(res[[inx]], "&")
          inx <- inx + 1
          }
          
        } else {
          res[[index]] <- c(res[[index]], "&")
        }
        newLine <- FALSE
      }
    }
    
  }
  return(res)
}

#'generate FBN interaction
#'
#'@param expressionString An expression string
#'@param genes The involved genes contains in the expression string
#'@return A type object of interaction
#'@examples
#' ##coming later
#' @export
generateFBNInteraction <- function(expressionString, genes) {
  res <- list()
  splitedexpression <- splitExpression(expressionString, 1, FALSE)
  geneinputs <- which(genes %in% splitedexpression)
  res$input <- geneinputs
  res$expression <- expressionString
  return(res)
}

#' An internal function to reconstrunct FBN functions
#' 
#' @param name The name of the interaction
#' @param expressionstring The string presentation of an expression.
#' @param genes The involve genes
#' @param error The error value
#' @param type The type of the interaction
#' @param probability The probability of the interaction
#' @param support The support threshold.
#' @param timestep The time step of the interaction to take effectiveness on. 
regenerateInteractions <- function(name, expressionstring, genes, error, type, probability = NA, support = NA, timestep = 1) {
  res <- list()
  if (mode(genes) == "list") {
    genelist <- unlist(genes)
  } else {
    genelist <- genes
  }
  split <- splitExpression(expressionstring, 1, FALSE)
  tree <- convertIntoExpressionTree(split)
  functions <- constructFBNFunctions(tree)
  for (i in seq_along(functions)) {
    thisfuc <- functions[[i]]
    fuc <- generateFBNInteraction(paste(unlist(thisfuc), collapse = ""), genelist)
    res[[i]] <- list()
    
    # gene input
    res[[i]][[1]] <- fuc$input
    names(res[[i]])[[1]] <- "input"
    
    # expression
    res[[i]][[2]] <- fuc$expression
    names(res[[i]])[[2]] <- "expression"
    
    # error need to conside when fuc$input==0 and expression =1 or 0
    if (is.element("error", names(fuc))) {
      res[[i]][[3]] <- fuc$error
    } else {
      if (is.null(error)) {
        res[[i]][[3]] <- "NA"
      } else {
        res[[i]][[3]] <- error
      }
      
    }
    names(res[[i]])[[3]] <- "error"
    
    # function type
    if (!is.null(type)) {
      if (!is.numeric(type) | (type != 1L & type != 0L)) {
        stop("The value of function type must be 0 (Inhibition) or 1 (Activation)")
      }
      res[[i]][[4]] <- type
    } else {
      res[[i]][[4]] <- 1L  #Activation by default
    }
    
    names(res[[i]])[[4]] <- "type"
    
    if (res[[i]][[4]] == 1) 
      typename <- "Activator" else typename <- "Inhibitor"
    
    
    # probability / confidence measure
    if (is.element("probability", names(fuc))) {
      res[[i]][[5]] <- fuc$probability
    } else {
      if (is.na(probability) | is.null(probability)) {
        res[[i]][[5]] <- "NA"
      } else {
        res[[i]][[5]] <- probability
      }
      
    }
    names(res[[i]])[[5]] <- "probability"
    
    
    # support valule / threshold
    if (is.element("support", names(fuc))) {
      res[[i]][[6]] <- fuc$support
    } else {
      if (is.na(support) | is.null(support)) {
        res[[i]][[6]] <- "NA"
      } else {
        res[[i]][[6]] <- support
      }
      
    }
    names(res[[i]])[[6]] <- "support"
    
    # time step
    if (is.element("timestep", names(fuc))) {
      res[[i]][[7]] <- fuc$timestep
    } else {
      if (is.na(timestep) | is.null(timestep)) {
        res[[i]][[7]] <- "NA"
      } else {
        res[[i]][[7]] <- timestep
      }
      
    }
    names(res[[i]])[[7]] <- "timestep"
    
    if (length(functions) == 1) {
      names(res)[[i]] <- paste(c(name, "_", typename), collapse = "")
    } else {
      names(res)[[i]] <- paste(c(name, "_", i, "_", typename), collapse = "")
    }
  }
  
  return(res)
}


#' An internal function to convert the mined result to FBN network
#' 
#' @param minerresult The result of mining.
#' @param genes The genes involved in the mining.
#' 
convertMinedResultToFBNNetwork <- function(minerresult, genes) {
  futile.logger::flog.info(sprintf("Enter convertMinedResultToFBNNetwork zone"))
  # code
  tryCatch({
    res <- list()
    entry <- list()
    res[[1]] <- genes
    names(res)[[1]] <- "genes"
    
    res[[2]] <- list()
    names(res)[[2]] <- "interactions"
    res$interactions <- list()
    
    res[[3]] <- vapply(genes, function(gene) -1, numeric(1))
    names(res)[[3]] <- "fixed"
    
    res[[4]] <- vapply(genes, function(gene) 1, numeric(1))
    names(res)[[4]] <- "timedecay"
    
    # lapply(network$interactions,function(interaction){
    
    for (name in genes) {
      entry[[length(entry) + 1]] <- list()
      names(entry)[[length(entry)]] <- name
      
      # entry[name][[1]]<-list()
      verifyItems <- unlist(minerresult[name])
      if (!is.null(verifyItems)) {
        interactionItems <- dissolve(minerresult[name])
        
        iniIndex <- 0L
        for (j in seq_along(interactionItems)) {
          expression <- interactionItems[[j]][[2]]
          error <- interactionItems[[j]][[5]]
          probability <- interactionItems[[j]][[6]]
          support <- interactionItems[[j]][[7]]
          type <- as.numeric(interactionItems[[j]][[3]])
          timestep <- interactionItems[[j]][[8]]
          
          interactions <- regenerateInteractions(paste(c(name, "_", j), collapse = ""), expression, unlist(genes), error, type, probability, 
          support, timestep)
          for (i in seq_along(interactions)) {
          iniIndex <- iniIndex + 1
          entryitem <- interactions[[i]]
          entry[[name]][[iniIndex]] <- entryitem
          names(entry[[name]])[[iniIndex]] <- names(interactions[i])[[1]]
          }
        }
      }
    }
    res$interactions <- entry
    class(res) <- c("FundamentalBooleanNetwork", class(res))
  }, error = function(e) {
    stop(sprintf("Error converting to FBN Network \"%s\": %s", expression, e$message))
  })
  futile.logger::flog.info(sprintf("Leave convertMinedResultToFBNNetwork zone"))
  return(res)
}

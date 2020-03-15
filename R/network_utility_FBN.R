
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

removeFistLevelBracket <- function(splitedExpression) {
    res <- splitedExpression
    if (isSubExpression(splitedExpression)) {
        res <- res[-1]
        res <- res[-length(res)]
    }
    
    return(res)
}

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

convertIntoExpressionTree <- function(splitedExpression) {
    splitexp <- removeFistLevelBracket(splitedExpression)
    Oprators <- c("|", "&")
    
    groupA <- c("(", "[", "{")
    groupB <- c(")", "]", "}")
    
    negation <- "!"
    
    isTop <- FALSE
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
#' coming later
#' @export
generateFBNInteraction <- function(expressionString, genes) {
    res <- list()
    splitedexpression <- splitExpression(expressionString, 1, FALSE)
    geneinputs <- which(genes %in% splitedexpression)
    res$input <- geneinputs
    res$expression <- expressionString
    return(res)
}

# internal functions
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



# bug!
#'Convert a boolean network object to fundamental boolean function object
#'
#'@param network A boolean network object
#'@return An object of FBN network
#'@examples
#' coming later
#' @export
convertToFBNNetwork <- function(network) {
    # validate network types
    if (!(inherits(network, "BooleanNetworkCollection"))) 
        stop("Network must be inherited from BooleanNetwork")
    
    # code
    tryCatch({
        res <- list()
        entry <- list()
        res[[1]] <- network$genes
        names(res)[[1]] <- "genes"
        
        res[[2]] <- list()
        names(res)[[2]] <- "interactions"
        res$interactions <- list()
        
        res[[3]] <- network$fixed
        names(res)[[3]] <- "fixed"
        
        res[[4]] <- sapply(network$genes, function(gene) gene = 1)
        names(res)[[4]] <- "timedecay"
        
        # lapply(network$interactions,function(interaction){
        
        for (name in names(network$interactions)) {
            entry[[length(entry) + 1]] <- list()
            names(entry)[[length(entry)]] <- name
            # entry[name][[1]]<-list()
            interactionItems <- network$interactions[name]
            if (length(interactionItems) > 1) 
                stop("duplicate gene list of %s found", name)
            
            # entry[name][[length(entry[name])+1]]<-list()
            item <- interactionItems[[1]]
            iniIndex <- 0L
            for (j in seq_along(item)) {
                expression <- item[[j]]$expression
                error <- item[[j]]$error
                probability <- 1
                timestep <- 1
                support <- 1
                if (is.element("type", names(item[[j]]))) {
                  type <- item[[j]]$type
                } else {
                  type = NULL
                }
                
                if (is.null(error) | length(error) == 0) {
                  error <- 0
                } else {
                  probability <- 1 - as.numeric(error)
                }
                
                interactions <- regenerateInteractions(paste(c(name, "_", j), collapse = ""), expression, unlist(network$genes), 
                  error, type, probability, support, timestep)
                # if(length(item)==1L) { entry[[name]][[length(entry[[name]])+1]]<-unlist(interactions,recursive=FALSE) }else {
                for (i in seq_along(interactions)) {
                  iniIndex <- iniIndex + 1
                  entryitem <- interactions[[i]]
                  entry[[name]][[iniIndex]] <- entryitem
                  names(entry[[name]])[[iniIndex]] <- names(interactions[i])[[1]]
                }
                # }
                
            }
            
        }
        res$interactions <- entry
        class(res) <- c("FundamentalBooleanNetwork", class(res))
    }, error = function(e) {
        stop(sprintf("Error converting to FBN Network \"%s\": %s", expression, e$message))
    })
    
    return(res)
}

#'convert the mined result to FBNNetwork object
#'
#'@param minerresult The mined result from FBNmine.R
#'@param genes The target genes in the output
#'@return A object of FBN network
#'@examples
#' mat1<-matrix(c('1','0','0','1','0','0','0','1','1'),3,3, dimnames=list(c('gene1','gene2','gene3'),c('1','2','3')))
#' mat2<-matrix(c('1','1','0','1','0','1','1','1','0'),3,3, dimnames=list(c('gene1','gene2','gene3'),c('1','2','3')))
#' listtest<-list(mat1,mat2)
#' cube<-constructFBNCube(c('gene1','gene2'),c('gene1','gene2','gene3'),listtest,4,1,FALSE)
#' res<-mineFBNNetworkCore(cube,c('gene1','gene2'))
#' res<-mineFBNNetworkStage2(res)
#' finalresult<-outputMinerResultToFBNNetwork(res,c('gene1','gene2'))
#' @export
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
        
        res[[3]] <- sapply(genes, function(gene) -1)
        names(res)[[3]] <- "fixed"
        
        res[[4]] <- sapply(genes, function(gene) 1)
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
                  
                  interactions <- regenerateInteractions(paste(c(name, "_", j), collapse = ""), expression, unlist(genes), 
                    error, type, probability, support, timestep)
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
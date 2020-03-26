#'This method is used to convert the FBN network into a FBN networkobject for graphs
#'
#'@param FBNnetwork A FBN network object
#'@param david_gene_list A data frame of DAVID biological annotation
#'@param show_decay Optional, if TRUE, will draw the connections for decay.
#'@return An network graphic object
#'@examples
#' ##coming later
ConvertToNetworkGraphicObject <- function(FBNnetwork,
                                          david_gene_list = NULL, 
                                          show_decay = FALSE) {
  if (!(inherits(FBNnetwork, "FundamentalBooleanNetwork"))) 
    stop("Network must be inherited from FundamentalBooleanNetwork")
  
  genes <- FBNnetwork$genes
  
  nodes <- list()
  nodes[[1]] <- genes
  names(nodes)[[1]] <- "id"
  
  nodes[[2]] <- rep("ellipse", length(genes))
  names(nodes)[[2]] <- "shape"
  
  nodes[[3]] <- rep("lightblue", length(genes))
  names(nodes)[[3]] <- "color"
  
  nodes[[4]] <- rep("gene", length(genes))
  names(nodes)[[4]] <- "type"
  
  nodes[[5]] <- rep(4, length(genes))
  names(nodes)[[5]] <- "size"
  
  nodes[[6]] <- genes
  names(nodes)[[6]] <- "label"
  
  nodes[[7]] <- rep(TRUE, length(genes))
  names(nodes)[[7]] <- "shadow"
  
  nodes[[8]] <- rep("Gene", length(genes))
  names(nodes)[[8]] <- "group"
  
  if (is.null(david_gene_list)) {
    nodes[[9]] <- genes
  } else {
    gene_details <- sapply(genes, function(gene, david_gene_list) {
      res <- as.character(david_gene_list[david_gene_list$Symbol %in% c(gene), ]$Gene.Name[1])
      if (is.null(res)) {
        gene
      } else {
        res
      }
      
    }, david_gene_list)
    
    nodes[[9]] <- gene_details
  }
  names(nodes)[[9]] <- "title"
  
  # converted into igraph type network
  links_1 <- c(NULL)
  links_2 <- c(NULL)
  links_3 <- c(NULL)
  links_4 <- c(NULL)
  links_5 <- c(NULL)
  links_6 <- c(NULL)
  links_7 <- c(NULL)
  links_8 <- c(NULL)
  links_9 <- c(NULL)
  links_10 <- c(NULL)
  links_11 <- c(NULL)
  links_12 <- c(NULL)
  links_13 <- c(NULL)
  links_14 <- c(NULL)
  links_15 <- c(NULL)
  links_16 <- c(NULL)
  for (i in seq_along(FBNnetwork$interactions)) {
    targetNode <- names(FBNnetwork$interactions)[[i]]
    if (show_decay) {
      # crate decay function create link part 1
      links_1 <- c(links_1, targetNode)
      links_2 <- c(links_2, targetNode)
      links_3 <- c(links_3, 0)
      links_4 <- c(links_4, targetNode)
      links_5 <- c(links_5, targetNode)
      links_6 <- c(links_6, "none")
      links_7 <- c(links_7, "odot")
      links_8 <- c(links_8, "grey")
      links_9 <- c(links_9, "longdash")
      links_10 <- c(links_10, 2)
      links_11 <- c(links_11, 0.2)
      links_12 <- c(links_12, "decay")
      links_13 <- c(links_13, "to")
      links_14 <- c(links_14, "")
      links_15 <- c(links_15, TRUE)
      links_16 <- c(links_16, TRUE)
    }
    
    for (j in seq_along(FBNnetwork$interactions[[i]])) {
      sourceNode <- names(FBNnetwork$interactions[[i]])[[j]]
      nodes[[1]] <- c(nodes[[1]], sourceNode)  #(length(nodes[[1]])+1))
      names(nodes)[[1]] <- "id"
      
      nodes[[2]] <- c(nodes[[2]], "box")
      names(nodes)[[2]] <- "shape"
      
      timestep <- FBNnetwork$interactions[[i]][[j]]$timestep
      if (as.numeric(FBNnetwork$interactions[[i]][[j]]$type) == 0) {
        arrowhead1 <- "odot"
        linecolor1 <- "darkred"
        nodes[[3]] <- c(nodes[[3]], "orange")
        names(nodes)[[3]] <- "color"
        linetype1 <- "longdash"
        links_14 <- c(links_14, "")
        nodes[[8]] <- c(nodes[[8]], "Inhibit Function")
        names(nodes)[[8]] <- "group"
        links_15 <- c(links_15, TRUE)
        nodes[[6]] <- c(nodes[[6]], paste0("-, ", timestep, sep = ""))
        names(nodes)[[6]] <- "label"
      } else {
        arrowhead1 <- "open"
        linecolor1 <- "darkblue"
        nodes[[3]] <- c(nodes[[3]], "LightGreen")
        names(nodes)[[3]] <- "color"
        linetype1 <- "solid"
        links_14 <- c(links_14, "")
        nodes[[8]] <- c(nodes[[8]], "Activate Function")
        names(nodes)[[8]] <- "group"
        links_15 <- c(links_15, FALSE)
        nodes[[6]] <- c(nodes[[6]], paste0("+, ", timestep, sep = ""))
        names(nodes)[[6]] <- "label"
      }
      
      nodes[[4]] <- c(nodes[[4]], "TF")
      names(nodes)[[4]] <- "type"
      
      nodes[[5]] <- c(nodes[[5]], 2)
      names(nodes)[[5]] <- "value"
      
      
      nodes[[7]] <- c(nodes[[7]], FALSE)
      names(nodes)[[7]] <- "shadow"
      
      
      nodes[[9]] <- c(nodes[[9]], FBNnetwork$interactions[[i]][[j]]$expression)
      names(nodes)[[9]] <- "title"
      
      arrowtail1 <- "none"
      
      support1 <- as.numeric(FBNnetwork$interactions[[i]][[j]]$support)
      
      # create link part 1
      links_1 <- c(links_1, sourceNode)
      links_2 <- c(links_2, targetNode)
      links_3 <- c(links_3, support1)
      links_4 <- c(links_4, targetNode)
      links_5 <- c(links_5, sourceNode)
      links_6 <- c(links_6, arrowtail1)
      links_7 <- c(links_7, arrowhead1)
      links_8 <- c(links_8, linecolor1)
      links_9 <- c(links_9, linetype1)
      links_10 <- c(links_10, 2)
      links_11 <- c(links_11, 0.2)
      links_12 <- c(links_12, "TF_to_Gene")
      links_13 <- c(links_13, "to")
      links_16 <- c(links_16, TRUE)
      
      for (n in seq_along(FBNnetwork$interactions[[i]][[j]]$input)) {
        nodeto <- names(FBNnetwork$interactions[[i]])[[j]]
        inputIndex <- FBNnetwork$interactions[[i]][[j]]$input[[n]]
        gene <- FBNnetwork$genes[[inputIndex]]
        nodefrom <- gene
        
        splitedexpression <- splitExpression(FBNnetwork$interactions[[i]][[j]]$expression, 1, FALSE)
        inputgeneindex <- which(splitedexpression == gene)
        
        arrowhead2 <- "none"
        if (length(inputgeneindex) == 0) {
          stop("The input gene index should not be zero")
        }
        if (inputgeneindex >= 2) {
          if (splitedexpression[inputgeneindex - 1] == "!") {
          arrowtail2 <- "tee"
          linecolor2 <- "red"
          links_14 <- c(links_14, "")
          links_15 <- c(links_15, TRUE)
          } else {
          arrowtail2 <- "none"
          linecolor2 <- "green"
          links_14 <- c(links_14, "")
          links_15 <- c(links_15, FALSE)
          }
        } else {
          arrowtail2 <- "none"
          linecolor2 <- "green"
          links_14 <- c(links_14, "")
          links_15 <- c(links_15, FALSE)
        }
        
        if (as.numeric(FBNnetwork$interactions[[i]][[j]]$type) == 0) {
          linetype2 <- "longdash"
        } else {
          linetype2 <- "solid"
        }
        
        support2 <- as.numeric(FBNnetwork$interactions[[i]][[j]]$support)
        
        # create link part 2
        links_1 <- c(links_1, nodefrom)
        links_2 <- c(links_2, nodeto)
        links_3 <- c(links_3, support2)
        links_4 <- c(links_4, targetNode)
        links_5 <- c(links_5, sourceNode)
        links_6 <- c(links_6, arrowtail2)
        links_7 <- c(links_7, arrowhead2)
        links_8 <- c(links_8, linecolor2)
        links_9 <- c(links_9, linetype2)
        links_10 <- c(links_10, 0)
        links_11 <- c(links_11, 0.2)
        links_12 <- c(links_12, "Gene_to_TF")
        links_13 <- c(links_13, "none")
        links_16 <- c(links_16, FALSE)
      }
    }
  }
  
  
  links <- list()
  links[[1]] <- links_1
  links[[2]] <- links_2
  links[[3]] <- links_3
  links[[4]] <- links_4
  links[[5]] <- links_5
  links[[6]] <- links_6
  links[[7]] <- links_7
  links[[8]] <- links_8
  links[[9]] <- links_9
  links[[10]] <- links_10
  links[[11]] <- links_11
  links[[12]] <- links_12
  links[[13]] <- links_13
  links[[14]] <- links_14
  links[[15]] <- links_15
  links[[16]] <- links_16
  
  names(links)[[1]] <- "from"
  names(links)[[2]] <- "to"
  names(links)[[3]] <- "support"
  names(links)[[4]] <- "targetNode"
  names(links)[[5]] <- "title"
  names(links)[[6]] <- "arrowtail"
  names(links)[[7]] <- "arrowhead"
  names(links)[[8]] <- "color"
  names(links)[[9]] <- "lty"
  names(links)[[10]] <- "arrow.mode"
  names(links)[[11]] <- "arrow.size"
  names(links)[[12]] <- "type"
  names(links)[[13]] <- "arrows"
  names(links)[[14]] <- "label"
  names(links)[[15]] <- "dashes"
  names(links)[[16]] <- "shadow"
  res <- list()
  res[[1]] <- as.data.frame(nodes)
  names(res)[[1]] <- "nodes"
  res[[2]] <- as.data.frame(links)
  names(res)[[2]] <- "edges"
  res
}

#' internal method to convert the timeseries data 
#' into the graphic object.
#' 
#' @param timeseries The timeseries data
#' @param FBNnetwork The fundamental Boolean network.
#' @param networkobject The fundamental Boolean network
#' objects.
#' @return  A dataframe object.
ConvertToDynamicNetworkGraphicObject <- function(timeseries, FBNnetwork, networkobject) {
  
  if (is.null(timeseries) | !is.matrix(timeseries)) {
    stop("The parameter timeseries must be not NULL and is a class of Matrix")
  }
  
  if (!(inherits(FBNnetwork, "FundamentalBooleanNetwork"))) {
    stop("Network must be inherited from FundamentalBooleanNetwork")
  }
  
  
  if (is.null(networkobject)) {
    stop("The networkobject is required")
  }
  
  timepoints <- colnames(timeseries)
  genes <- rownames(timeseries)
  interactions <- FBNnetwork$interactions
  
  
  onset <- c()
  terminus <- c()
  tail <- c()
  head <- c()
  onsetcensored <- c()
  terminuscensored <- c()
  duration <- c()
  edgeid <- c()
  fromGeneState <- c()
  temporarymat <- list()
  functypes <- c()
  for (i in seq_along(timepoints)) {
    decayIndex <- rep(1L, length(genes))
    names(decayIndex) <- genes
    
    if (i > 1) {
      premat <- timeseries[, i - 1]
      
      for (j in seq_along(genes)) {
        gene <- genes[[j]]
        ini <- premat[[gene]] == 1L
        decay <- FBNnetwork$timedecay[[gene]]
        genefunctions <- interactions[[gene]]
        # find all activators' probabilities
        condOfActivation <- sapply(genefunctions, function(activator) activator$type == 1L)
        funcOfActivators <- genefunctions[condOfActivation]
        # find all inhibitors' probabilities
        condOfInhibitors <- sapply(genefunctions, function(inhibitor) inhibitor$type == 0L)
        funcOfInhibitors <- genefunctions[condOfInhibitors]
        
        prFA <- FALSE
        selectedActivationFunctions <- list()
        
        if (length(funcOfActivators) > 0) {
          for (fa in seq_along(funcOfActivators)) {
          fbnName <- names(funcOfActivators)[[fa]]
          
          pregeneInput <- dissolve(lapply(funcOfActivators[[fa]]$input, function(geneindex) {
            res <- list()
            res[[1]] <- premat[[geneindex]]
            names(res)[[1]] <- genes[[geneindex]]
            return(res)
          }))
          
          probability <- getProbabilityFromFunctionInput(1, funcOfActivators[[fa]]$expression, funcOfActivators[[fa]]$probability, pregeneInput)
          
          if (!length(probability) == 0) {
            result <- randomSelection(probability)
            prFA <- prFA | result
            if (result) {
            temp_res <- list()
            temp_res[[1]] <- fbnName
            temp_res[[2]] <- pregeneInput
            names(temp_res)[[1]] <- "fuc"
            names(temp_res)[[2]] <- "input"
            selectedActivationFunctions[[length(selectedActivationFunctions) + 1]] <- temp_res
            }
          }
          }
        }
        
        prFD <- FALSE
        selectedInhibitionFunctions <- list()
        if (length(funcOfInhibitors) > 0) {
          for (fd in seq_along(funcOfInhibitors)) {
          fbnName <- names(funcOfInhibitors)[[fd]]
          pregeneInput <- dissolve(lapply(funcOfInhibitors[[fd]]$input, function(geneindex) {
            res <- list()
            res[[1]] <- premat[[geneindex]]
            names(res)[[1]] <- genes[[geneindex]]
            return(res)
          }))
          
          probability <- getProbabilityFromFunctionInput(0, funcOfInhibitors[[fd]]$expression, funcOfInhibitors[[fd]]$probability, pregeneInput)
          
          if (!length(probability) == 0) {
            result <- randomSelection(probability)
            prFD <- prFD | result
            if (result) {
            temp_res <- list()
            temp_res[[1]] <- fbnName
            temp_res[[2]] <- pregeneInput
            names(temp_res)[[1]] <- "fuc"
            names(temp_res)[[2]] <- "input"
            selectedInhibitionFunctions[[length(selectedInhibitionFunctions) + 1]] <- temp_res
            }
          }
          }
        }
        
        # nothing happen then deactivate when the decay time is due
        decayFunctions <- list()
        if (decay > 0) {
          if (prFA | prFD) {
          decayIndex[[gene]] <- 1L
          } else {
          # check decay
          if (as.numeric(decayIndex[[gene]]) >= as.numeric(decay)) {
            
            decayIndex[[gene]] <- 1L
            input <- list(ini)
            names(input)[[1]] <- gene
            temp_res <- list()
            temp_res[[1]] <- "decay"
            temp_res[[2]] <- input
            names(temp_res)[[1]] <- "fuc"
            names(temp_res)[[2]] <- "input"
            decayFunctions[[length(decayFunctions) + 1]] <- temp_res
            # ini<-FALSE
          } else {
            decayIndex[[gene]] <- as.numeric(decayIndex[[gene]]) + 1L
          }
          }
        }
        
        if (length(selectedInhibitionFunctions) > 0) {
          selectedActivationFunctions <- list()
        }
        # activation
        if (length(selectedActivationFunctions) > 0) {
          for (a in seq_along(selectedActivationFunctions)) {
          selectedActivationFunction <- selectedActivationFunctions[[a]]$fuc
          edges <- networkobject$edges[which(networkobject$edges$title == selectedActivationFunction & networkobject$edges$targetNode == 
            gene), ]
          functioninput <- selectedActivationFunctions[[a]]$input
          for (e in seq_along(rownames(edges))) {
            rowvalue <- edges[e, ]
            rowid <- rownames(rowvalue)
            fromid <- rownames(networkobject$nodes[which(networkobject$nodes$id == rowvalue$from), ])[[1]]
            toid <- rownames(networkobject$nodes[which(networkobject$nodes$id == rowvalue$to), ])[[1]]
            
            fromState <- NA
            if (!is.null(functioninput[[rowvalue$from]])) {
            fromState <- as.numeric(functioninput[[rowvalue$from]])
            }
            if (is.null(temporarymat[[rowid]])) {
            newindx <- length(temporarymat) + 1
            temporarymat[[newindx]] <- list()
            names(temporarymat)[[newindx]] <- rowid
            
            temporarymat[[newindx]][[1]] <- c(i - 1, i, as.numeric(fromid), as.numeric(toid), fromState, selectedActivationFunction, 
              "activate")
            names(temporarymat[[newindx]])[[1]] <- i
            } else {
            newindex <- length(temporarymat[[rowid]]) + 1
            temporarymat[[rowid]][[newindex]] <- c(i - 1, i, as.numeric(fromid), as.numeric(toid), fromState, selectedActivationFunction, 
              "activate")
            names(temporarymat[[rowid]])[[newindex]] <- i
            }
            
          }
          }
        }
        
        # inhibition
        if (length(selectedInhibitionFunctions) > 0) {
          for (d in seq_along(selectedInhibitionFunctions)) {
          selectedInhibitionFunction <- selectedInhibitionFunctions[[d]]$fuc
          edges <- networkobject$edges[which(networkobject$edges$title == selectedInhibitionFunction & networkobject$edges$targetNode == 
            gene), ]
          functioninput <- selectedInhibitionFunctions[[d]]$input
          for (e in seq_along(rownames(edges))) {
            rowvalue <- edges[e, ]
            rowid <- rownames(rowvalue)
            
            fromid <- rownames(networkobject$nodes[which(networkobject$nodes$id == rowvalue$from), ])[[1]]
            toid <- rownames(networkobject$nodes[which(networkobject$nodes$id == rowvalue$to), ])[[1]]
            
            fromState <- NA
            if (!is.null(functioninput[[rowvalue$from]])) {
            fromState <- as.numeric(functioninput[[rowvalue$from]])
            }
            
            if (is.null(temporarymat[[rowid]])) {
            newindx <- length(temporarymat) + 1
            temporarymat[[newindx]] <- list()
            names(temporarymat)[[newindx]] <- rowid
            
            temporarymat[[newindx]][[1]] <- c(i - 1, i, as.numeric(fromid), as.numeric(toid), fromState, selectedInhibitionFunction, 
              "inhibit")
            names(temporarymat[[newindx]])[[1]] <- i
            } else {
            newindex <- length(temporarymat[[rowid]]) + 1
            temporarymat[[rowid]][[newindex]] <- c(i - 1, i, as.numeric(fromid), as.numeric(toid), fromState, selectedInhibitionFunction, 
              "inhibit")
            names(temporarymat[[rowid]])[[newindex]] <- i
            }
          }
          }
        }
        
        # decay
        
        if (length(decayFunctions) > 0) {
          for (d in seq_along(decayFunctions)) {
          selecteddecayFunction <- decayFunctions[[d]]$fuc
          edges <- networkobject$edges[which(networkobject$edges$type == selecteddecayFunction & networkobject$edges$targetNode == gene), 
            ]
          functioninput <- decayFunctions[[d]]$input
          for (e in seq_along(rownames(edges))) {
            rowvalue <- edges[e, ]
            rowid <- rownames(rowvalue)
            
            fromid <- rownames(networkobject$nodes[which(networkobject$nodes$id == rowvalue$from), ])[[1]]
            toid <- rownames(networkobject$nodes[which(networkobject$nodes$id == rowvalue$to), ])[[1]]
            
            fromState <- NA
            if (!is.null(functioninput[[rowvalue$from]])) {
            fromState <- as.numeric(functioninput[[rowvalue$from]])
            }
            
            if (is.null(temporarymat[[rowid]])) {
            newindx <- length(temporarymat) + 1
            temporarymat[[newindx]] <- list()
            names(temporarymat)[[newindx]] <- rowid
            
            temporarymat[[newindx]][[1]] <- c(i - 1, i, as.numeric(fromid), as.numeric(toid), fromState, selecteddecayFunction, "decay")
            names(temporarymat[[newindx]])[[1]] <- i
            } else {
            newindex <- length(temporarymat[[rowid]]) + 1
            temporarymat[[rowid]][[newindex]] <- c(i - 1, i, as.numeric(fromid), as.numeric(toid), fromState, selecteddecayFunction, 
              "decay")
            names(temporarymat[[rowid]])[[newindex]] <- i
            }
          }
          }
        }
      }
    }
  }
  
  
  edgeids <- names(temporarymat)
  for (i in seq_along(edgeids)) {
    tempid <- edgeids[[i]]
    timesteps <- temporarymat[[tempid]]
    firstrecord <- timesteps[[1]]
    
    if (length(timesteps) == 1) {
      onset <- c(onset, firstrecord[[1]])
      terminus <- c(terminus, firstrecord[[2]])
      tail <- c(tail, as.numeric(firstrecord[[3]]))
      head <- c(head, as.numeric(firstrecord[[4]]))
      fromGeneState <- c(fromGeneState, firstrecord[[5]])
      functypes <- c(functypes, firstrecord[[7]])
      onsetcensored <- c(onsetcensored, FALSE)
      terminuscensored <- c(terminuscensored, FALSE)
      duration <- c(duration, (as.numeric(firstrecord[[2]]) - as.numeric(firstrecord[[1]])))
      edgeid <- c(edgeid, tempid)
    } else {
      index <- 2
      while (index <= length(timesteps)) {
        nextrecord <- timesteps[[index]]
        united <- FALSE
        if (nextrecord[[1]] == firstrecord[[2]]) {
          firstrecord[[2]] <- nextrecord[[2]]
          united <- TRUE
        } else {
          onset <- c(onset, firstrecord[[1]])
          terminus <- c(terminus, firstrecord[[2]])
          tail <- c(tail, as.numeric(firstrecord[[3]]))
          head <- c(head, as.numeric(firstrecord[[4]]))
          fromGeneState <- c(fromGeneState, firstrecord[[5]])
          onsetcensored <- c(onsetcensored, FALSE)
          terminuscensored <- c(terminuscensored, FALSE)
          duration <- c(duration, (as.numeric(firstrecord[[2]]) - as.numeric(firstrecord[[1]])))
          functypes <- c(functypes, firstrecord[[7]])
          edgeid <- c(edgeid, tempid)
          firstrecord <- nextrecord
        }
        
        if (index == length(timesteps) && united) {
          onset <- c(onset, firstrecord[[1]])
          terminus <- c(terminus, firstrecord[[2]])
          tail <- c(tail, as.numeric(firstrecord[[3]]))
          head <- c(head, as.numeric(firstrecord[[4]]))
          fromGeneState <- c(fromGeneState, firstrecord[[5]])
          onsetcensored <- c(onsetcensored, FALSE)
          terminuscensored <- c(terminuscensored, FALSE)
          duration <- c(duration, (as.numeric(firstrecord[[2]]) - as.numeric(firstrecord[[1]])))
          functypes <- c(functypes, firstrecord[[7]])
          edgeid <- c(edgeid, tempid)
          firstrecord <- nextrecord
        }
        index <- index + 1
      }
    }
    
  }
  
  dynamicnet <- list()
  dynamicnet[[1]] <- onset
  names(dynamicnet)[[1]] <- "onset"
  dynamicnet[[2]] <- terminus
  names(dynamicnet)[[2]] <- "terminus"
  dynamicnet[[3]] <- tail
  names(dynamicnet)[[3]] <- "tail"
  dynamicnet[[4]] <- head
  names(dynamicnet)[[4]] <- "head"
  dynamicnet[[5]] <- onsetcensored
  names(dynamicnet)[[5]] <- "onset.censored"
  dynamicnet[[6]] <- terminuscensored
  names(dynamicnet)[[6]] <- "terminus.censored"
  dynamicnet[[7]] <- duration
  names(dynamicnet)[[7]] <- "duration"
  dynamicnet[[8]] <- edgeid
  names(dynamicnet)[[8]] <- "edge.id"
  dynamicnet[[9]] <- fromGeneState
  names(dynamicnet)[[9]] <- "fromState"
  dynamicnet[[10]] <- functypes
  names(dynamicnet)[[10]] <- "functype"
  return(as.data.frame(dynamicnet))
}


#'Display a dynamic network in a slice
#'
#' @param networkobject A net work object
#' @param timepoint The target time point
#' @param dynamicNetworkGraphicObject A dynamic network graphic object
#' @return No return
#' @examples
#' ##coming later
#' @export
StaticNetworkInSlice <- function(networkobject, timepoint, dynamicNetworkGraphicObject) {
  
  if (is.null(dynamicNetworkGraphicObject)) {
    stop("The dynamicNetworkGraphicObject is required")
  }
  
  
  if (is.null(networkobject)) {
    stop("The networkobject is required")
  }
  
  
  if (timepoint > 0) {
    filterddynamic <- dynamicNetworkGraphicObject[which(as.numeric(dynamicNetworkGraphicObject$onset) < timepoint & timepoint <= as.numeric(dynamicNetworkGraphicObject$terminus)), 
      ]
  } else {
    filterddynamic <- dynamicNetworkGraphicObject
  }
  
  
  newedges <- networkobject$edges[which(rownames(networkobject$edges) %in% filterddynamic$edge.id), ]
  
  filterednode <- networkobject$nodes[which(networkobject$nodes$type == "gene"), ]
  filterednode2 <- networkobject$nodes[which(networkobject$nodes$id %in% newedges$from | networkobject$nodes$id %in% newedges$to), ]
  
  newnodes <- unique(rbind.data.frame(filterednode, filterednode2))  # nodes data.frame for legend
  lnodes <- data.frame(label = c("Gene", "Activate Function (+)", "Inhibit Function (-)"), shape = c("ellipse", "box", "box"), color = c("lightblue", 
    "lightgreen", "orange"), shadow = c(TRUE, FALSE, FALSE))
  
  # edges data.frame for legend
  ledges <- data.frame(color = c("darkblue", "darkred", "green", "red"), label = c("activate", "inhibit", "activated input", "deactivated input"), 
    arrows = c("to", "to", "none", "none"), dashes = c(FALSE, TRUE, FALSE, TRUE), shadow = c(TRUE, TRUE, FALSE, FALSE))
  
  graph <- visNetwork::visNetwork(newnodes, newedges, main = paste("Fundamental Boolean Networks in the time step of ", " ", timepoint, sep = "", 
    collapse = ""))
  graph <- visNetwork::visLegend(graph, addEdges = ledges, addNodes = lnodes, useGroups = FALSE)
  visNetwork::visInteraction(graph, dragNodes = TRUE, dragView = TRUE, zoomView = TRUE)
}

#'Display a static network in a slice
#'
#'@param networkobject A FBN network object
#'@param david_gene_list the list object of DAVID annotation gene list
#'@return No return
#'@examples
#' ##coming later
#' @export
StaticNetwork <- function(networkobject, david_gene_list = NULL) {
  
  if (is.null(networkobject)) {
    stop("The networkobject is required")
  }
  
  newedges <- networkobject$edges
  newnodes <- networkobject$nodes
  
  lnodes <- data.frame(label = c("Gene", "Activate Function (+, Timestep)", "Inhibit Function (-, Timestep)"), shape = c("ellipse", "box", "box"), 
    color = c("lightblue", "lightgreen", "orange"), shadow = c(TRUE, FALSE, FALSE))
  
  # edges data.frame for legend
  ledges <- data.frame(color = c("darkblue", "darkred", "green", "red"), label = c("activate", "inhibit", "activated input", "deactivated input"), 
    arrows = c("to", "to", "none", "none"), dashes = c(FALSE, TRUE, FALSE, TRUE), shadow = c(TRUE, TRUE, FALSE, FALSE))
  
  
  graph <- visNetwork::visNetwork(newnodes, newedges, main = paste("Fundamental Boolean Networks", sep = "", collapse = ""))
  graph <- visNetwork::visLegend(graph, addEdges = ledges, addNodes = lnodes, useGroups = FALSE)
  visNetwork::visInteraction(graph, dragNodes = TRUE, dragView = TRUE, zoomView = TRUE)
}

#' An Internal method to generate dynamic network graphic objects
#' 
#' @param networkobject Fundamental Boolean Network objects.
#' @param fromtimepoint A start point of the timeseries matrix.
#' @param totimepoint The end point of the timeseries matrix.
#' @param dynamicNetworkGraphicObject dynamicNetworkGraphicObject.
#' @param timeseries The timeseries matrix
GenerateDynamicNetworkGraphicObject <- function(networkobject, fromtimepoint, totimepoint, dynamicNetworkGraphicObject, timeseries) {
  
  if (is.null(timeseries) | !is.matrix(timeseries)) {
    stop("The parameter timeseries must be not NULL and is a class of Matrix")
  }
  
  if (is.null(dynamicNetworkGraphicObject)) {
    stop("The dynamicNetworkGraphicObject is required")
  }
  
  if (is.null(networkobject)) {
    stop("The networkobject is required")
  }
  
  
  newedges <- NULL
  newnodes <- NULL
  index <- fromtimepoint
  while (index < totimepoint) {
    curstate <- timeseries[, index]
    nextindex <- index + 1
    nextstate <- timeseries[, nextindex]
    
    # fromfilterddynamic<-dynamicNetworkGraphicObject[which(as.numeric(dynamicNetworkGraphicObject$onset)<index &
    # index<=as.numeric(dynamicNetworkGraphicObject$terminus)),]
    
    tofilterddynamic <- dynamicNetworkGraphicObject[which(as.numeric(dynamicNetworkGraphicObject$onset) < nextindex & nextindex <= as.numeric(dynamicNetworkGraphicObject$terminus)), 
      ]
    
    
    tonewedges <- networkobject$edges[which(rownames(networkobject$edges) %in% tofilterddynamic$edge.id), ]
    
    # from
    filterednodeFrom <- networkobject$nodes[which(networkobject$nodes$id %in% tonewedges$from), ]
    filterednodeFromGene <- filterednodeFrom[which(filterednodeFrom$type == "gene"), ]
    filterednodeFromGene[which(curstate[filterednodeFromGene$id] == 0), ]$color = "pink"
    filterednodeFromGene[which(curstate[filterednodeFromGene$id] == 1), ]$color = "lightblue"
    
    tonewedges[which(tonewedges$from %in% filterednodeFromGene$id), ]$from = paste(tonewedges[which(tonewedges$from %in% filterednodeFromGene$id), 
      ]$from, "_", index, sep = "")
    filterednodeFromGene$id <- paste(filterednodeFromGene$id, "_", index, sep = "")
    filterednodeFromGene$label <- paste(filterednodeFromGene$label, "_", index, sep = "")
    
    # to
    filterednodeTo <- networkobject$nodes[which(networkobject$nodes$id %in% tonewedges$to), ]
    filterednodeToGene <- filterednodeTo[which(filterednodeTo$type == "gene"), ]
    filterednodeToGene[which(nextstate[filterednodeToGene$id] == 0), ]$color = "pink"
    filterednodeToGene[which(nextstate[filterednodeToGene$id] == 1), ]$color = "lightblue"
    
    tonewedges[which(tonewedges$to %in% filterednodeToGene$id), ]$to = paste(tonewedges[which(tonewedges$to %in% filterednodeToGene$id), ]$to, 
      "_", nextindex, sep = "")
    filterednodeToGene$id <- paste(filterednodeToGene$id, "_", nextindex, sep = "")
    filterednodeToGene$label <- paste(filterednodeToGene$label, "_", nextindex, sep = "")
    
    # TF
    filterednodeFromTF <- filterednodeFrom[which(filterednodeFrom$type == "TF"), ]
    tonewedges[which(tonewedges$from %in% filterednodeFromTF$id), ]$from = paste(tonewedges[which(tonewedges$from %in% filterednodeFromTF$id), 
      ]$from, "_", nextindex, sep = "")
    filterednodeFromTF$id <- paste(filterednodeFromTF$id, "_", nextindex, sep = "")
    
    filterednodeToTF <- filterednodeTo[which(filterednodeTo$type == "TF"), ]
    tonewedges[which(tonewedges$to %in% filterednodeToTF$id), ]$to = paste(tonewedges[which(tonewedges$to %in% filterednodeToTF$id), ]$to, "_", 
      nextindex, sep = "")
    filterednodeToTF$id <- paste(filterednodeToTF$id, "_", nextindex, sep = "")
    
    newnodes <- unique(rbind.data.frame(newnodes, filterednodeToGene, filterednodeFromGene, filterednodeFromTF, filterednodeToTF))  # nodes data.frame for legend
    newedges <- unique(rbind.data.frame(newedges, tonewedges))
    index <- index + 1
  }
  
  
  lnodes <- data.frame(label = c("Activated Gene", "Inhibited Gene", "Activate Function (+)", "Inhibit Function (-)"), shape = c("ellipse", "ellipse", 
    "box", "box"), color = c("lightblue", "pink", "lightgreen", "orange"), shadow = c(TRUE, TRUE, FALSE, FALSE))
  
  # edges data.frame for legend
  ledges <- data.frame(color = c("darkblue", "darkred", "grey", "green", "red"), label = c("activate", "inhibit", "decay", "activated input", 
    "deactivated input"), arrows = c("to", "to", "to", "none", "none"), dashes = c(FALSE, TRUE, TRUE, FALSE, TRUE), shadow = c(TRUE, TRUE, TRUE, 
    FALSE, FALSE))
  
  
  
  graph <- visNetwork::visNetwork(newnodes, newedges, main = paste("Dynamic Fundamental Boolean Networks from the time point of", " ", fromtimepoint, 
    " ", "to ", totimepoint, sep = "", collapse = ""))
  graph <- visNetwork::visLegend(graph, addEdges = ledges, addNodes = lnodes, useGroups = FALSE)
  
  visNetwork::visInteraction(graph, dragNodes = TRUE, dragView = TRUE, zoomView = TRUE)
}


#' An Internal method to draw atractor internally
#' 
#' @param networkobject Fundamental Boolean Network objects.
#' @param dynamicNetworkGraphicObject dynamicNetworkGraphicObject.
#' @param matrix The timeseries matrix
DrawingAttractorInternal <- function(networkobject, dynamicNetworkGraphicObject, matrix) {
  if (is.null(matrix) | !is.matrix(matrix)) {
    stop("The parameter matrix must be not NULL and is a class of Matrix")
  }
  
  if (is.null(dynamicNetworkGraphicObject)) {
    stop("The dynamicNetworkGraphicObject is required")
  }
  
  if (is.null(networkobject)) {
    stop("The networkobject is required")
  }
  
  genes <- rownames(matrix)
  cNames <- colnames(matrix)
  newedges <- NULL
  newnodes <- NULL
  index <- 1
  startlevel <- 1
  while (index < length(cNames)) {
    curstate <- matrix[, index]
    nextindex <- index + 1
    nextstate <- matrix[, nextindex]
    
    # fromfilterddynamic<-dynamicNetworkGraphicObject[which(as.numeric(dynamicNetworkGraphicObject$onset)<index &
    # index<=as.numeric(dynamicNetworkGraphicObject$terminus)),]
    
    tofilterddynamic <- dynamicNetworkGraphicObject[which(as.numeric(dynamicNetworkGraphicObject$onset) < nextindex & nextindex <= as.numeric(dynamicNetworkGraphicObject$terminus)), 
      ]
    
    
    tonewedges <- networkobject$edges[which(rownames(networkobject$edges) %in% tofilterddynamic$edge.id), ]
    
    # from
    filterednodeFrom <- networkobject$nodes[which(networkobject$nodes$id %in% tonewedges$from), ]
    filterednodeFromGene <- filterednodeFrom[which(filterednodeFrom$type == "gene"), ]
    filterednodeFromGene[which(curstate[filterednodeFromGene$id] == 0), ]$color = "pink"
    filterednodeFromGene[which(curstate[filterednodeFromGene$id] == 1), ]$color = "lightblue"
    
    tonewedges[which(tonewedges$from %in% filterednodeFromGene$id), ]$from = paste(tonewedges[which(tonewedges$from %in% filterednodeFromGene$id), 
      ]$from, "_", index, sep = "")
    filterednodeFromGene$id <- paste(filterednodeFromGene$id, "_", index, sep = "")
    filterednodeFromGene$label <- paste(filterednodeFromGene$label, "_", index, sep = "")
    filterednodeFromGene <- cbind(filterednodeFromGene, level = rep(startlevel, length(rownames(filterednodeFromGene))))
    
    involvedgenes <- unique(filterednodeFromGene$title)
    # no connection node
    independeGenes <- genes[which(!genes %in% involvedgenes)]
    # individualnode<-networkobject$nodes[which(!filterednodeFromGene$nodes$title%in%genes),] individualnode$color='grey'
    individualfromnode <- list()
    if (length(independeGenes) > 0) {
      for (k in seq_along(independeGenes)) {
        gene <- independeGenes[[k]]
        genestate <- matrix[independeGenes[[k]], index]
        if (as.numeric(genestate == 0)) {
          makeupnode <- list(id = paste(gene, "_", index, sep = ""), shape = "ellipse", color = "pink", type = "gene", value = 4, label = paste(gene, 
          "_", index, sep = ""), shadow = TRUE, group = "Gene", title = gene, level = startlevel)
        } else {
          makeupnode <- list(id = paste(gene, "_", index, sep = ""), shape = "ellipse", color = "lightblue", type = "gene", value = 4, label = paste(gene, 
          "_", index, sep = ""), shadow = TRUE, group = "Gene", title = gene, level = startlevel)
        }
        
        individualfromnode[[length(individualfromnode) + 1]] <- makeupnode
      }
      individualfromnode <- do.call(rbind.data.frame, individualfromnode)
      colnames(individualfromnode) <- c("id", "shape", "color", "type", "value", "label", "shadow", "group", "title", "level")
    }
    
    # TF from
    filterednodeFromTF <- filterednodeFrom[which(filterednodeFrom$type == "TF"), ]
    tonewedges[which(tonewedges$from %in% filterednodeFromTF$id), ]$from = paste(tonewedges[which(tonewedges$from %in% filterednodeFromTF$id), 
      ]$from, "_", nextindex, sep = "")
    filterednodeFromTF$id <- paste(filterednodeFromTF$id, "_", nextindex, sep = "")
    filterednodeFromTF <- cbind(filterednodeFromTF, level = rep(startlevel + 1, length(rownames(filterednodeFromTF))))
    
    
    
    # to
    filterednodeTo <- networkobject$nodes[which(networkobject$nodes$id %in% tonewedges$to), ]
    filterednodeToGene <- filterednodeTo[which(filterednodeTo$type == "gene"), ]
    filterednodeToGene[which(nextstate[filterednodeToGene$id] == 0), ]$color = "pink"
    filterednodeToGene[which(nextstate[filterednodeToGene$id] == 1), ]$color = "lightblue"
    
    tonewedges[which(tonewedges$to %in% filterednodeToGene$id), ]$to = paste(tonewedges[which(tonewedges$to %in% filterednodeToGene$id), ]$to, 
      "_", nextindex, sep = "")
    filterednodeToGene$id <- paste(filterednodeToGene$id, "_", nextindex, sep = "")
    filterednodeToGene$label <- paste(filterednodeToGene$label, "_", nextindex, sep = "")
    filterednodeToGene <- cbind(filterednodeToGene, level = rep(startlevel + 2, length(rownames(filterednodeToGene))))
    
    involvedgenes <- unique(filterednodeToGene$title)
    # no connection node
    independeGenes <- genes[which(!genes %in% involvedgenes)]
    # individualnode<-networkobject$nodes[which(!filterednodeFromGene$nodes$title%in%genes),] individualnode$color='grey'
    individualtonode <- list()
    if (length(independeGenes) > 0) {
      for (k in seq_along(independeGenes)) {
        gene <- independeGenes[[k]]
        genestate <- matrix[independeGenes[[k]], nextindex]
        if (as.numeric(genestate == 0)) {
          makeupnode <- list(id = paste(gene, "_", nextindex, sep = ""), shape = "ellipse", color = "pink", type = "gene", value = 4, label = paste(gene, 
          "_", nextindex, sep = ""), shadow = TRUE, group = "Gene", title = gene, level = startlevel + 2)
        } else {
          makeupnode <- list(id = paste(gene, "_", nextindex, sep = ""), shape = "ellipse", color = "lightblue", type = "gene", value = 4, 
          label = paste(gene, "_", nextindex, sep = ""), shadow = TRUE, group = "Gene", title = gene, level = startlevel + 2)
        }
        
        individualtonode[[length(individualtonode) + 1]] <- makeupnode
      }
      individualtonode <- do.call(rbind.data.frame, individualtonode)
      colnames(individualtonode) <- c("id", "shape", "color", "type", "value", "label", "shadow", "group", "title", "level")
    }
    
    # to TF
    
    filterednodeToTF <- filterednodeTo[which(filterednodeTo$type == "TF"), ]
    tonewedges[which(tonewedges$to %in% filterednodeToTF$id), ]$to = paste(tonewedges[which(tonewedges$to %in% filterednodeToTF$id), ]$to, "_", 
      nextindex, sep = "")
    filterednodeToTF$id <- paste(filterednodeToTF$id, "_", nextindex, sep = "")
    filterednodeToTF <- cbind(filterednodeToTF, level = rep(startlevel + 1, length(rownames(filterednodeToTF))))
    newnodes <- unique(rbind.data.frame(newnodes, individualfromnode, individualtonode, filterednodeToGene, filterednodeFromGene, filterednodeFromTF, 
      filterednodeToTF))  # nodes data.frame for legend
    # newnodes<-newnodes[order(newnodes$title),]
    newedges <- unique(rbind.data.frame(newedges, tonewedges))
    index <- index + 1
    startlevel <- startlevel + 2
  }
  
  
  lnodes <- data.frame(label = c("Activated Gene", "Inhibited Gene", "Activate Function (+)", "Inhibit Function (-)"), shape = c("ellipse", "ellipse", 
    "box", "box"), color = c("lightblue", "pink", "lightgreen", "orange"), shadow = c(TRUE, TRUE, FALSE, FALSE))
  
  # edges data.frame for legend
  ledges <- data.frame(color = c("darkblue", "darkred", "grey", "green", "red"), label = c("activate", "inhibit", "decay", "activated input", 
    "deactivated input"), arrows = c("to", "to", "to", "none", "none"), dashes = c(FALSE, TRUE, TRUE, FALSE, TRUE), shadow = c(TRUE, TRUE, TRUE, 
    FALSE, FALSE))
  
  
  
  graph <- visNetwork::visNetwork(newnodes, newedges, main = paste0("Dynamic Fundamental Boolean Networks from the time point of", " ", 1, " ", 
    "to ", length(cNames)))
  graph <- visNetwork::visLegend(graph, addEdges = ledges, addNodes = lnodes, useGroups = FALSE, width = 0.1)
  graph <- visNetwork::visInteraction(graph, dragNodes = TRUE, dragView = TRUE, zoomView = TRUE)
  graph <- visNetwork::visEdges(graph, arrows = "from")
  visNetwork::visHierarchicalLayout(graph, direction = "LR", levelSeparation = 200)
}



# type=static | dynamic | animation | staticSlice
#'Executing a list of processes in parallel with decreasing list
#'
#'@param fbnNetwork A FBN Network object
#'@param type A type of object and its value is type=static | dynamic | staticSlice
#'@param timeseriesMatrix A sample timeseries matrix
#'@param fromTimePoint The time point at the beginning
#'@param toTimePoint The time point at the end
#'@param networkobject An object created by the funciton
#' \code{ConvertToNetworkGraphicObject}.
#'@return No return
#'@examples
#' ##coming later
#' @export
FBNNetwork.Graph <- function(fbnNetwork, type = "static", timeseriesMatrix = NULL, fromTimePoint = 1, toTimePoint = 5, networkobject = NULL) {
  
  if (length(fbnNetwork$interactions) == 0) 
    return(NULL)
  
  utils::data("DAVID_gene_list", overwrite = TRUE)
  
  if (is.null(networkobject)) {
    if (type == "dynamic") 
      networkobject <- ConvertToNetworkGraphicObject(fbnNetwork, DAVID_gene_list, show_decay = TRUE) 
    else 
      networkobject <- ConvertToNetworkGraphicObject(fbnNetwork, DAVID_gene_list, show_decay = FALSE)
  }
  
  
  dynamicnetworkobject <- NULL
  if (!is.null(timeseriesMatrix)) {
    if (!is.matrix(timeseriesMatrix)) {
      stop("The parameter timeseriesMatrix should be a class of matrix!")
    }
    
    dynamicnetworkobject <- ConvertToDynamicNetworkGraphicObject(timeseriesMatrix, fbnNetwork, networkobject)
  }
  
  switch(type, 
         static = StaticNetwork(networkobject, DAVID_gene_list), 
         staticSlice = StaticNetworkInSlice(networkobject, toTimePoint, dynamicnetworkobject), 
         dynamic = GenerateDynamicNetworkGraphicObject(networkobject, fromTimePoint, toTimePoint, 
    dynamicnetworkobject, timeseriesMatrix))
}

#' A method to draw attractors in the form of dynamic Network
#' 
#' @param fbnNetwork An object of FBNNetwork
#' @param FBMAttractors A list of attractors extracted via 
#' \code{searchForAttractors}.
#' @param index The indiex of an attractor in the list 
#' \code{FBMAttractors}.
#'@examples
#' require(BoolNet)
#' data('ExampleNetwork')
#' initialStates<-generateAllCombinationBinary(ExampleNetwork$genes)
#' trainingseries<-genereateBoolNetTimeseries(ExampleNetwork,
#'                                            initialStates,
#'                                            43,
#'                                            type='synchronous')
#' cube<-constructFBNCube(target_genes = ExampleNetwork$genes,
#'                        conditional_genes = ExampleNetwork$genes,
#'                        timeseriesCube = trainingseries,
#'                        maxK = 4,
#'                        temporal = 1,
#'                        useParallel = FALSE)
#' NETWORK2<-mineFBNNetwork(cube,ExampleNetwork$genes)
#' attractor<-searchForAttractors(NETWORK2,
#'                                initialStates,
#'                                ExampleNetwork$genes)
#' print(attractor)
#' FBNNetwork.Graph.DrawAttractor(NETWORK2,attractor,2)
#' @export
FBNNetwork.Graph.DrawAttractor <- function(fbnNetwork, FBMAttractors, index = 1) {
  attractor <- FBMAttractors$Attractors[[index]]
  genes <- FBMAttractors$Genes
  networkobject <- ConvertToNetworkGraphicObject(fbnNetwork, show_decay = FALSE)
  if (!is.null(attractor)) {
    matrix <- do.call(cbind, attractor)
    rownames(matrix) <- genes
    len <- length(attractor)
    colnames(matrix) <- c(1:len)
    
    
    dynamicnetworkobject <- ConvertToDynamicNetworkGraphicObject(matrix,
                                                                 fbnNetwork,
                                                                 networkobject)
  }
  DrawingAttractorInternal(networkobject,
                           dynamicnetworkobject,
                           matrix)
}

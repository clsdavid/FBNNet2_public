#'Load FBN network data into a FBN network object
#'
#'@param file A file that contains FBN data
#'@param bodySeparator A specific text separator
#'@param lowercaseGenes if all genes are in lower case
#'@return An object of FBN network
#'@examples
#' coming later
#' @export
loadFBNNetwork <- function(file, bodySeparator = ",", lowercaseGenes = FALSE) {
    # internal function
    matchNames <- function(rule) {
        regexpr <- "([_a-zA-Z][_a-zA-Z0-9]*)[,| |\\)|\\||\\&|\\[]"
        rule <- paste(gsub(" ", "", rule, fixed = TRUE), " ", sep = "")
        res <- unique(unname(sapply(regmatches(rule, gregexpr(regexpr, rule))[[1]], function(m) {
            sapply(regmatches(m, regexec(regexpr, m)), function(x) x[2])
        })))
        
        # remove operators
        isOp <- sapply(res, function(x) {
            tolower(x) %in% c("all", "any", "sumis", "sumgt", "sumlt", "maj", "timegt", "timelt", "timeis")
        })
        
        return(res[!isOp])
    }
    
    func <- readLines(file, -1)
    func <- gsub("#.*", "", trimws(func))
    func <- func[nchar(func) > 0]
    if (length(func) == 0) 
        stop("Header expected!")
    
    header <- func[1]
    header <- tolower(trimws(strsplit(header, bodySeparator)[[1]]))
    
    if (length(header) < 3 || header[1] != "targets" || !(header[2] %in% c("functions", "factors")) || header[3] != "type") 
        stop(paste("Invalid header:", func[1]))
    
    func <- func[-1]
    if (lowercaseGenes) 
        func <- tolower(func)
    func <- gsub("[^\\[\\]a-zA-Z0-9_\\|\\&!\\(\\) \t\\-+=.,]+", "_", func, perl = TRUE)
    tmp <- unname(lapply(func, function(x) {
        bracketCount <- 0
        lastIdx <- 1
        chars <- strsplit(x, split = "")[[1]]
        res <- c()
        if (length(chars) > 0) {
            for (i in seq_along(chars)) {
                if (chars[i] == "(") 
                  bracketCount <- bracketCount + 1 else if (chars[i] == ")") 
                  bracketCount <- bracketCount - 1 else if (chars[i] == bodySeparator && bracketCount == 0) {
                  res <- c(res, trimws(paste(chars[lastIdx:(i - 1)], collapse = "")))
                  lastIdx <- i + 1
                }
            }
            res <- c(res, trimws(paste(chars[lastIdx:length(chars)], collapse = "")))
        }
        return(res)
    }))
    
    targets <- sapply(tmp, function(rule) trimws(rule[1]))
    for (target in targets) {
        if (regexec("^[a-zA-Z_][a-zA-Z0-9_]*$", target)[[1]] == -1) 
            stop(paste("Invalid gene name:", target))
    }
    
    factors <- sapply(tmp, function(rule) trimws(rule[2]))
    
    types <- sapply(tmp, function(rule) as.numeric(rule[3]))
    
    factors.tmp <- lapply(factors, matchNames)
    
    genes <- unique(c(targets, unname(unlist(factors.tmp))))
    
    
    fixed <- rep(-1, length(genes))
    names(fixed) <- genes
    interactions <- list()
    
    
    for (i in seq_along(targets)) {
        target <- targets[i]
        interaction <- generateFBNInteraction(factors[i], genes)
        if (length(interaction$func) == 1) {
            fixed[target] <- interaction$func
        }
        interaction[[length(interaction) + 1]] <- types[i]
        names(interaction)[[length(interaction)]] <- "type"
        interactions[[target]][[length(interactions[[target]]) + 1]] <- interaction
    }
    onlyInputs <- setdiff(genes, targets)
    if (length(onlyInputs) > 0) {
        for (gene in onlyInputs) {
            warning(paste("There is no transition function for gene \"", gene, "\"! Assuming an input!", sep = ""))
            
            interactions[[gene]] = list(input = length(interactions) + 1, func = c(0, 1), expression = gene)
        }
    }
    
    res <- list(interactions = interactions, genes = genes, fixed = fixed)
    
    class(res) <- c("BooleanNetworkCollection")
    
    # return(convertToFBNNetwork(res))
    return(res)
}

#' @export
convertToBooleanNetworkCollection <- function(network) {
    # validate network types
    if (!inherits(network, "BooleanNetwork")) 
        stop("Network1 must be inherited from BooleanNetwork")
    genes1 <- network$genes
    genes2 <- network$genes
    totalGenes <- network$genes
    merges <- mergeInteraction(lapply(network$interactions, convertInteraction), lapply(network$interactions, convertInteraction), genes1, genes2, totalGenes)
    res <- list(interactions = merges, genes = unique(c(network$genes, network$genes)), fixed = network$fixed, timedecay = rep(1, length(network$genes)))
    
    class(res) <- "BooleanNetworkCollection"
    return(res)
}

# need to fix the problem of genes and fixed's order when merge or filter
#' @export
mergeNetwork <- function(network1, network2) {
    # validate network types
    if (!(inherits(network1, "BooleanNetworkCollection")) & !(inherits(network1, "FundamentalBooleanNetwork"))) 
        stop("Network1 must be inherited from FundamentalBooleanNetwork or BooleanNetworkCollection")
    
    if (!(inherits(network2, "BooleanNetworkCollection")) & !(inherits(network2, "FundamentalBooleanNetwork"))) 
        stop("Network2 must be inherited from FundamentalBooleanNetwork or BooleanNetworkCollection")
    genes1 <- network1$genes
    genes2 <- network2$genes
    totalGenes <- unique(c(network1$genes, network2$genes))
    merges <- mergeInteraction(network1$interactions, network2$interactions, genes1, genes2, totalGenes)
    
    if (length(network1$fixed) < length(genes1)) {
        network1$fixed <- rep(-1, length(genes1))
    }
    
    fixed1 <- network1$fixed
    names(fixed1) <- genes1
    
    if (length(network2$fixed) < length(genes2)) {
        network2$fixed <- rep(-1, length(genes2))
    }
    fixed2 <- network2$fixed
    names(fixed2) <- genes2
    
    newfixed <- rep(-1, length(totalGenes))
    names(newfixed) <- totalGenes
    newfixed[which(totalGenes %in% names(fixed1))] <- fixed1
    newfixed[which(totalGenes %in% names(fixed2))] <- fixed2
    
    if (length(network1$timedecay) < length(genes1)) {
        network1$timedecay <- rep(-1, length(genes1))
    }
    
    timedecay1 <- network1$timedecay
    names(timedecay1) <- genes1
    
    if (length(network2$timedecay) < length(genes2)) {
        network2$timedecay <- rep(-1, length(genes2))
    }
    
    timedecay2 <- network2$timedecay
    names(timedecay2) <- genes2
    
    newtimedecay <- rep(-1, length(totalGenes))
    names(newtimedecay) <- totalGenes
    newtimedecay[which(totalGenes %in% names(timedecay1))] <- timedecay1
    newtimedecay[which(totalGenes %in% names(timedecay2))] <- timedecay2
    
    res <- list(interactions = merges, genes = totalGenes, fixed = newfixed, timedecay = newtimedecay)
    
    class(res) <- "FundamentalBooleanNetwork"
    return(res)
}

#' @export
convertInteraction <- function(interaction) {
    res <- list()
    if (mode(interaction) == "list") {
        res <- interaction
    } else if (mode(interaction) == "pairlist") {
        res <- list(list(input = interaction$input, expression = interaction$expression, error == interaction$error, type = NA, probability = (1 - as.numeric(interaction$error)), 
            support = NA, timestep = 1))
    }
    return(res)
}

# need to revise this as it produce duliplcates
mergeInteraction <- function(interactions1, interactions2, genes1, genes2, mergedgenes = NULL) {
    if (is.null(mergedgenes)) 
        mergedgenes <- sort(unique(c(genes1, genes2)))
    
    res <- list()
    ## loog for interactions in the first group

    for (name1 in names(interactions1)) {
        unique_express_1 <- c()
        unique_express_0 <- c()
        a_index = 1
        i_index = 1
        
        if (!(name1 %in% names(res))) {
            index <- length(res) + 1
            res[[index]] <- list()
            for (j in seq(interactions1[[name1]])) {
                inputgene1 <- genes1[interactions1[[name1]][[j]]$input]
                if (length(inputgene1) == 0) {
                    next
                }
                
                type <- interactions1[[name1]][[j]]$type
                expression <- interactions1[[name1]][[j]]$expression

                temp_exp <- paste(sort(splitExpression(expression, 2, TRUE)), sep = "", collapse = "")
                if (type == 1) {
                    if (temp_exp %in% unique_express_1) {
                        next
                    }
  
                    unique_express_1 <- c(unique_express_1, temp_exp)
                } else {
                    if (temp_exp %in% unique_express_0) {
                        next
                    }

                    unique_express_0 <- c(unique_express_0, temp_exp)
                }
                
                subindex <- length(res[[index]]) + 1
                newinput <- which(mergedgenes %in% inputgene1)
                if (length(newinput) == 0) {
                  next
                }
                res[[index]][[subindex]] <- list(input = newinput, 
                                                 expression = interactions1[[name1]][[j]]$expression, 
                                                 error = interactions1[[name1]][[j]]$error, 
                                                 type = interactions1[[name1]][[j]]$type, 
                                                 probability = interactions1[[name1]][[j]]$probability, 
                                                 support = interactions1[[name1]][[j]]$support, 
                                                 timestep = interactions1[[name1]][[j]]$timestep)
                if (type == 1) {
                    names(res[[index]])[[subindex]] <- paste(name1, "_", a_index, "_", "Activator", sep = "", collapse = "")
                    a_index = a_index + 1
                } else {
                    names(res[[index]])[[subindex]] <- paste(name1, "_", i_index, "_", "Inhibitor", sep = "", collapse = "")
                    i_index = i_index + 1
                }
                
            }
            names(res)[[index]] <- name1
        }
        
        ## if the name also in the second group
        if (name1 %in% names(interactions2)) {
            for (j in seq(interactions2[[name1]])) {
                inputgene2 <- genes2[interactions2[[name1]][[j]]$input]
                if (length(inputgene2) == 0) {
                    next
                }
                
                type <- interactions2[[name1]][[j]]$type
                expression <- interactions2[[name1]][[j]]$expression
                temp_exp <- paste(sort(splitExpression(expression, 2, TRUE)), sep = "", collapse = "")
                if (type == 1) {
                    if (temp_exp %in% unique_express_1) {
                        next
                    }
                    unique_express_1 <- c(unique_express_1, temp_exp)
                } else {
                    if (temp_exp %in% unique_express_0) {
                        next
                    }
                    
                    unique_express_0 <- c(unique_express_0, temp_exp)
                }
   

              subindex <- length(res[[name1]]) + 1
              res[[name1]][[subindex]] <- list(input = which(mergedgenes %in% inputgene2), 
                                               expression = interactions2[[name1]][[j]]$expression, 
                                               error = interactions2[[name1]][[j]]$error, 
                                               type = interactions2[[name1]][[j]]$type, 
                                               probability = interactions2[[name1]][[j]]$probability, 
                                               support = interactions2[[name1]][[j]]$support, 
                                               timestep = interactions2[[name1]][[j]]$timestep)
              if (type == 1) {
                  names(res[[name1]])[[subindex]] <- paste(name1, "_", a_index, "_", "Activator", sep = "", collapse = "")
                  a_index = a_index + 1
              } else {
                  names(res[[name1]])[[subindex]] <- paste(name1, "_", i_index, "_", "Inhibitor", sep = "", collapse = "")
                  i_index = i_index + 1
              }
   
            }
        }
    }
    
    names_2 <- names(interactions2)[!names(interactions2) %in% names(interactions1)]
    
    for (name2 in names_2) {
        unique_express_1 <- c()
        unique_express_0 <- c()

        if (!(name2 %in% names(res))) {
            index <- length(res) + 1
            res[[index]] <- list()
            for (j in seq(interactions2[[name2]])) {
                inputgene2 <- genes2[interactions2[[name2]][[j]]$input]
                type <- interactions2[[name2]][[j]]$type
                expression <- interactions2[[name2]][[j]]$expression
                temp_exp <- paste(sort(splitExpression(expression, 2, TRUE)), sep = "", collapse = "")
                if (length(inputgene2) == 0) {
                  next
                }

                if (type == 1) {
                    if (temp_exp %in% unique_express_1)
                        next

                    unique_express_1 <- c(unique_express_1, temp_exp)
                } else {
                    if (temp_exp %in% unique_express_0)
                        next

                    unique_express_0 <- c(unique_express_0, temp_exp)
                }

                subindex <- length(res[[index]]) + 1
                newinput <- which(mergedgenes %in% inputgene2)
                if (length(newinput) == 0) {
                  next
                }
                res[[index]][[subindex]] <- list(input = newinput, expression = interactions2[[name2]][[j]]$expression, error = interactions2[[name2]][[j]]$error,
                  type = interactions2[[name2]][[j]]$type, probability = interactions2[[name2]][[j]]$probability, support = interactions2[[name2]][[j]]$support,
                  timestep = interactions2[[name2]][[j]]$timestep)
                names(res[[index]])[[subindex]] <- names(interactions2[[name2]])[[j]]
            }
            names(res)[[index]] <- name2
        }
    }
    res
}

#' @export
mergeClusterNetworks <- function(clusteredFBNCube, threshold_error, maxFBNRules) {
    network_cores <- list()
    # validate network types
    if (!(inherits(clusteredFBNCube, "ClusteredFBNCube"))) 
        stop("clusteredFBNCube must be inherited from ClusteredFBNCube")
    
    len_cube <- length(clusteredFBNCube)
    i <- 1
    initial_cores <- clusteredFBNCube[[1]]$NetworkCores
    genes <- clusteredFBNCube[[1]]$Genes
    network_cores <- initial_cores
    
    i <- i + 1
    while (i <= len_cube) {
        gene_names <- names(network_cores)
        next_cores <- clusteredFBNCube[[i]]$NetworkCores
        genes2 <- clusteredFBNCube[[i]]$Genes
        gene_names_new <- names(next_cores)
        common_genes <- gene_names_new[gene_names_new %in% gene_names]
        un_common_genes <- gene_names_new[!gene_names_new %in% gene_names]
        genes <- unique(c(genes, genes2))
        combine_cores <- list(network_cores, next_cores[un_common_genes])
        combine_cores <- unlist(combine_cores,
                                recursive = FALSE)
        for (t_genes in common_genes){
            identities <- sapply(combine_cores[[t_genes]], function(entry)entry["identity.Identity"])
            len_next <- length(next_cores[[t_genes]])
            for (j in seq_len(len_next)){
                entry <- next_cores[[t_genes]][[j]]
                if (!entry["identity.Identity"] %in% identities) {
                    len_c <- length(combine_cores[[t_genes]]) + 1
                    combine_cores[[t_genes]][[len_c]] <- entry
                }
            }
        }
        network_cores <- combine_cores
        i <- i + 1
    }
    mineFBNNetworkWithCores(searchFBNNetworkCore = network_cores, 
                            genes = genes, 
                            threshold_error = threshold_error, 
                            maxFBNRules = maxFBNRules)
}

findAllInputGenes <- function(networkinteractions, genes) {
    # todo, the for has a problem and need to address carefully
    res <- c()
    for (i in seq(networkinteractions)) {
        for (j in seq(networkinteractions[[i]])) {
            res <- c(res, networkinteractions[[i]][[j]]$input)
        }
    }
    unique(genes[res])
}


findAllTargetGenes <- function(networkinteractions) {
    # todo, the for has a problem and need to address carefully
    res <- c()
    for(name in names(networkinteractions)) {
        if(length(networkinteractions[[name]]) > 0) {
            res <- c(res, name)
        }
    }
    unique(res)
}


#'filtering out non regulate genes
#'
#'@param networks The FBN networks that are going to be filtered
#'@return filtered networks
#' @export
filterNetworkConnections <- function(networks) {
    # need to check the index for all inputs
    res <- networks
    genes <- res$genes
    filterednetworks <- res$interactions[lapply(res$interactions, length) > 0]
    
    regulategenes <- names(filterednetworks)
    filteredinputgenes <- findAllInputGenes(filterednetworks, genes)
    
    mixedgenes <- unique(c(regulategenes, filteredinputgenes))
    
    extranetworks <- res$interactions[names(res$interactions) %in% mixedgenes]
    
    merges <- mergeInteraction(extranetworks, extranetworks, genes, genes, mixedgenes)
    
    if (length(res$fixed) < length(genes)) {
        res$fixed <- rep(-1, length(genes))
    }
    
    fixed1 <- res$fixed
    names(fixed1) <- genes
    
    newfixed <- rep(-1, length(mixedgenes))
    names(newfixed) <- mixedgenes
    newfixed[which(mixedgenes %in% names(fixed1))] <- fixed1[which(names(fixed1) %in% mixedgenes)]
    
    
    if (length(res$timedecay) < length(genes)) {
        res$timedecay <- rep(-1, length(genes))
    }
    
    timedecay1 <- res$timedecay
    names(timedecay1) <- genes
    
    newtimedecay <- rep(-1, length(mixedgenes))
    names(newtimedecay) <- mixedgenes
    newtimedecay[which(mixedgenes %in% names(timedecay1))] <- timedecay1[which(names(timedecay1) %in% mixedgenes)]
    
    res <- list(interactions = merges, genes = mixedgenes, fixed = newfixed, timedecay = newtimedecay)
    
    class(res) <- "FundamentalBooleanNetwork"
    return(res)
}

#'filtering networks with specific genes
#'
#'@param networks The FBN networks that are going to be filtered
#'@return filtered networks
#' @export
filterNetworkConnectionsByGenes <- function(networks, genelist = c(), exclusive = TRUE, expand = FALSE) {
    # individualFilter2<-filterNetworkConnections(individualNetworks[[2]]) FBNNetwork.Graph(individualFilter2)
    # individualFilter2<-filterNetworkConnectionsByGenes(individualFilter2,c('OR7A5','GFOD1','ANGPT2','LCN2','MCTP1','PNMT', 'SLC22A23'))
    # FBNNetwork.Graph(individualFilter2) individualFilter2 originalgenes2<-individualNetworks[[2]]$genes
    # diffgenes2<-originalgenes2[!originalgenes2%in%individualFilter2$genes] diffgenes2
    if (length(genelist) == 0) {
        stop("The genelist is empty")
    }
    
    res <- networks
    genes <- res$genes
    if(exclusive) {
        filterednetworks <- res$interactions[!names(res$interactions) %in% genelist]
    } else {
        filterednetworks <- res$interactions[names(res$interactions) %in% genelist]
    }

    regulategenes <- names(filterednetworks)

    filteredinputgenes <- findAllInputGenes(filterednetworks, genes)
    
    if(!expand) {
        extranetworks <- filterednetworks
        mixedgenes <- filteredinputgenes
    }else {
        mixedgenes <- unique(c(regulategenes, filteredinputgenes))
        extranetworks <- res$interactions[names(res$interactions) %in% mixedgenes]
        mixedgenes <- findAllInputGenes(extranetworks, genes)
    }
    
    merges <- mergeInteraction(extranetworks, extranetworks, genes, genes, mixedgenes)
    
    if (length(res$fixed) < length(genes)) {
        res$fixed <- rep(-1, length(genes))
    }
    
    fixed1 <- res$fixed
    names(fixed1) <- genes
    
    newfixed <- rep(-1, length(mixedgenes))
    names(newfixed) <- mixedgenes
    newfixed[which(mixedgenes %in% names(fixed1))] <- fixed1[which(names(fixed1) %in% mixedgenes)]
    
    
    if (length(res$timedecay) < length(genes)) {
        res$timedecay <- rep(-1, length(genes))
    }
    
    timedecay1 <- res$timedecay
    names(timedecay1) <- genes
    
    newtimedecay <- rep(-1, length(mixedgenes))
    names(newtimedecay) <- mixedgenes
    newtimedecay[which(mixedgenes %in% names(timedecay1))] <- timedecay1[which(names(timedecay1) %in% mixedgenes)]
    
    res <- list(interactions = merges, genes = mixedgenes, fixed = newfixed, timedecay = newtimedecay)
    
    class(res) <- "FundamentalBooleanNetwork"
    filterNetworkConnections(res)
}



#'internal
#' @noRd
findAllBackwardRelatedGenes <- function(networks, 
                                        target_gene, 
                                        regulationType = NULL, 
                                        target_type = NULL,  
                                        maxDeep = 1, 
                                        next_level_mix_type = FALSE) {
    
    prepare_network <- filterNetworkConnectionsByGenes(networks,target_gene, exclusive = FALSE, expand = FALSE)
    genes <- names(prepare_network$interactions)
    remove_index_gene <- list()

    len <- length(genes)
    for (i in seq_len(len)) {
        interactions <- prepare_network$interactions[[i]]  #gene level
        len_i <- length(interactions)
        remove_index_function <- c()
        for (j in seq_len(len_i)) {
            interaction <- interactions[[j]]  #function level
            type <- interaction$type

            if (!is.null(regulationType) && type != regulationType) {
                remove_index_function <- c(remove_index_function, j)
                next
            }
        }
        
        gene_name <- names(prepare_network$interactions)[[i]]
        remove_index_gene[[gene_name]] <- unique(remove_index_function)
    }
    
    len_new <- length(remove_index_gene)
    for (i in seq_len(len_new)) {
        name <- names(remove_index_gene)[[i]]
        remove_indexes <- remove_index_gene[[i]]
        prepare_network$interactions[[name]] <- prepare_network$interactions[[name]][-c(remove_indexes)]
    }
   
    prepare_network <- filterNetworkConnections(prepare_network)
    if (maxDeep > 1 && length(prepare_network$interactions) > 0) {
        # cond1 <- sapply(prepare_network$interactions, function(interaction) length(interaction) > 0)
        # filteredNetworkInteractions <- prepare_network$interactions[cond1]
        # expand_genes <- names(filteredNetworkInteractions)
        expand_genes <- prepare_network$genes[!prepare_network$genes%in%target_gene]
        maxDeep <- maxDeep - 1

        for (this_target_gene in expand_genes) {
          if (next_level_mix_type) {
            
          } else {
            
          }
            new_nextwork <- findAllBackwardRelatedGenes(networks = networks, 
                                                        target_gene = this_target_gene, 
                                                        regulationType = regulationType,
                                                        target_type = target_type,
                                                        maxDeep = maxDeep)
            prepare_network <- filterNetworkConnections(mergeNetwork(prepare_network, new_nextwork))
        }
    }
    prepare_network
}

#' Find Backward Network By Genes
#' 
#' @param networks The Fundamental Boolean Network
#' @param genelist The target genes
#' @export
findBackwardRelatedNetworkByGenes <- function(networks,
                                       target_gene_list = c(), 
                                       regulationType = NULL, 
                                       maxDeep = 1,
                                       next_level_mix_type = FALSE) {
    
    
    if (length(target_gene_list) == 0) {
        stop("The genelist is empty")
    }
    
    prepare_network <- NULL
    for (target_gene in target_gene_list) {
        if (regulationType == 0) {
            target_type <- 0
        }else{
            target_type <- 1
        }
        new_nextwork <- findAllBackwardRelatedGenes(networks = networks, 
                                                    target_gene = target_gene, 
                                                    regulationType = regulationType, 
                                                    target_type = target_type,
                                                    maxDeep = maxDeep,
                                                    next_level_mix_type = next_level_mix_type)
        if (is.null(prepare_network)) {
            prepare_network <- new_nextwork
        } else {
            prepare_network <- mergeNetwork(prepare_network, new_nextwork)
        }
        prepare_network <- filterNetworkConnections(prepare_network)
    }
    prepare_network
}

#'internal
#' @noRd
findAllForwardRelatedGenes <- function(networks, 
                                       target_gene, 
                                       regulationType = NULL, 
                                       target_type = NULL, 
                                       main_target_gene = NULL,
                                       main_target_type = NULL,
                                       maxDeep = 1,
                                       next_level_mix_type = next_level_mix_type) {
    genes <- networks$genes
    prepare_network <- networks
    remove_index_gene <- list()
    len <- length(genes)
    for (i in seq_len(len)) {
        interactions <- networks$interactions[[i]]  #gene level
        len_i <- length(interactions)
        remove_index_function <- c()
        for (j in seq_len(len_i)) {
            interaction <- interactions[[j]]  #function level
            input_Genes <- genes[interaction$input]
            type <- interaction$type
            expression <- interaction$expression

            
            if (!target_gene %in% input_Genes) {
                remove_index_function <- c(remove_index_function, j)
                next
            }
            
            if (!is.null(regulationType) && type != regulationType) {
                remove_index_function <- c(remove_index_function, j)
                next
            }
            thistarget_type <- 1
            
            if (stringr::str_detect(expression, paste0("!", target_gene))) {
                thistarget_type <- 0
            }
            
            thismaintarget_type <- 1
            if (stringr::str_detect(expression, paste0("!", main_target_gene))) {
                thismaintarget_type <- 0
            }

            if(main_target_gene != target_gene && main_target_gene %in% input_Genes && thismaintarget_type != main_target_type) {
                remove_index_function <- c(remove_index_function, j)
                next
            }
            
            if (!is.null(target_type) && thistarget_type != target_type) {
                remove_index_function <- c(remove_index_function, j)
                next
            }
        }
        
        gene_name <- names(networks$interactions)[[i]]
        remove_index_gene[[gene_name]] <- unique(remove_index_function)
    }
    
    len_new <- length(remove_index_gene)
    for (i in seq_len(len_new)) {
        name <- names(remove_index_gene)[[i]]
        remove_indexes <- remove_index_gene[[i]]
        prepare_network$interactions[[name]] <- prepare_network$interactions[[name]][-c(remove_indexes)]
    }
    
    prepare_network <- filterNetworkConnections(prepare_network)
    if (maxDeep > 1 && length(prepare_network$interactions) > 0) {
        cond1 <- sapply(prepare_network$interactions, function(interaction) length(interaction) > 0)
        filteredNetworkInteractions <- prepare_network$interactions[cond1]
        #expand_genes <- names(filteredNetworkInteractions)
        #expand_genes <- prepare_network$genes[!prepare_network$genes%in%target_gene]
        expand_genes <- findAllTargetGenes(filteredNetworkInteractions)
        expand_genes <- expand_genes[!expand_genes %in% target_gene]
        maxDeep <- maxDeep - 1
        if (regulationType == 1) {
            target_type <- 1
        }else {
            target_type <- 0
        }
        
        for (this_target_gene in expand_genes) {
            if (next_level_mix_type) {
                new_nextwork <- findAllForwardRelatedGenes(networks = networks, 
                                                           target_gene = this_target_gene, 
                                                           regulationType = 1, 
                                                           target_type = target_type, 
                                                           main_target_gene = main_target_gene,
                                                           main_target_type = main_target_type,
                                                           maxDeep = maxDeep,
                                                           next_level_mix_type = next_level_mix_type)
                prepare_network <- filterNetworkConnections(mergeNetwork(prepare_network, new_nextwork))
                
                new_nextwork <- findAllForwardRelatedGenes(networks = networks, 
                                                           target_gene = this_target_gene, 
                                                           regulationType = 0, 
                                                           target_type = target_type, 
                                                           main_target_gene = main_target_gene,
                                                           main_target_type = main_target_type,
                                                           maxDeep = maxDeep,
                                                           next_level_mix_type = next_level_mix_type)
                prepare_network <- filterNetworkConnections(mergeNetwork(prepare_network, new_nextwork))
            } else {
                new_nextwork <- findAllForwardRelatedGenes(networks = networks, 
                                                           target_gene = this_target_gene, 
                                                           regulationType = regulationType, 
                                                           target_type = target_type, 
                                                           main_target_gene = main_target_gene,
                                                           main_target_type = main_target_type,
                                                           maxDeep = maxDeep,
                                                           next_level_mix_type = next_level_mix_type)
                prepare_network <- filterNetworkConnections(mergeNetwork(prepare_network, new_nextwork))
            }
        }
    }
    prepare_network
}

#' Find Forward related Network By Genes
#' 
#' @param networks The Fundamental Boolean Network
#' @param target_gene_list The target genes
#' @param regulationType The type of regulation either in 1 (Up regulation) or 0 (down regulation)
#' @param target_type The boolean type of target genes either in 1 (Activation) or 0 (Inhibition)
#' @export
findForwardRelatedNetworkByGenes <- function(networks, 
                                             target_gene_list = c(), 
                                             regulationType = NULL, 
                                             target_type = NULL, 
                                             maxDeep = 1,
                                             next_level_mix_type = FALSE) {
    if (length(target_gene_list) == 0) {
        stop("The genelist is empty")
    }
    
    prepare_network <- NULL
    for (target_gene in target_gene_list) {
        new_nextwork <- findAllForwardRelatedGenes(networks = networks, 
                                                   target_gene = target_gene, 
                                                   regulationType = regulationType, 
                                                   target_type = target_type, 
                                                   main_target_gene = target_gene,
                                                   main_target_type = target_type,
                                                   maxDeep = maxDeep,
                                                   next_level_mix_type = next_level_mix_type)
        if (is.null(prepare_network)) {
            prepare_network <- new_nextwork
        } else {
            prepare_network <- mergeNetwork(prepare_network, new_nextwork)
        }
        prepare_network <- filterNetworkConnections(prepare_network)
    }
    prepare_network
}
# use sink to output result to text file sink('analysis-output.txt') attactor7 sink()

# new feature to plot all type of networks in one go, not finish yet.
#' A function to plot the Fundamental boolean networks
#' 
#' @param FBNNetwork The FBN networks
#' @param target_genes A list of target genes
#' @param type The visualization type
#' @param expand_level The level of expandation
#' @param output_network Optional, if TRUE, then output the FBN Network
#' @param timeseries_matrix The timeseries matrix that used to plot dynamic 
#' networks.
#' @param start_time_point The start time point
#' @param end_time_point The end time point
#' @param target_time_point The target time point
#' @export
plotNetwork <- function(FBNNetwork, target_genes = c(), type = c("static", "staticSlice", "dynamic", "forward_1a", "forward_2a", "forward_3a", "forward_4a", "backward_1a", 
    "backward_2a", "forward_1b", "forward_2b", "forward_3b", "forward_4b", "backward_1b", "backward_2b"), expand_level = 2, output_network = FALSE, timeseries_matrix = NULL, 
    start_time_point = 1, end_time_point = 1, target_time_point = 1) {
    if (!(inherits(FBNNetwork, "FundamentalBooleanNetwork"))) 
        stop("Network must be inherited from FundamentalBooleanNetwork")
    
    switch(type, static = {
        if (length(target_genes) > 0) {
            FBNNetwork <- filterNetworkConnectionsByGenes(FBNNetwork, genelist = target_genes, exclusive = FALSE, expand = FALSE)
        }
        FBNNetwork.Graph(FBNNetwork, type = "static")
    }, staticSlice = {
        if (length(target_genes) > 0) {
            FBNNetwork <- filterNetworkConnectionsByGenes(FBNNetwork, genelist = target_genes, exclusive = FALSE, expand = FALSE)
        }
        FBNNetwork.Graph(FBNNetwork, type = "staticSlice", timeseriesMatrix = timeseries_matrix, toTimePoint = target_time_point)
    }, dynamic = {
        if (length(target_genes) > 0) {
            FBNNetwork <- filterNetworkConnectionsByGenes(FBNNetwork, genelist = target_genes, exclusive = FALSE, expand = FALSE)
        }
        FBNNetwork.Graph(FBNNetwork, type = "dynamic", timeseriesMatrix = timeseries_matrix, fromTimePoint = start_time_point, toTimePoint = end_time_point)
    }, forward_1a = {
        FBNNetwork <- findForwardRelatedNetworkByGenes(FBNNetwork, target_gene_list = target_genes, regulationType = 0, target_type = 0, maxDeep = 2, next_level_mix_type = FALSE)
        FBNNetwork.Graph(FBNNetwork)
    }, forward_2a = {
        FBNNetwork <- findForwardRelatedNetworkByGenes(FBNNetwork, target_gene_list = target_genes, regulationType = 0, target_type = 1, maxDeep = 2, next_level_mix_type = FALSE)
        FBNNetwork.Graph(FBNNetwork)
    }, forward_3a = {
        FBNNetwork <- findForwardRelatedNetworkByGenes(FBNNetwork, target_gene_list = target_genes, regulationType = 1, target_type = 0, maxDeep = 2, next_level_mix_type = FALSE)
        FBNNetwork.Graph(FBNNetwork)
    }, forward_4a = {
        FBNNetwork <- findForwardRelatedNetworkByGenes(FBNNetwork, target_gene_list = target_genes, regulationType = 1, target_type = 1, maxDeep = 2, next_level_mix_type = FALSE)
        FBNNetwork.Graph(FBNNetwork)
    }, forward_1b = {
        FBNNetwork <- findForwardRelatedNetworkByGenes(FBNNetwork, target_gene_list = target_genes, regulationType = 0, target_type = 0, maxDeep = 2, next_level_mix_type = TRUE)
        FBNNetwork.Graph(FBNNetwork)
    }, forward_2b = {
        FBNNetwork <- findForwardRelatedNetworkByGenes(FBNNetwork, target_gene_list = target_genes, regulationType = 0, target_type = 1, maxDeep = 2, next_level_mix_type = TRUE)
        FBNNetwork.Graph(FBNNetwork)
    }, forward_3b = {
        FBNNetwork <- findForwardRelatedNetworkByGenes(FBNNetwork, target_gene_list = target_genes, regulationType = 1, target_type = 0, maxDeep = 2, next_level_mix_type = TRUE)
        FBNNetwork.Graph(FBNNetwork)
    }, forward_4b = {
        FBNNetwork <- findForwardRelatedNetworkByGenes(FBNNetwork, target_gene_list = target_genes, regulationType = 1, target_type = 1, maxDeep = 2, next_level_mix_type = TRUE)
        FBNNetwork.Graph(FBNNetwork)
    }, backward_1a = {
        FBNNetwork <- findBackwardRelatedNetworkByGenes(FBNNetwork, target_gene_list = target_genes, regulationType = 0, maxDeep = 2, next_level_mix_type = FALSE)
        FBNNetwork.Graph(FBNNetwork)
    }, backward_2a = {
        FBNNetwork <- findBackwardRelatedNetworkByGenes(FBNNetwork, target_gene_list = target_genes, regulationType = 1, maxDeep = 2, next_level_mix_type = FALSE)
        FBNNetwork.Graph(FBNNetwork)
    }, backward_1b = {
        FBNNetwork <- findBackwardRelatedNetworkByGenes(FBNNetwork, target_gene_list = target_genes, regulationType = 0, maxDeep = 2, next_level_mix_type = TRUE)
        FBNNetwork.Graph(FBNNetwork)
    }, backward_2b = {
        FBNNetwork <- findBackwardRelatedNetworkByGenes(FBNNetwork, target_genes, regulationType = 1, maxDeep = 2, next_level_mix_type = TRUE)
        FBNNetwork.Graph(FBNNetwork)
    })
    if (output_network) {
        return(FBNNetwork)
    }
    NULL
}

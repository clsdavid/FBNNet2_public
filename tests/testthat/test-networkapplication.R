context("networkapplication")
setupdata <- function() {
  with(ExampleNetwork, {
    initialStates <- generateAllCombinationBinary(ExampleNetwork$genes)
    trainingseries <- genereateBoolNetTimeseries(ExampleNetwork, initialStates, 43, type = "synchronous")
    
    cube <- constructFBNCube(ExampleNetwork$genes, ExampleNetwork$genes, trainingseries, 5, 1, FALSE)
    NETWORK2 <- mineFBNNetwork(cube, ExampleNetwork$genes)
    
    return(list(network = NETWORK2, initialStates = initialStates, timeseries = trainingseries))
  })
}
network <- setupdata()$network
describe("test  filterNetworkConnectionsByGenes", {
  
  it("target = Gene1, with exclusive = TRUe, Expand = FALSE", {
     tt <- filterNetworkConnectionsByGenes(network, "Gene1", exclusive = TRUE, expand = FALSE)
     expect_true(is.null(tt$interactions[["Gene1"]]))
  })
  
  it("target = Gene1, with exclusive = FALSE, Expand = FALSE", {
    tt <- filterNetworkConnectionsByGenes(network, "Gene1", exclusive = FALSE, expand = FALSE)
    expect_true(length(tt$interactions)==1)
  })
})

describe("test  findForwardRelatedNetworkByGenes", {
  
  it("target = Gene1, with 0,0, 1", {
    tt <- findForwardRelatedNetworkByGenes(network, "Gene1", 0,0, 1)
    expect_true(length(tt$interactions) == 2)
    expect_true(length(tt$interactions[["Gene1"]]) == 1)
    expect_true(length(tt$interactions[["Gene2"]]) == 1)

    expect_true(tt$interactions[["Gene1"]]$Gene1_1_Inhibitor$expression == "!Gene1")
    expect_true(tt$interactions[["Gene2"]]$Gene2_1_Inhibitor$expression == "!Gene1")
  })
  
  it("target = Gene1, with exclusive = FALSE, Expand = FALSE", {
    tt <-  findForwardRelatedNetworkByGenes(network, "Gene1", 0,1, 1)
    expect_true(length(tt$interactions)==3)
    expect_true(tt$interactions[["Gene4"]]$Gene4_1_Inhibitor$expression == "Gene1&Gene5")
  })
})
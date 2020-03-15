context("networkapplication")
setupdata <- function() {
  network <- BoolNet::loadNetwork("D:\\Dropbox/Dropbox/FBNNet/tempdata/Example.txt")
  # trainingseries<-FBNDataReduction(generateTimeSeries(network,1000,43,noiseLevel=0.0))
  initialStates <- generateAllCombinationBinary(network$genes)
  trainingseries <- genereateBoolNetTimeseries(network, initialStates, 43, type = "synchronous")
  
  cube <- constructFBNCube(network$genes, network$genes, trainingseries, 5, 1, FALSE)
  NETWORK2 <- mineFBNNetwork(cube, network$genes)
  
  return(list(network = NETWORK2, initialStates = initialStates, timeseries = trainingseries, cube = cube))
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
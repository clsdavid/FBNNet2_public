context("constructFBNCube")

describe("test with the example gene list", {
    data(ExampleNetwork)
    network <- ExampleNetwork
    print(network)
    initialStates <- generateAllCombinationBinary(network$genes)
    trainingseries <- genereateBoolNetTimeseries(network, initialStates, 43, type = "synchronous")
    it("Test", {
        cube <- expect_error(constructFBNCube(network$genes, network$genes, trainingseries, 5, 1, FALSE), NA)
        
        expect_true(!is.null(cube))
        NETWORK2 <- expect_error(mineFBNNetwork(cube, network$genes), NA)
        print(NETWORK2)
        expect_error(FBNNetwork.Graph(NETWORK2), NA)
        # training dataset
        resultfile <- expect_error(reconstructTimeseries(NETWORK2, initialStates, type = "synchronous", maxTimepoints = 43, useParallel = FALSE), NA)
        
        similarreport <- expect_error(generateSimilaryReport(trainingseries, resultfile), NA)
        expect_equal(similarreport$ErrorRate, 0)
        expect_equal(similarreport$AccurateRate, 1)
        expect_equal(similarreport$MissMatchedRate, 0)
        expect_equal(similarreport$PerfectMatchedRate, 1)
    })
})

describe("test with the sub cube", {
    data(ExampleNetwork)
    network <- ExampleNetwork
    print(network)
    initialStates <- generateAllCombinationBinary(network$genes)
    trainingseries <- genereateBoolNetTimeseries(network, initialStates, 43, type = "synchronous")
    it("Test", {
        NETWORK2 <- generateFBMNetwork(timeseries_data = trainingseries, 
                                                  maxK = 4, 
                                                  max_deep_temporal = 1, 
                                                  useParallel = FALSE)
        
        print(NETWORK2)
        expect_error(FBNNetwork.Graph(NETWORK2), NA)
        # training dataset
        resultfile <- expect_error(reconstructTimeseries(NETWORK2, initialStates, type = "synchronous", maxTimepoints = 43, useParallel = FALSE), NA)
        
        similarreport <- expect_error(generateSimilaryReport(trainingseries, resultfile), NA)
        expect_equal(similarreport$ErrorRate, 0)
        expect_equal(similarreport$AccurateRate, 1)
        expect_equal(similarreport$MissMatchedRate, 0)
        expect_equal(similarreport$PerfectMatchedRate, 1)
    })
})

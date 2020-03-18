context("attractor")
require(BoolNet)
setupdata <- function() {
    with(ExampleNetwork, {
        initialStates <- generateAllCombinationBinary(ExampleNetwork$genes)
        trainingseries <- genereateBoolNetTimeseries(ExampleNetwork, initialStates, 43, type = "synchronous")
        
        cube <- constructFBNCube(ExampleNetwork$genes, ExampleNetwork$genes, trainingseries, 5, 1, FALSE)
        NETWORK2 <- mineFBNNetwork(cube, ExampleNetwork$genes)
        
        return(list(network = NETWORK2, initialStates = initialStates, timeseries = trainingseries))
    })

}

describe("run synchronous should succeed", {
    test_info <- setup(setupdata())
    it("reconstruct timeseries", {
        resultfile <- expect_error(reconstructTimeseries(test_info$network,
                                                         test_info$initialStates,
                                                         type = "synchronous",
                                                         maxTimepoints = 43,
                                                         useParallel = FALSE), 
            NA)
        
        similarreport <- expect_error(generateSimilaryReport(test_info$timeseries, resultfile), NA)
        expect_equal(similarreport$ErrorRate, 0)
        expect_equal(similarreport$AccurateRate, 1)
        expect_equal(similarreport$MissMatchedRate, 0)
        expect_equal(similarreport$PerfectMatchedRate, 1)
        
        
        attractor <- expect_error(searchForAttractors(test_info$network, test_info$initialStates, test_info$network$genes), NA)
        expect_error(print(attractor), NA)
        # display the dynamic trajectory of the attactor 2
        expect_error(FBNNetwork.Graph.DrawAttractor(test_info$network, attractor, 2), NA)
    })
})

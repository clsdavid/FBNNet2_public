context("constructFBNCube")
setupdata <- function() {
    genesInput <- c("CycD", "p27", "CycE", "E2F")
    testseries <- list()
    testseries[[1]] <- matrix(c(1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0), nrow = 4, ncol = 6, byrow = FALSE, dimnames = list(genesInput, 
        c("1", "2", "3", "4", "5", "6")))
    testseries[[2]] <- matrix(c(1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1), nrow = 4, ncol = 6, byrow = FALSE, dimnames = list(genesInput, 
        c("1", "2", "3", "4", "5", "6")))
    testseries[[3]] <- matrix(c(1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0), nrow = 4, ncol = 6, byrow = FALSE, dimnames = list(genesInput, 
        c("1", "2", "3", "4", "5", "6")))
    
    
    # Test 1
    getCurrentStates <- list()
    getpreviousStates <- list()
    getCurrentStates_c <- list()
    getpreviousStates_c <- list()
    index <- 3
    while (index > 0) {
        
        getCurrentStates[[index]] <- extractGeneStateFromTimeSeriesCube(testseries, index)
        getpreviousStates[[index]] <- getCurrentStates[[index]]
        getCurrentStates_c[[index]] <- getCurrentStates[[index]]
        getpreviousStates_c[[index]] <- getCurrentStates[[index]]
        index <- index - 1
    }
    
    mainParameters <- new.env(parent = globalenv())
    mainParameters$currentStates <- getCurrentStates
    mainParameters$previousStates <- getpreviousStates
    mainParameters$currentStates_c <- getCurrentStates_c
    mainParameters$previousStates_c <- getpreviousStates_c
    mainParameters$timeseries <- testseries
    return(mainParameters)
}

# describe('uccessful',{ mainParameters<-setup(setupdata()) genesInput<-c('CycD','p27','CycE','E2F') it('run a normal set',{
# #constructFBNCube<-function(targetgenes,allgenes,timeseriesCube,maxK=5,temporal=1,useParallel=FALSE,optimizedTemoral=FALSE)
# cube<-expect_error(constructFBNCube(targetgenes=genesInput,allgenes=genesInput,timeseriesCube=mainParameters$timeseries,maxK=5,temporal=1,useParallel=FALSE,optimizedTemoral=FALSE),NA)
# }) })

generate_test_example_file <- function() {
    sink("example.bn")
    cat("targets, factors\n")
    cat("Gene1, Gene1\n")
    cat("Gene2, Gene1 & Gene5 & !Gene4\n")
    cat("Gene3, Gene3\n")
    cat("Gene4, Gene3 & !(Gene1 & Gene5)\n")
    cat("Gene5, !Gene2\n")
    sink()
}
describe("test with the example gene list", {
    generate_test_example_file()
    network <- BoolNet::loadNetwork("example.bn")
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
    generate_test_example_file()
    network <- BoolNet::loadNetwork("example.bn")
    print(network)
    initialStates <- generateAllCombinationBinary(network$genes)
    trainingseries <- genereateBoolNetTimeseries(network, initialStates, 43, type = "synchronous")
    it("Test", {
        NETWORK2 <- generateFBMNetwork(timeseries_data = trainingseries, 
                                                  maxK = 4, 
                                                  max_deep_temporal = 1, 
                                                  useParallel = FALSE,
                                                  maxGenesForSingleCube = 1)
        
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

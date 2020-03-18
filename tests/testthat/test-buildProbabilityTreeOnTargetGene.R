context("buildProbabilityTreeOnTargetGene")
require(BoolNet)
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

generate_test_example_file <- function() {
    require(BoolNet)
    sink("example.bn")
    cat("targets, factors\n")
    cat("Gene1, Gene1\n")
    cat("Gene2, Gene1 & Gene5 & !Gene4\n")
    cat("Gene3, Gene3\n")
    cat("Gene4, Gene3 & !(Gene1 & Gene5)\n")
    cat("Gene5, !Gene2\n")
    sink()
    
    network <- BoolNet::loadNetwork("example.bn")
    print(network)
    initialStates <- generateAllCombinationBinary(network$genes)
    trainingseries <- genereateBoolNetTimeseries(network, initialStates, 43, type = "synchronous")
    
    getCurrentStates <- list()
    getpreviousStates <- list()
    getCurrentStates_c <- list()
    getpreviousStates_c <- list()
    index <- 1
    while (index > 0) {
        
        getCurrentStates[[index]] <- extractGeneStateFromTimeSeriesCube(trainingseries, index)
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
    mainParameters$timeseries <- trainingseries
    return(mainParameters)
}
describe("Test with example genes", {
    mainParameters <- generate_test_example_file()
    
    it("test", {
        first <- mainParameters$timeseries[[1]]
        genes <- rownames(first)
        cube <- lapply(genes, function(gene, mainParameters, genesInput) {
            print(paste0("Test: target=", gene, " condition genes=", genesInput, collapse = ","))
            cube <- expect_error(buildProbabilityTreeOnTargetGene(gene, mainParameters, genesInput, NULL, NULL, 4, 1), NA)
            return(cube)
        }, mainParameters, genes)
        
    })
})

describe("CycD Successful", {
    mainParameters <- setup(setupdata())
    genesInput <- c("CycD", "p27", "CycE", "E2F")

    it("Run with temporal =1 with optimizedTemoral=T, successful", {
        cube <- expect_error(buildProbabilityTreeOnTargetGene("CycD", mainParameters, genesInput, NULL, NULL, 4, 1), NA)
    })
    
    it("Run with temporal =2 with optimizedTemoral=T, successful", {
        cube <- expect_error(buildProbabilityTreeOnTargetGene("CycD", mainParameters, genesInput, NULL, NULL, 4, 2), NA)
    })
    
    it("Run with temporal =3 with optimizedTemoral=T, successful", {
        cube <- expect_error(buildProbabilityTreeOnTargetGene("CycD", mainParameters, genesInput, NULL, NULL, 4, 3), NA)
    })
    
})

describe("p27 Successful", {
    mainParameters <- setup(setupdata())
    genesInput <- c("CycD", "p27", "CycE", "E2F")
    it("Run with temporal =1 with optimizedTemoral=T, successful", {
        cube <- expect_error(buildProbabilityTreeOnTargetGene("p27", mainParameters, genesInput, NULL, NULL, 4, 1), NA)
    })
    
    it("Run with temporal =2 with optimizedTemoral=T, successful", {
        cube <- expect_error(buildProbabilityTreeOnTargetGene("p27", mainParameters, genesInput, NULL, NULL, 4, 2), NA)
    })
    
    it("Run with temporal =3 with optimizedTemoral=T, successful", {
        cube <- expect_error(buildProbabilityTreeOnTargetGene("p27", mainParameters, genesInput, NULL, NULL, 4, 3), NA)
    })
    
})

describe("CycE Successful", {
    mainParameters <- setup(setupdata())
    genesInput <- c("CycD", "p27", "CycE", "E2F")
    it("Run with temporal =1 with optimizedTemoral=T, successful", {
        cube <- expect_error(buildProbabilityTreeOnTargetGene("CycE", mainParameters, genesInput, NULL, NULL, 4, 1), NA)
    })
    
    it("Run with temporal =2 with optimizedTemoral=T, successful", {
        cube <- expect_error(buildProbabilityTreeOnTargetGene("CycE", mainParameters, genesInput, NULL, NULL, 4, 2), NA)
    })
    
    it("Run with temporal =3 with optimizedTemoral=T, successful", {
        cube <- expect_error(buildProbabilityTreeOnTargetGene("CycE", mainParameters, genesInput, NULL, NULL, 4, 3), NA)
    })
    
})

describe("E2F Successful", {
    mainParameters <- setup(setupdata())
    genesInput <- c("CycD", "p27", "CycE", "E2F")

    it("Run with temporal =1 with optimizedTemoral=T, successful", {
        cube <- expect_error(buildProbabilityTreeOnTargetGene("E2F", mainParameters, genesInput, NULL, NULL, 4, 1), NA)
    })
    
    it("Run with temporal =2 with optimizedTemoral=T, successful", {
        cube <- expect_error(buildProbabilityTreeOnTargetGene("E2F", mainParameters, genesInput, NULL, NULL, 4, 2), NA)
    })
    
    it("Run with temporal =3 with optimizedTemoral=T, successful", {
        cube <- expect_error(buildProbabilityTreeOnTargetGene("E2F", mainParameters, genesInput, NULL, NULL, 4, 3), NA)
    })
    
})

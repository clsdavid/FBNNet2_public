context("getgeneprababilities_advanced")

setupdata <- function() {
    genesInput <- c("CycD", "p27", "CycE", "E2F")
    testseries <- list()
    testseries[[1]] <- matrix(c(1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1), nrow = 4, ncol = 4, byrow = FALSE, dimnames <- list(genesInput, c("1", "2", "3", 
        "4")))
    testseries[[2]] <- matrix(c(1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1), nrow = 4, ncol = 4, byrow = FALSE, dimnames <- list(genesInput, c("1", "2", "3", 
        "4")))
    testseries[[3]] <- matrix(c(1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1), nrow = 4, ncol = 4, byrow = FALSE, dimnames <- list(genesInput, c("1", "2", "3", 
        "4")))
    
    
    # Test 1
    getCurrentStates <- list()
    getpreviousStates <- list()
    getCurrentStates_c <- list()
    getpreviousStates_c <- list()
    index <- 2
    while (index > 0) {
        
        getCurrentStates[[index]] <- extractGeneStateFromTimeSeriesCube(testseries, index)
        getpreviousStates[[index]] <- getCurrentStates[[index]]
        getCurrentStates_c[[index]] <- getCurrentStates[[index]]
        getpreviousStates_c[[index]] <- getCurrentStates[[index]]
        index <- index - 1
    }
    
    total_timepoints <- sum(vapply(testseries, function(timeshet) ncol(timeshet),
                                   integer(1)));
    total_samples <- length(testseries)
    all_gene_names <- rownames(testseries[[1]])
    
    
    mainParameters <- new.env(parent = globalenv())
    mainParameters$currentStates <- getCurrentStates
    mainParameters$previousStates <- getpreviousStates
    mainParameters$currentStates_c <- getCurrentStates_c
    mainParameters$previousStates_c <- getpreviousStates_c
    mainParameters$total_samples <- total_samples
    mainParameters$all_gene_names <- all_gene_names
    mainParameters$total_timepoints <- total_timepoints
    
    return(mainParameters)
}

describe("get gene probabilities temporal should succeed", {
    mainParameters <- setup(setupdata())
    it("getGenePrababilities(mainParameters,NULL,\"CycD\",\"p27\",1)", {
        basic_measures <- getGenePrababilities_basic(mainParameters, NULL, "CycD", "p27", 1, NULL)
        probability <- getGenePrababilities_advanced(basic_measures)$getBestFitP
        expect_equal(round(probability$TT, 3), 0.5)
        expect_equal(round(probability$FT, 3), 0.5)
        expect_equal(round(probability$TF, 3), 1)
        expect_equal(round(probability$FF, 3), 0)
        
        # test counter
        expect_equal(round(probability$TT_c, 3), 0.833)
        expect_equal(round(probability$FT_c, 3), 0.167)
        expect_equal(round(probability$TF_c, 3), 0.667)
        expect_equal(round(probability$FF_c, 3), 0.333)
    })
    
    it("getGenePrababilities(mainParameters,NULL,\"p27\",\"CycD\",1)", {
        basic_measures <- getGenePrababilities_basic(mainParameters, NULL, "p27", "CycD", 1, NULL)
        probability <- getGenePrababilities_advanced(basic_measures)$getBestFitP
        expect_equal(round(probability$TT, 3), 0.833)
        expect_equal(round(probability$FT, 3), 0.167)
        expect_equal(round(probability$TF, 3), 0.667)
        expect_equal(round(probability$FF, 3), 0.333)
    })
    
    it("getGenePrababilities(mainParameters,NULL,\"p27\",\"CycE\",1)", {
        basic_measures <- getGenePrababilities_basic(mainParameters, NULL, "p27", "CycE", 1, NULL)
        probability <- getGenePrababilities_advanced(basic_measures)$getBestFitP
        expect_equal(round(probability$TT, 3), 0.5)
        expect_equal(round(probability$FT, 3), 0.5)
        expect_equal(round(probability$TF, 3), 0.857)
        expect_equal(round(probability$FF, 3), 0.143)
    })
    
    it("getGenePrababilities(mainParameters,list(\"p27\"=1,\"CycD\"=0),\"CycE\",\"E2F\",1)", {
        basic_measures <- getGenePrababilities_basic(mainParameters, list(p27 = 1, CycD = 0), "CycE", "E2F", 1, NULL)
        probability <- getGenePrababilities_advanced(basic_measures)$getBestFitP
        expect_equal(round(probability$TT, 3), 0.333)
        expect_equal(round(probability$FT, 3), 0.667)
        expect_equal(round(probability$TF, 3), 0)
        expect_equal(round(probability$FF, 3), 0)
        
    })
})

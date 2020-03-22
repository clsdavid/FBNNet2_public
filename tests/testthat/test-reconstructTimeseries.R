context("reconstructTimeseries")
setupdata <- function() {
  with(ExampleNetwork, {
    initialStates <- generateAllCombinationBinary(ExampleNetwork$genes)
    trainingseries <- genereateBoolNetTimeseries(ExampleNetwork,
                                                 initialStates,
                                                 43,
                                                 type = "synchronous")

      
    cube <- constructFBNCube(target_genes = ExampleNetwork$genes,
                             conditional_genes = ExampleNetwork$genes,
                             timeseriesCube = trainingseries,
                             maxK = 5,
                             temporal = 1, 
                             useParallel = FALSE)
    NETWORK2 <- mineFBNNetwork(cube,
                               ExampleNetwork$genes)
    
    return(list(network = NETWORK2,
                initialStates = initialStates,
                timeseries = trainingseries))
  })
}

setupTestData <- function() {
    genesInput <- c("CycD", "p27", "CycE", "E2F")
    testseries <- list()
    testseries[[1]] <- matrix(c(1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1),
                              nrow = 4,
                              ncol = 4,
                              byrow = FALSE,
                              dimnames = list(genesInput, c("1", "2", "3", "4")))
    testseries[[2]] <- matrix(c(1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1),
                              nrow = 4,
                              ncol = 4, 
                              byrow = FALSE,
                              dimnames = list(genesInput, c("1", "2", "3", "4")))
    testseries[[3]] <- matrix(c(1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1), 
                              nrow = 4,
                              ncol = 4,
                              byrow = FALSE,
                              dimnames = list(genesInput, c("1", "2", "3", "4")))
    return(list(genes = genesInput,
                timeseries = testseries))
}

describe("run synchronous should succeed", {
    test_info <- setup(setupdata())
    it("reconstruct timeseries with parallel", {

        resultfile <- expect_error(reconstructTimeseries(test_info$network,
                                                         test_info$initialStates,
                                                         type = "synchronous",
                                                         maxTimepoints = 43,
                                                         useParallel = FALSE), NA)

        similarreport <- expect_error(generateSimilaryReport(test_info$timeseries,
                                                             resultfile), NA)
        expect_equal(similarreport$ErrorRate, 0)
        expect_equal(similarreport$AccurateRate, 1)
        expect_equal(similarreport$MissMatchedRate, 0)
        expect_equal(similarreport$PerfectMatchedRate, 1)
    })

    it("reconstruct timeseries", {
        resultfile <- expect_error(reconstructTimeseries(
          test_info$network,
          test_info$initialStates,
          type = "synchronous", 
          maxTimepoints = 43,
          useParallel = FALSE), NA)

        similarreport <- expect_error(generateSimilaryReport(
          test_info$timeseries,
          resultfile), NA)
        expect_equal(similarreport$ErrorRate, 0)
        expect_equal(similarreport$AccurateRate, 1)
        expect_equal(similarreport$MissMatchedRate, 0)
        expect_equal(similarreport$PerfectMatchedRate, 1)

    })

    it("second reconstruct timeseries", {


        reconstructed_timeseries <- expect_error(reconstructTimeseries(
          test_info$network, 
          test_info$initialStates,
          type = "synchronous", maxTimepoints = 43, useParallel = FALSE), NA)
        

        cube <- expect_error(constructFBNCube(target_genes = test_info$network$genes,
                                              conditional_genes = test_info$network$genes,
                                              timeseriesCube = reconstructed_timeseries,
                                              maxK = 5, 
                                              temporal = 1, 
                                              useParallel = FALSE),
                             NA)
        NETWORK2 <- expect_error(mineFBNNetwork(cube, 
                                                test_info$network$genes),
                                 NA)
        expect_error(FBNNetwork.Graph(NETWORK2), NA)
        print(test_info$network)
        print(NETWORK2)

        expect_true(identical(test_info$network, NETWORK2))

        resultfile <- expect_error(reconstructTimeseries(NETWORK2, 
                                                         test_info$initialStates,
                                                         type = "synchronous",
                                                         maxTimepoints = 43,
                                                         useParallel = FALSE), NA)

        similarreport <- expect_error(generateSimilaryReport(reconstructed_timeseries,
                                                             resultfile), NA)
        expect_equal(similarreport$ErrorRate, 0)
        expect_equal(similarreport$AccurateRate, 1)
        expect_equal(similarreport$MissMatchedRate, 0)
        expect_equal(similarreport$PerfectMatchedRate, 1)

    })
    
    it("second reconstruct timeseries with decay ==2", {

        test_info$network$timedecay[c(1, 2, 3, 4, 5)] <- 2
        reconstructed_timeseries <- expect_error(
          reconstructTimeseries(test_info$network,
                                test_info$initialStates,
                                type = "synchronous",
                                maxTimepoints = 43,
                                useParallel = FALSE),
            NA)
        cube <- expect_error(constructFBNCube(test_info$network$genes,
                                              test_info$network$genes,
                                              reconstructed_timeseries, 
                                              5, 1, FALSE), NA)
        NETWORK2 <- expect_error(mineFBNNetwork(cube, 
                                                test_info$network$genes), 
                                 NA)
        expect_error(FBNNetwork.Graph(NETWORK2), NA)
        print(test_info$network)
        print(NETWORK2)

        resultfile <- expect_error(reconstructTimeseries(NETWORK2,
                                                         test_info$initialStates,
                                                         type = "synchronous", 
                                                         maxTimepoints = 43,
                                                         useParallel = FALSE), NA)

        similarreport <- expect_error(generateSimilaryReport(reconstructed_timeseries,
                                                             resultfile), NA)
        expect_equal(similarreport$ErrorRate, 0)
        expect_equal(similarreport$AccurateRate, 1)
        expect_equal(similarreport$MissMatchedRate, 0)
        expect_equal(similarreport$PerfectMatchedRate, 1)

    })

    it("getProbabilityFromFunctionInput should succeed", {
        fns <- test_info$network$interactions[[test_info$network$genes[[1]]]]
        # find all activators' probabilities
        condOfActivation <- sapply(fns, function(activator) activator$type == 1L)
        funcOfActivators <- fns[condOfActivation]
        # find all inhibitors' probabilities
        condOfInhibitors <- sapply(fns, function(inhibitor) inhibitor$type == 0L)
        funcOfInhibitors <- fns[condOfInhibitors]
        pregeneInput <- dissolve(lapply(funcOfActivators[[1]]$input,
                                        function(geneindex) {
            res <- list()
            res[[1]] <- test_info$initialStates[[1]][[geneindex]]
            names(res)[[1]] <- test_info$network$genes[[geneindex]]
            return(res)
        }))

        expect_error(getProbabilityFromFunctionInput(1, 
                                                     funcOfActivators[[1]]$expression,
                                                     funcOfActivators[[1]]$probability,
                                                     pregeneInput), NA)
    })
})

describe("run synchronous with different timestep should succeed", {
    test_info <- setup(setupdata())
    it("reconstruct cube with timporal = 2", {
        reconstructed_timeseries <- expect_error(
          reconstructTimeseries(test_info$network, 
                                test_info$initialStates, 
                                type = "synchronous", 
                                maxTimepoints = 43, 
                                useParallel = FALSE),
            NA)
        cube <- expect_error(constructFBNCube(target_genes = test_info$network$genes,
                                              conditional_genes = test_info$network$genes,
                                              timeseriesCube = reconstructed_timeseries,
                                              maxK = 5, 
                                              temporal = 2, 
                                              useParallel = FALSE), NA)
        NETWORK2 <- expect_error(mineFBNNetwork(cube,
                                                test_info$network$genes), NA)
        expect_error(FBNNetwork.Graph(NETWORK2), NA)
        print(test_info$network)
        print(NETWORK2)

        # expect_true(identical(test_info$network,NETWORK2))

        resultfile <- expect_error(reconstructTimeseries(NETWORK2,
                                                         test_info$initialStates, 
                                                         type = "synchronous",
                                                         maxTimepoints = 43, 
                                                         useParallel = FALSE), NA)

        similarreport <- expect_error(generateSimilaryReport(reconstructed_timeseries,
                                                             resultfile), NA)
        expect_equal(similarreport$ErrorRate, 0)
        expect_equal(similarreport$AccurateRate, 1)
        expect_equal(similarreport$MissMatchedRate, 0)
        expect_equal(similarreport$PerfectMatchedRate, 1)
    })

    it("reconstruct cube with timporal = 2, decay=2", {
        test_info$network$timedecay[c(1, 2, 3, 4, 5)] <- 2
        reconstructed_timeseries <- expect_error(
          reconstructTimeseries(test_info$network,
                                test_info$initialStates,
                                type = "synchronous", 
                                maxTimepoints = 43, 
                                useParallel = FALSE),
            NA)
        cube <- expect_error(constructFBNCube(test_info$network$genes,
                                              test_info$network$genes,
                                              reconstructed_timeseries, 5, 2, FALSE),
                             NA)
        NETWORK2 <- expect_error(mineFBNNetwork(cube,
                                                test_info$network$genes), NA)
        expect_error(FBNNetwork.Graph(NETWORK2), NA)
        print(test_info$network)
        print(NETWORK2)

        # expect_true(identical(test_info$network,NETWORK2))

        resultfile <- expect_error(reconstructTimeseries(NETWORK2,
                                                         test_info$initialStates,
                                                         type = "synchronous",
                                                         maxTimepoints = 43, 
                                                         useParallel = FALSE),
                                   NA)

        similarreport <- expect_error(generateSimilaryReport(reconstructed_timeseries,
                                                             resultfile), NA)
        expect_equal(similarreport$ErrorRate, 0)
        expect_equal(similarreport$AccurateRate, 1)
        expect_equal(similarreport$MissMatchedRate, 0)
        expect_equal(similarreport$PerfectMatchedRate, 1)
    })

    # it("reconstruct timeseries with alerted networks, with max timeporal =2", {
    #     # test timporal =3
    #     test_info$network$interactions$Gene1$Gene1_1_Activator$timestep <- 2
    #     test_info$network$interactions$Gene5$Gene5_1_Activator$timestep <- 2
    #     reconstructed_timeseries <- expect_error(
    #       reconstructTimeseries(test_info$network,
    #                             test_info$initialStates,
    #                             type = "synchronous", 
    #                             maxTimepoints = 43, 
    #                             useParallel = FALSE),
    #         NA)
    #     cube <- expect_error(constructFBNCube(test_info$network$genes,
    #                                           test_info$network$genes,
    #                                           reconstructed_timeseries, 5, 2, FALSE)
    #                          , NA)
    #     NETWORK2 <- expect_error(mineFBNNetwork(cube, test_info$network$genes), NA)
    #     expect_error(FBNNetwork.Graph(NETWORK2), NA)
    #     print(test_info$network)
    #     print(NETWORK2)
    # 
    #     # expect_true(identical(test_info$network,NETWORK2))
    # 
    #     resultfile <- expect_error(reconstructTimeseries(NETWORK2,
    #                                                      test_info$initialStates,
    #                                                      type = "synchronous",
    #                                                      maxTimepoints = 43,
    #                                                      useParallel = FALSE),
    #                                NA)
    # 
    #     similarreport <- expect_error(generateSimilaryReport(reconstructed_timeseries,
    #                                                          resultfile), NA)
    #     expect_equal(similarreport$ErrorRate, 0)
    #     expect_equal(similarreport$AccurateRate, 1)
    #     expect_equal(similarreport$MissMatchedRate, 0)
    #     expect_equal(similarreport$PerfectMatchedRate, 1)
    # })
# 
#     it("reconstruct timeseries with alerted networks, with max timeporal =, decay =2", {
#         # test timporal =3
#         test_info$network$timedecay[c(1, 2, 3, 4, 5)] <- 2
# 
#         test_info$network$interactions$Gene1$Gene1_1_Activator$timestep <- 2
#         test_info$network$interactions$Gene5$Gene5_1_Activator$timestep <- 2

# 
#         reconstructed_timeseries <- expect_error(reconstructTimeseries(test_info$network, test_info$initialStates, type = "synchronous", maxTimepoints = 43, useParallel = TRUE),
#             NA)
#         cube <- expect_error(constructFBNCube(test_info$network$genes, test_info$network$genes, reconstructed_timeseries, 5, 2, TRUE), NA)
# 
#         NETWORK2 <- expect_error(mineFBNNetwork(cube, test_info$network$genes), NA)
#         expect_error(FBNNetwork.Graph(NETWORK2), NA)
#         print(test_info$network)
#         print(NETWORK2)
# 
#         # expect_true(identical(test_info$network,NETWORK2))
# 
#         resultfile <- expect_error(reconstructTimeseries(NETWORK2, test_info$initialStates, type = "synchronous", maxTimepoints = 43, useParallel = TRUE), NA)
# 
#         similarreport <- expect_error(generateSimilaryReport(reconstructed_timeseries, resultfile), NA)
#         expect_equal(similarreport$ErrorRate, 0)
#         expect_equal(similarreport$AccurateRate, 1)
#         expect_equal(similarreport$MissMatchedRate, 0)
#         expect_equal(similarreport$PerfectMatchedRate, 1)
#     })
# 
#     it("reconstruct cube with timporal = 3", {
#         reconstructed_timeseries <- expect_error(reconstructTimeseries(test_info$network, test_info$initialStates, type = "synchronous", maxTimepoints = 43, useParallel = TRUE),
#             NA)
#         cube <- expect_error(constructFBNCube(test_info$network$genes, test_info$network$genes, reconstructed_timeseries, 5, 3, TRUE), NA)
# 
#         NETWORK2 <- expect_error(mineFBNNetwork(cube, test_info$network$genes), NA)
#         expect_error(FBNNetwork.Graph(NETWORK2), NA)
#         print(test_info$network)
#         print(NETWORK2)
# 
#         # expect_true(identical(test_info$network,NETWORK2))
# 
#         resultfile <- expect_error(reconstructTimeseries(NETWORK2, test_info$initialStates, type = "synchronous", maxTimepoints = 43, useParallel = TRUE), NA)
# 
#         similarreport <- expect_error(generateSimilaryReport(reconstructed_timeseries, resultfile), NA)
#         expect_equal(similarreport$ErrorRate, 0)
#         expect_equal(similarreport$AccurateRate, 1)
#         expect_equal(similarreport$MissMatchedRate, 0)
#         expect_equal(similarreport$PerfectMatchedRate, 1)
#     })
#     it("reconstruct cube with timporal = 3, decay=3", {
#         test_info$network$timedecay[c(1, 2, 3, 4, 5)] <- 3
#         reconstructed_timeseries <- expect_error(reconstructTimeseries(test_info$network, test_info$initialStates, type = "synchronous", maxTimepoints = 43, useParallel = TRUE),
#             NA)
#         cube <- expect_error(constructFBNCube(test_info$network$genes, test_info$network$genes, reconstructed_timeseries, 5, 3, TRUE), NA)
#         NETWORK2 <- expect_error(mineFBNNetwork(cube, test_info$network$genes), NA)
#         expect_error(FBNNetwork.Graph(NETWORK2), NA)
#         print(test_info$network)
#         print(NETWORK2)
# 
#         # expect_true(identical(test_info$network,NETWORK2))
# 
#         resultfile <- expect_error(reconstructTimeseries(NETWORK2, test_info$initialStates, type = "synchronous", maxTimepoints = 43, useParallel = TRUE), NA)
# 
#         similarreport <- expect_error(generateSimilaryReport(reconstructed_timeseries, resultfile), NA)
#         expect_equal(similarreport$ErrorRate, 0)
#         expect_equal(similarreport$AccurateRate, 1)
#         expect_equal(similarreport$MissMatchedRate, 0)
#         expect_equal(similarreport$PerfectMatchedRate, 1)
#     })
# 
#     it("reconstruct timeseries with alerted networks to timestep =2, with max timeporal =3", {
#         # test timporal =3
#         test_info$network$interactions$Gene1$Gene1_1_Activator$timestep <- 2
#         test_info$network$interactions$Gene5$Gene5_1_Activator$timestep <- 2
# 
#         # test_info$network$timedecay[c(1,2,3,4,5)]<-3
# 
#         reconstructed_timeseries <- expect_error(reconstructTimeseries(test_info$network, test_info$initialStates, type = "synchronous", maxTimepoints = 43, useParallel = TRUE),
#             NA)
#         cube <- expect_error(constructFBNCube(test_info$network$genes, test_info$network$genes, reconstructed_timeseries, 5, 3, TRUE), NA)
# 
#         NETWORK2 <- expect_error(mineFBNNetwork(cube, test_info$network$genes), NA)
#         expect_error(FBNNetwork.Graph(NETWORK2), NA)
#         print(test_info$network)
#         print(NETWORK2)
# 
#         # expect_true(identical(test_info$network,NETWORK2))
# 
#         resultfile <- expect_error(reconstructTimeseries(NETWORK2, test_info$initialStates, type = "synchronous", maxTimepoints = 43, useParallel = TRUE), NA)
# 
#         similarreport <- expect_error(generateSimilaryReport(reconstructed_timeseries, resultfile), NA)
#         expect_equal(similarreport$ErrorRate, 0)
#         expect_equal(similarreport$AccurateRate, 1)
#         expect_equal(similarreport$MissMatchedRate, 0)
#         expect_equal(similarreport$PerfectMatchedRate, 1)
#     })
# 
#     it("reconstruct timeseries with alerted networks to timestep =3, with max timeporal =3", {
#         # test timporal =3
#         test_info$network$interactions$Gene1$Gene1_1_Activator$timestep <- 3
#         test_info$network$interactions$Gene5$Gene5_1_Activator$timestep <- 3

# 
#         # test_info$network$timedecay[c(1,2,3,4,5)]<-3
# 
#         reconstructed_timeseries <- expect_error(reconstructTimeseries(test_info$network, test_info$initialStates, type = "synchronous", maxTimepoints = 43, useParallel = TRUE),
#             NA)
#         cube <- expect_error(constructFBNCube(test_info$network$genes, test_info$network$genes, reconstructed_timeseries, 5, 3, TRUE), NA)
# 
#         NETWORK2 <- expect_error(mineFBNNetwork(cube, test_info$network$genes), NA)
#         expect_error(FBNNetwork.Graph(NETWORK2), NA)
#         print(test_info$network)
#         print(NETWORK2)
# 
#         # expect_true(identical(test_info$network,NETWORK2))
# 
#         resultfile <- expect_error(reconstructTimeseries(NETWORK2, test_info$initialStates, type = "synchronous", maxTimepoints = 43, useParallel = TRUE), NA)
# 
#         similarreport <- expect_error(generateSimilaryReport(reconstructed_timeseries, resultfile), NA)
#         expect_equal(similarreport$ErrorRate, 0)
#         expect_equal(similarreport$AccurateRate, 1)
#         expect_equal(similarreport$MissMatchedRate, 0)
#         expect_equal(similarreport$PerfectMatchedRate, 1)
#     })
# 
#     it("reconstruct timeseries with alerted networks to various timesteps, with max timeporal =4", {
#         # test timporal =4
#         test_info$network$interactions$Gene1$Gene1_1_Activator$timestep <- 3
#         test_info$network$interactions$Gene5$Gene5_1_Activator$timestep <- 4
#         test_info$network$interactions$Gene4$Gene4_1_Activator$timestep <- 4
#         # test_info$network$timedecay[c(1,2,3,4,5)]<-3
# 
#         reconstructed_timeseries <- expect_error(reconstructTimeseries(test_info$network, test_info$initialStates, type = "synchronous", maxTimepoints = 43, useParallel = TRUE),
#             NA)
#         cube <- expect_error(constructFBNCube(test_info$network$genes, test_info$network$genes, reconstructed_timeseries, 5, 4, TRUE), NA)
# 
#         NETWORK2 <- expect_error(mineFBNNetwork(cube, test_info$network$genes), NA)
#         expect_error(FBNNetwork.Graph(NETWORK2), NA)
#         print(test_info$network)
#         print(NETWORK2)
# 
#         # expect_true(identical(test_info$network,NETWORK2))
# 
#         resultfile <- expect_error(reconstructTimeseries(NETWORK2, test_info$initialStates, type = "synchronous", maxTimepoints = 43, useParallel = TRUE), NA)
# 
#         similarreport <- expect_error(generateSimilaryReport(reconstructed_timeseries, resultfile), NA)
#         expect_equal(similarreport$ErrorRate, 0)
#         expect_equal(similarreport$AccurateRate, 1)
#         expect_equal(similarreport$MissMatchedRate, 0)
#         expect_equal(similarreport$PerfectMatchedRate, 1)
#     })
# 
# 
#     it("reconstruct timeseries with alerted networks to timestep =4, with max timeporal =4", {
#         # test timporal =4
#         test_info$network$interactions$Gene1$Gene1_1_Activator$timestep <- 4
#         test_info$network$interactions$Gene5$Gene5_1_Activator$timestep <- 4
#         test_info$network$interactions$Gene4$Gene4_1_Activator$timestep <- 4
# 
#         # test_info$network$timedecay[c(1,2,3,4,5)]<-3
# 
#         reconstructed_timeseries <- expect_error(reconstructTimeseries(test_info$network, test_info$initialStates, type = "synchronous", maxTimepoints = 43, useParallel = TRUE),
#             NA)
#         cube <- expect_error(constructFBNCube(test_info$network$genes, test_info$network$genes, reconstructed_timeseries, 5, 4, TRUE), NA)
# 
#         NETWORK2 <- expect_error(mineFBNNetwork(cube, test_info$network$genes), NA)
#         expect_error(FBNNetwork.Graph(NETWORK2), NA)
#         print(test_info$network)
#         print(NETWORK2)
# 
#         # expect_true(identical(test_info$network,NETWORK2))
# 
#         resultfile <- expect_error(reconstructTimeseries(NETWORK2, test_info$initialStates, type = "synchronous", maxTimepoints = 43, useParallel = TRUE), NA)
# 
#         similarreport <- expect_error(generateSimilaryReport(reconstructed_timeseries, resultfile), NA)
#         expect_equal(similarreport$ErrorRate, 0)
#         expect_equal(similarreport$AccurateRate, 1)
#         expect_equal(similarreport$MissMatchedRate, 0)
#         expect_equal(similarreport$PerfectMatchedRate, 1)
#     })
# 
# })
# 
# describe("Run asynchronous timeseries should succeed", {
#     test_info <- setup(setupdata())
#     it("FBN asynchronous with max timeporal =5", {
# 
#         reconstructed_timeseries <- expect_error(reconstructTimeseries(test_info$network, test_info$initialStates, type = "asynchronous", maxTimepoints = 43, useParallel = TRUE),
#             NA)
#         cube <- expect_error(constructFBNCube(test_info$network$genes, test_info$network$genes, reconstructed_timeseries, 5, 5, TRUE), NA)
# 
#         NETWORK2 <- expect_error(mineFBNNetwork(cube, test_info$network$genes), NA)
#         expect_error(FBNNetwork.Graph(NETWORK2), NA)
#         print(test_info$network)
#         print(NETWORK2)
# 
#         # expect_true(identical(test_info$network,NETWORK2))
# 
#         resultfile <- expect_error(reconstructTimeseries(NETWORK2, test_info$initialStates, type = "synchronous", maxTimepoints = 43, useParallel = TRUE), NA)
# 
#         similarreport <- expect_error(generateSimilaryReport(reconstructed_timeseries, resultfile), NA)
#         expect_equal(similarreport$ErrorRate, 0)
#         expect_equal(similarreport$AccurateRate, 1)
#         expect_equal(similarreport$MissMatchedRate, 0)
#         expect_equal(similarreport$PerfectMatchedRate, 1)
#     })
# 
#     it("BoolNet asynchronou with max timeporal =5", {
#         require(BoolNet)
# 
#         network <- loadNetwork("D:\\Dropbox/Dropbox/FBNNet/tempdata/Example.txt")
#         print(network)
#         initialStates <- generateAllCombinationBinary(network$genes)
#         trainingseries <- genereateBoolNetTimeseries(network, initialStates, 43, type = "asynchronous")
#         # trainingseries<-FBNDataReduction(BoolNet::generateTimeSeries(network,1042,43,type='asynchronous'))
# 
#         # targetgenes,allgenes,timeseriesCube,maxK=5,temporal=1,useParallel=FALSE,optimizedTemoral=FALSE
#         cube <- constructFBNCube(network$genes, network$genes, trainingseries, 5, 5, TRUE)
#         NETWORK2 <- mineFBNNetwork(cube, network$genes)
#         FBNNetwork.Graph(NETWORK2)
#         print(NETWORK2)
# 
#         resultfile <- reconstructTimeseries(NETWORK2, initialStates, type = "synchronous", maxTimepoints = 43, useParallel = TRUE)
#         similarreport <- expect_error(generateSimilaryReport(trainingseries, resultfile), NA)
#         # expect_equal(similarreport$ErrorRate,0) expect_equal(similarreport$AccurateRate,1) expect_equal(similarreport$MissMatchedRate,0)
#         # expect_equal(similarreport$PerfectMatchedRate,1)
# 
#         print("Benchmark result->FBNNet")
#         print(paste("ErrorRate=", similarreport$ErrorRate, sep = "", collapse = ""))
#         print(paste("AccurateRate=", similarreport$AccurateRate, sep = "", collapse = ""))
#         print(paste("MissMatchedRate=", similarreport$MissMatchedRate, sep = "", collapse = ""))
#         print(paste("PerfectMatchedRate=", similarreport$PerfectMatchedRate, sep = "", collapse = ""))
# 
#         # print('++++++++using boolean reconstructnetwork on short time series+++++++++++++++++++++++++++++')
#         # boolnet_net<-reconstructNetwork(trainingseries,readableFunctions = TRUE) boolnet_net<-BoolNet::chooseNetwork(boolnet_net,rep(1, length(boolnet_net$genes)))
#         # #converted to BooleanNetwork.  print(boolnet_net) trainingseries_Boolnet<-genereateBoolNetTimeseries(boolnet_net,initialStates,43, type='synchronous')
#         # report<-generateSimilaryReport(trainingseries_Boolnet,trainingseries) print('Benchmark result->BoolNet')
#         # print(paste('ErrorRate=',report$ErrorRate,sep='',collapse = '')) print(paste('AccurateRate=',report$AccurateRate,sep='',collapse = ''))
#         # print(paste('MissMatchedRate=',report$MissMatchedRate,sep='',collapse = '')) print(paste('PerfectMatchedRate=',report$PerfectMatchedRate,sep='',collapse =
#         # ''))
#     })


})

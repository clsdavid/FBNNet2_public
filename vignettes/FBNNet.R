## ---- out.width = "100%"------------------------------------------------------
knitr::include_graphics("images/evaluation.jpg")

## ---- out.width = "100%"------------------------------------------------------
library(BoolNet)
library(utils)
library(FBNNet)
data("ExampleNetwork")
ExampleNetwork
print(ExampleNetwork)
## show BoolNet type of network graph
BoolNet::plotNetworkWiring(ExampleNetwork)

## ---- out.width = "100%"------------------------------------------------------
   initialStates <- generateAllCombinationBinary(ExampleNetwork$genes)
   trainingseries <- genereateBoolNetTimeseries(ExampleNetwork,
                                             initialStates,43, 
                                             type = "synchronous")
   FBNcellcyclenetwork <- generateFBMNetwork(timeseries_data = trainingseries, 
                                 maxK = 4, 
                                 max_deep_temporal = 1, 
                                 useParallel = FALSE,
                                 verbose = FALSE)
   print(FBNcellcyclenetwork)

   resultfile <- reconstructTimeseries(FBNcellcyclenetwork,
                                    initialStates,
                                    type = "synchronous", 
                                    maxTimepoints = 43, 
                                    useParallel = FALSE)

   similarreport <- generateSimilaryReport(trainingseries,resultfile)
   print(paste("ErrorRate=",similarreport$ErrorRate,sep = "",collapse = ""))
   print(paste("AccurateRate=",similarreport$AccurateRate,sep = "",collapse = ""))
   print(paste("MissMatchedRate=",similarreport$MissMatchedRate,sep = "",collapse = ""))
   print(paste("PerfectMatchedRate=",similarreport$PerfectMatchedRate,sep = "",collapse = ""))

   #get attractors
   genes <- rownames(trainingseries[[1]])

   attractor <- searchForAttractors(FBNcellcyclenetwork,initialStates,genes)
   print(attractor)
   #display the dynamic trajectory of the attactor 2
   FBNNetwork.Graph.DrawAttractor(FBNcellcyclenetwork,attractor,2)


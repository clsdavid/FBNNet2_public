---
title: "Fundamental Boolean Network"
pagetitle: "Fundamental Boolean Network"
author: "Leshi Chen"
affiliation: Lincoln University, New Zealand
date: "2020-08-05"
url: https://www.frontiersin.org/articles/10.3389/fphys.2018.01328/full
output:
  rmarkdown::html_document:
    citation_package: natbib
    fig_caption: true
    theme: null
    lib_dir: libs
#output: rmarkdown::pdf_document
bibliography: master.bib
header-includes:
  -  \usepackage{hyperref}
always_allow_html: true
keywords: "Fundamental Boolean Modelling, FBNNet, FBN, FBM"
geometry: margin=1in
fontsize: 11pt
spacing: double
endnote: no
graphics: yes
abstract: "A Boolean model is a simple, discrete and dynamic model without the need to consider the effects at the intermediate levels and is powerful in qualitatively describing large-scale system dynamics. However, little effort has been made into constructing activation and inhibition networks, which could indicate the direct roles of a gene (or its synthesised protein) as an activator or inhibitor of a target gene. The major reason for this is that the hypotheses of the current Boolean models do not provide an intuitive way to identify the effects of individual activation or inhibition pathways on the target gene. Therefore, we propose to focus on the general Boolean functions at the subfunction level and further split the subfunctions into the activation and inhibition domains. As a consequence, we developed a novel data-driven Boolean model; namely, the Fundamental Boolean Model (FBM), to draw insights into gene activation and inhibition. This novel Boolean model provided an intuitive definition of activation and inhibition pathways and included mechanisms to handle protein decay issues as well as introducing uncertainty into the Boolean functions. Apart from this new concept, we also proposed a new data mining method to extract fundamental Boolean networks (FBNs) from time series data. To prove the concept of the novel model, we implemented a platform using R language, called FBNNet, based on the proposed novel Boolean model and a novel data mining technique. Our experimental results showed that the proposed FBM can explicitly display the internal connections between genes separated into the connection types of activation and inhibition, and the novel Boolean model was shown to infer GRNs from the time series data generated from an original cell cycle network with a high degree of accuracy. Moreover, the method we proposed to infer the gene regulatory networks for the novel Boolean model can be run in parallel and; hence, the computation cost is affordable. To summarise, the novel Boolean model and related FBNs could show significant trajectories in genes to reveal how genes regulated each other over a given time period. This new feature could facilitate further research about drug interactions to detect the side effects of the use of a newly-proposed drug. (source:https://www.frontiersin.org/articles/10.3389/fphys.2018.01328/full)"
vignette: >
  %\VignetteIndexEntry{P}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


## Introduction
Fundamental Boolean modelling(FBM) has been proposed to draw insights into gene activation, inhibition, and protein decay. This novel Boolean model provides an intuitive definition of activation and inhibition pathways and includes mechanisms to handle protein decay issues. To prove the concept of the novel model, we implemented a platform using R language, called FBNNet. Our experimental results show that the proposed FBM could explicitly display the internal
connections of the mammalian cell cycle between genes separated into the connection
types of activation, inhibition and protein decay [@fbnnet-2018].

## Experiment
The experiments conducted and described here intend to prove the concept of the new Boolean Model, i.e., the FBM. To verify the results, we apply the general processes described in Figure 1 as a benchmark to compare the results generated via BoolNet [@boolnet-2010], with these consequencely reconstructed from the new R package, FBNNet.


```r
knitr::include_graphics("IMAGES/evaluation.jpg")
```

<img src="IMAGES/evaluation.jpg" title="plot of chunk unnamed-chunk-1" alt="plot of chunk unnamed-chunk-1" width="100%" />
```
Figure 1 Evaluation assessment of Fundamental Boolean network inference using BoolNet. The blue arrows represent the processes using BoolNet and brown arrows represent the processes using our FBNNet R package. The green arrows represent the evaluation process. (A) We use the BoolNet script loadNetwork.R to load pre-defined networks from files and then generate the time series and networks; (B) We use the time series generated from BoolNet and the new R package, FBNNet, to generate FBNs; (C) We reconstruct the time series via the FBM; this process can be used to expand the short time series data; (D) To evaluate the FBM, we rebuild the BoolNet type network based on the reconstructed time series; and (E) We evaluate the FBN inference methods by comparing the generated time series and the generated BoolNet type of network with the original time series and network that were generated in step A.
```
## Data

```r
library(knitr)
library(BoolNet)
library(utils)
library(FBNNet)
library(visNetwork)
data("ExampleNetwork")
ExampleNetwork
#> Boolean network with 5 genes
#> 
#> Involved genes:
#> Gene1 Gene2 Gene3 Gene4 Gene5
#> 
#> Transition functions:
#> Gene1 = Gene1
#> Gene2 = Gene1 & Gene5 & !Gene4
#> Gene3 = Gene3
#> Gene4 = Gene3 & !(Gene1 & Gene5)
#> Gene5 = !Gene2
```

## Extract the Fundamental Boolean Network

```r
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
#> Fundamental Boolean Network with  5 genes
#> Genes involved:
#> Gene1, Gene2, Gene3, Gene4, Gene5
#> 
#> Networks:
#> Multiple Transition Functions for Gene1 with decay value = 1:
#> Gene1_1_Activator: Gene1 = Gene1 (Confidence: 1, TimeStep: 1)
#> Gene1_1_Inhibitor: Gene1 = !Gene1 (Confidence: 1, TimeStep: 1)
#> 
#> Multiple Transition Functions for Gene2 with decay value = 1:
#> Gene2_1_Activator: Gene2 = Gene1&!Gene4&Gene5 (Confidence: 1, TimeStep: 1)
#> Gene2_1_Inhibitor: Gene2 = !Gene1 (Confidence: 1, TimeStep: 1)
#> Gene2_2_Inhibitor: Gene2 = Gene4 (Confidence: 1, TimeStep: 1)
#> Gene2_3_Inhibitor: Gene2 = !Gene5 (Confidence: 1, TimeStep: 1)
#> 
#> Multiple Transition Functions for Gene3 with decay value = 1:
#> Gene3_1_Activator: Gene3 = Gene3 (Confidence: 1, TimeStep: 1)
#> Gene3_1_Inhibitor: Gene3 = !Gene3 (Confidence: 1, TimeStep: 1)
#> 
#> Multiple Transition Functions for Gene4 with decay value = 1:
#> Gene4_1_Activator: Gene4 = !Gene1&Gene3 (Confidence: 1, TimeStep: 1)
#> Gene4_2_Activator: Gene4 = Gene3&!Gene5 (Confidence: 1, TimeStep: 1)
#> Gene4_1_Inhibitor: Gene4 = !Gene3 (Confidence: 1, TimeStep: 1)
#> Gene4_2_Inhibitor: Gene4 = Gene1&Gene5 (Confidence: 1, TimeStep: 1)
#> 
#> Multiple Transition Functions for Gene5 with decay value = 1:
#> Gene5_1_Activator: Gene5 = !Gene2 (Confidence: 1, TimeStep: 1)
#> Gene5_1_Inhibitor: Gene5 = Gene2 (Confidence: 1, TimeStep: 1)
```



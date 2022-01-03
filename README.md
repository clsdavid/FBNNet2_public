# FBNNet
### Background
Fundamental Boolean Model (FBM), published in Chen et al. (2018) <https://doi.org/10.3389/fphys.2018.01328>, provides an intuitive definition of activation and inhibition pathways and includes mechanisms to handle protein decay issues. To prove the concept of the novel model, we implemented an R package, called FBNNet. Our experimental results show that the proposed FBM could explicitly display the internal connections of the mammalian cell cycle between genes separated into the connection types of activation, inhibition and protein decay. Moreover, the method we proposed to infer the gene regulatory networks for the novel Boolean model can be run in parallel and; hence, the computation cost is affordable. Finally, the novel Boolean model and related Fundamental Boolean Networks (FBNs) could show significant trajectories in genes to reveal how genes regulated each other over a given period. This new feature could facilitate further research on drug interventions to detect the side effects of a newly-proposed drug.

### Introduction
This package adopted the concepts of fundamental Boolean modelling and networks to provide mechanisms for extracting the fundamental Boolean networks from the microarray timeseires data. The methodologies implemented in this package are documented in Chen et al. (2018) <https://doi.org/10.3389/fphys.2018.01328>.

### Main functions
* `generateFBMNetwork`: is the main entry of the package FBNNet that can be used to mine the gene regulatory network.
* `plotNetwork`: plot the Fundamental boolean networks

### Data

#### Main Example Data
* `BoolNet_CellCycle_Network`: Data generated via BoolNet::BoolNet::loadNetwork
* `ExampleNetwork`: Data generated for testing purpose using the traditional Boolean network via BoolNet::BoolNet::loadNetwork. This Example described in the paper 'A Novel Data-Driven Boolean Model for Genetic Regulatory Networks' (Chen, L., etc., 2018)
* `ExampleTimeseriesData`: Data generated via this package and is in the concept of FBM. This Example described in the paper 'A Novel Data-Driven Boolean Model for Genetic Regulatory Networks' (Chen, L., etc., 2018)
* `FBNExampleNetworks`: Data generated via this package and is in the concept of FBM. This Example described in the paper 'A Novel Data-Driven Boolean Model for Genetic Regulatory Networks' (Chen, L., etc., 2018)
* `FBNcellcycleNetwork`: Data generated via this package and is in the concept of FBM. This Cell Cycle Network of FBN described in the paper 'A Novel Data-Driven Boolean Model for Genetic Regulatory Networks' (Chen, L., etc., 2018)


#### Experiment Data
* `Common_Genes_Leukeamia`: Data contains the 286 common genes generated via this package for the study of Fundamental Boolean Network Modelling for Childhood Acute Lymphoblastic Leukaemia Pathways.
* `Leukeamia_Timeseries`: Data contains the timeseries data for the common genes generated via this package for the study of Fundamental Boolean Network Modelling for Childhood Acute Lymphoblastic Leukaemia Pathways.
* `Leukeamia_Networks`: The fundamental Boolean networks mined from the Leukeamia_Timeseries data via this package for the study of Fundamental Boolean Network Modelling for Childhood Acute Lymphoblastic Leukaemia Pathways.
* `DAVID_Gene_List`: Data extracted from DAVID tools for mapping probeset ids with gene names

### Installation
#### Install development version `FBNNet` from GitHub:
devtools::install_github("clsdavid/FBNNet2_public")

#### Documentation
https://clsdavid.github.io/FBNNet2_public/

### SAMPLE CODE
#### TFBM
##### find forward related genes with FAA
TFBM_FAA_CDC42EP3_Networks <- findForwardRelatedNetworkByGenes(networks = TFBM_Leukeamia_Networks, target_gene_list = "CDC42EP3", regulationType = 1, target_type = 1, maxDeep = 1)
FBNNetwork.Graph(TFBM_FAA_CDC42EP3_Networks)

##### find forward related genes with FAI
TFBM_FIA_CDC42EP3_Networks <- findForwardRelatedNetworkByGenes(networks = TFBM_Leukeamia_Networks, target_gene_list = "CDC42EP3", regulationType = 0, target_type = 1, maxDeep = 1)
FBNNetwork.Graph(TFBM_FAA_CDC42EP3_Networks)

##### find forward related genes with FAI with 2 levels
TFBM_FAI_CDC42EP3_Networks_2 <- findForwardRelatedNetworkByGenes(networks = TFBM_Leukeamia_Networks, target_gene_list = "CDC42EP3", regulationType = 1, target_type = 0, maxDeep = 2, next_level_mix_type = TRUE)
FBNNetwork.Graph(TFBM_FAI_CDC42EP3_Networks_2)

##### find backward related genes with BA
TFBM_FAA_CDC42EP3_Networks <- findAllBackwardRelatedGenes(networks = TFBM_Leukeamia_Networks, target_gene = "CDC42EP3", regulationType = 0, target_type = 1, maxDeep = 1)
FBNNetwork.Graph(TFBM_FAA_CDC42EP3_Networks)

#### FBM
##### find forward related genes with FAA
FBM_FAA_CDC42EP3_Networks <- findForwardRelatedNetworkByGenes(networks = FBM_Leukeamia_Networks, target_gene_list = "CDC42EP3", regulationType = 1, target_type = 1, maxDeep = 1)
FBNNetwork.Graph(FBM_FAA_CDC42EP3_Networks)

##### find forward related genes with FAI
FBM_FAA_CDC42EP3_Networks <- findForwardRelatedNetworkByGenes(networks = FBM_Leukeamia_Networks, target_gene_list = "CDC42EP3", regulationType = 0, target_type = 1, maxDeep = 1)
FBNNetwork.Graph(FBM_FAA_CDC42EP3_Networks)

##### find forward related genes with FAI with 2 levels
FBM_FAI_CDC42EP3_Networks_2 <- findForwardRelatedNetworkByGenes(networks = FBM_Leukeamia_Networks, target_gene_list = "CDC42EP3", regulationType = 0, target_type = 1, maxDeep = 2, next_level_mix_type = TRUE)
FBNNetwork.Graph(FBM_FAI_CDC42EP3_Networks_2)

##### find backward related genes with BA
FBM_FAA_CDC42EP3_Networks <- findAllBackwardRelatedGenes(networks = FBM_Leukeamia_Networks, target_gene = "CDC42EP3", regulationType = 0, target_type = 1, maxDeep = 1)
FBNNetwork.Graph(FBM_FAA_CDC42EP3_Networks)


---
__Citation__
Chen, L., D. Kulasiri and S. Samarasinghe (2018). A Novel Data-Driven Boolean Model for Genetic Regulatory Networks. Front Physiol 9: 1328.

__Copyright and Licensing__
This package and related novel concepts were originally proposed and developed by Leshi Chen <https://doi.org/10.3389/fphys.2018.01328>, under the supervision of Don Kulasiri and Sandhya Samarasinghe during the PHD study at Lincoln University, New Zealand. Currently, the package is licensed under the MIT License for public.



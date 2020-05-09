# FBNNet
### Background
Fundamental Boolean Model (FBM), published in Chen et al. (2018) <https://doi.org/10.3389/fphys.2018.01328>, provides an intuitive definition of activation and inhibition pathways and includes mechanisms to handle protein decay issues. To prove the concept of the novel model, we implemented an R package, called FBNNet. Our experimental results show that the proposed FBM could explicitly display the internal connections of the mammalian cell cycle between genes separated into the connection types of activation, inhibition and protein decay. Moreover, the method we proposed to infer the gene regulatory networks for the novel Boolean model can be run in parallel and; hence, the computation cost is affordable. Finally, the novel Boolean model and related Fundamental Boolean Networks (FBNs) could show significant trajectories in genes to reveal how genes regulated each other over a given period. This new feature could facilitate further research on drug interventions to detect the side effects of a newly-proposed drug.

### Introduction
This package adopted the concepts of fundamental Boolean modelling and networks to provide mechanisms for extracting the fundamental Boolean networks from the microarray timeseires data. The methodologies implemented in this package are documented in Chen et al. (2018) <https://doi.org/10.3389/fphys.2018.01328>.

### Main functions
* `generateFBMNetwork`: is the main entry of the package FBNNet that can be used to mine the gene regulatory network.
* `plotNetwork`: plot the Fundamental boolean networks

### Main Example Data
* `BoolNet_CellCycle_Network`: Data generated via BoolNet::BoolNet::loadNetwork
* `ExampleNetwork`: Data generated for testing purpose using the traditional Boolean network via BoolNet::BoolNet::loadNetwork. This Example described in the paper 'A Novel Data-Driven Boolean Model for Genetic Regulatory Networks' (Chen, L., etc., 2018)
* `FBNSampleData`: Data generated via this package and is in the concept of FBM. This Example described in the paper 'A Novel Data-Driven Boolean Model for Genetic Regulatory Networks' (Chen, L., etc., 2018)

### Installation
#### Install development version `FBNNet` from GitHub:
devtools::install_github("clsdavid/FBNNet2_public")

#### Documentation
coming soon


---

__Copyright and Licensing__
This package and related novel concepts were originally proposed and implemented by Leshi Chen <https://doi.org/10.3389/fphys.2018.01328>, under the supervision of Don Kulasiri and Sandhya Samarasinghe during the PHD study at Lincoln University, New Zealand. Currently the package is licensed under the MIT License for public.



# MOSClip
### Multi Omics Survival Clip

_MOSClip_ is a new method that combines survival analysis and graphical model theory to test the survival association of pathways or of their connected components in a multi-omic framework.
Multi-omic gene measurements are tested as covariates of a Cox proportional hazard model after dimensionality reduction of data.
The final goal is to find the biological processes impacting the patient's survival. 

MOSClip has a modular structure, allowing the use of one or multiple different omics, as well as different data reduction strategies and tests.

Furthermore, in MOSClip multiple efforts have been dedicated to the implementation of specific graphical tools to browse,
manage and provide help in the interpretation of results.

To install MOSClip please run 

devtools::install_github("cavei/MOSClip")

from your R prompt.


Please follow this [tutorial](https://cavei.github.io/). The tutorial reproduce the analysis performed by MOSClip in the TCGA Ovarian Cancer dataset. You can download the initial dataset [here](https://github.com/cavei/example-datasets).


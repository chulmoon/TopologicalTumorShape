# Persistent Homology Topological Features for Medical Images

Tumor shape is a key factor that affects tumor growth and metastasis. In [Moon et al.](https://arxiv.org/abs/2012.12102), 
we propose a topological feature computed by persistent homology to characterize tumor progression from digital pathology and radiology images and examines its effect on the time-to-event data.
Two case studies are conducted using consecutive from lung cancer pathology imaging data from the National Lung Screening Trial (NLST) and brain tumor radiology imaging data from the The Cancer Imaging Archive (TCIA). The results of both studies show that the topological features predict survival prognosis after adjusting clinical variables, and the predicted high-risk groups have significantly (at the level of 0.001) worse survival outcomes than the low-risk groups.
Also, the topological shape features found to be positively associated with survival hazards are irregular and heterogeneous shape patterns, which are known to be related to tumor progression. 
    
# Data
1. Image files:
* For the NLST datasets, the raw imaging data are available at https://biometry.nci.nih.gov/cdas/datasets/nlst/. The users need to fill out data request form and apply for permission in order to download the data. We only provided three example datasets.
* For the TCIA datasets, the imaging data are available at https://github.com/lorinanthony/SECT. The users need to download image data from the public repository.

2. Distance transform:
The size of distance transform data is too large (about 1.8 GB in total), so we only include three examples for the NLST and the TCIA datasets, respectively. However, the same distance transform images could be computed given binary or trenary images. Due to the size issue, the simulated pixel-rearranged images are not provided.

3. Persistence diagrams:
Persistence diagram data files are provided in this repository. Persistence diagrams computed by [GUDHI](http://gudhi.gforge.inria.fr/) are given as txt files and given ./data/lung/PersistenceDiagram/ for the NLST data and ./data/brain/PersistenceDiagram/ for the TCIA data in the repository. The persistence diagram txt files have three columns of dimension, birth, and death.

4. Clinical information:
Clinical information data files of the NLST and TCIA datasets are provided in `clinical_info_lung.Rdata` and `clinical_info_brain.Rdata`, respectively. 

5. Selected smoothing parameter:
The smoothing parameters selected by leave-one-out cross-validations are given in `brain_loocv_result.Rdata` and `lung_loocv_result.Rdata`.

# Workflow
We provide the necessary code to distance transform the data (using Python), compute persistent homology with a cubical complex (using Python), represent persistence diagrams in a functional space (using R), fit the functional Cox proportional-hazards model and the Gaussian Process model of [Crawford et al.](https://doi.org/10.1080/01621459.2019.1671198) (using R), and analyze the results (using R).

## Model Estimation and Validation

1. Distance transform and compute persistent homology – using files in “./TopologicalTumorShape/code/Python/”:
* `PersistentHomologyBrain.py`: SEDT-2 transform and compute persistent homology for the TCIA brain tumor radiology images
* `PersistentHomologyLung.py`: SEDT-3 transform and compute persistent homology for the NLST lung cancer pathology images

2. Represent computed persistence diagrams as persistence functions and fit a survival model – using files in “./TopologicalTumorShape/code/R/”: 
* `model.brain.R`: Fit the Cox proportional-hazards and the functional Cox proportional-hazards models for the TCIA brain tumor images
* `model.brain.SECT.R`: Fit the Gaussian Process model of Crawford et al. (2020) for the TCIA brain tumor images
* `model.lung.R`: Fit the Cox proportional-hazards and the functional Cox proportional-hazards models for the NLST lung cancer pathology images

## Simulation Studies with Rearranged Pixels

1. Distance transform and compute persistent homology – using files in “./TopologicalTumorShape/code/Python/”:
* `ComputeDistance_brain_rearranged.py`: Rearrange pixels and SEDT-2 transform
* `ComputeDistance_lung_rearranged.py`: Rearrange pixels and SEDT-3 transform
* `ComputePD_brain_rearranged.py`: Compute persistent homology from rearranged brain tumor radiology images
* `ComputePD_lung_rearranged.py`: Compute persistent homology from rearranged lung cancer pathology images

2. Represent computed persistence diagrams as persistence functions and fit a survival model – using files in “./TopologicalTumorShape/code/R/”: 
* `model.brain.rearranged.R`: Fit the Cox proportional-hazards and the functional Cox proportional-hazards models for the pixel-rearranged TCIA brain tumor images
* `model.lung.rearranged.R`: Fit the Cox proportional-hazards and the functional Cox proportional-hazards models for the pixel-rearranged NLST lung cancer pathology images

3. Plot results – using files in “./TopologicalTumorShape/code/R/”: 
* `rearranged.result.R`: Draw summary plots using data files `prescv_brain.Rdata` and `prescv_lung.Rdata`.

# Relevant Citations


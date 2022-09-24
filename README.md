# Persistent Homology Topological Features for Medical Images

Tumor shape is a key factor that affects tumor growth and metastasis. In [Moon et al.](https://arxiv.org/abs/2012.12102), 
we propose a topological feature computed by persistent homology to characterize tumor progression from digital pathology and radiology images and examines its effect on the time-to-event data.
Two case studies are conducted using consecutive from lung cancer pathology imaging data from the National Lung Screening Trial (NLST) and brain tumor radiology imaging data from the The Cancer Imaging Archive (TCIA). The results of both studies show that the topological features predict survival prognosis after adjusting clinical variables, and the predicted high-risk groups have significantly (at the level of 0.001) worse survival outcomes than the low-risk groups.
Also, the topological shape features found to be positively associated with survival hazards are irregular and heterogeneous shape patterns, which are known to be related to tumor progression. Persistence homology over the cubical complex is computed by [GUDHI](http://gudhi.gforge.inria.fr/). The persistence diagram txt files have three columns of dimension, birth, and death.
    
# Simulation – in “./simulation”:
* `simulation_1_data_generation.ipynb`: Generate binary tumor images, SEDT-2 transform, and compute persistent homology for scenario 1
* `simulation_1_result.R`: Fit the Cox and the functional Cox proportional-hazards models and draw summary plots for scenario 1
* `simulation_2_data_generation.ipynb`: Generate binary tumor images, SEDT-2 transform, and compute persistent homology for scenario 2 
* `simulation_2_result.R`: Fit the Cox and the functional Cox proportional-hazards models and draw summary plots for scenario 2

# Lung cancer application – in “./lung”:
* For the NLST datasets, the raw imaging data are available at https://biometry.nci.nih.gov/cdas/datasets/nlst/. The users need to fill out data request form and apply for permission in order to download the data. We only provided three example datasets.
* The size of distance transform data is too large, so we only include three examples for the NLST dataset. However, the same distance transform images could be computed given images.
* `lung_sedt3_persistent_homology.ipynb`: Generate binary tumor images, SEDT-3 transform, and compute persistent homology for the NLST lung cancer pathology images
* `lung_fcoxph_functions.R`: Functions of the Cox and the functional Cox proportional-hazards models for the NLST lung cancer pathology images
* `lung_fcoxph_main.R`: Fit the Cox and the functional Cox proportional-hazards models and draw summary plots for the NLST lung cancer pathology images
* `lung_cv_mdw_function.R`: Functions for finding smoothing paramters for the maximum distance weight using cross validations
* `lung_cv_linear_function.R`: Functions for finding smoothing paramters for the linear weight using cross validations
* `lung_cv_pwgk_function.R`: Functions for finding smoothing paramters for the persistence weighted Gaussian kernel weight using cross validations
* `lung_cv_main.R`: Find smoothing parameters under the three weights
* `clinical_info_lung.Rdata`: clinical data for for the NLST lung cancer patients

# Lung simulation – in “./lung_simulation”:
* Due to the size issue, the simulated pixel-rearranged images are not provided.
* `lung_simulation_data_generation.ipynb`: Generate the lung cancer pathology images with false shape information, SEDT-3 transform, and compute persistent homology
* `lung_simulation_functions.R`: Functions of the Cox and the functional Cox proportional-hazards models for the the pixel-rearranged lung cancer pathology images
* `lung_simulation_main.R`: Fit the Cox and the functional Cox proportional-hazards models and draw summary plots for the the pixel-rearranged lung cancer pathology images
* `clinical_info_lung.Rdata`: clinical data for for the NLST lung cancer patients

# Brain tumor application - in "./brain"
* For the TCIA datasets, the imaging data are available at [the public repository](https://github.com/lorinanthony/SECT). The users need to download image data from the repository.
* The size of distance transform data is too large, so we only include three examples for the TCIA dataset. However, the same distance transform images could be computed given images.
* `brain_sedt2_persistent_homology.ipynb`: Generate binary tumor images, SEDT-3 transform, and compute persistent homology for the TCIA brain tumor images
* `brain_fcoxph_functions.R`: Functions of the Cox and the functional Cox proportional-hazards models for the TCIA brain tumor images
* `brain_fcoxph_main.R`: Fit the Cox and the functional Cox proportional-hazards models and draw summary plots for the TCIA brain tumor images
* `brain_fcoxph_sect.R`: Fit the functional Cox proportional-hazards model and draw summary plots for the TCIA brain tumor images using the Smooth Euler Characteristic Transform (SECT) [Crawford et al.](https://doi.org/10.1080/01621459.2019.1671198)
* `brain_cv_mdw_function.R`: Functions for finding smoothing paramters for the maximum distance weight using cross validations
* `brain_cv_linear_function.R`: Functions for finding smoothing paramters for the linear weight using cross validations
* `brain_cv_pwgk_function.R`: Functions for finding smoothing paramters for the persistence weighted Gaussian kernel weight using cross validations
* `brain_cv_main.R`: Find smoothing parameters under the three weights
* `clinical_data_brain.csv`: clinical data for for the TCIA brain tumor patients


# Relevant Citations


library(tidyverse)

# p-value obtained by LOOCV under various smoothing parameters

## lung cancer 
load("./data/lung/prescv_lung.Rdata")
prescv_lung=prescv
sig0cv = sig1cv = seq(0.1,3,by=0.1)
sigcv = expand.grid(sig0cv,sig1cv)
sigcv[prescv_lung==min(prescv_lung),]
table_lung = sigcv[prescv_lung==min(prescv_lung),] %>%
  arrange(Var1)
table_lung # list of smoothing parameter combinations

## brain tumor
load("./data/brain/prescv_brain.Rdata")
prescv_brain=prescv
sig0cv = sig1cv = seq(0.1,3,by=0.1)
sigcv = expand.grid(sig0cv,sig1cv)
sigcv[prescv_brain==min(prescv_brain),]
table_brain = sigcv[prescv_brain==min(prescv_brain),] %>%
  arrange(Var1)
table_brain # list of smoothing parameter combinations

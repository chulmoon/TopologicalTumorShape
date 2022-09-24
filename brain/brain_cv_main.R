library(snow)
library(parallel)
library(tidyverse)

getCurrentFileLocation <-  function()
{
  this_file <- commandArgs() %>% 
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file)==0)
  {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}

source(paste0(getCurrentFileLocation(),"/brain_fcoxph_functions.R"))

# slide ids
load(paste0(getCurrentFileLocation(),"/slideid100.Rdata"))
# it can be obtained from the code below
# set the directory path pulled from https://github.com/lorinanthony/SECT
#SECTdirectory=".../SECT"
#idres=getid(SECTdirectory)
#slideid100=idres$slideid100
#ids=idres$ids

pdrangeres=pdrange(slideid100)
filename=pdrangeres$filename
minmax0=pdrangeres$minmax0
minmax1=pdrangeres$minmax1
patientid=pdrangeres$patientid
ntumor=pdrangeres$ntumor
pd0.total=pdrangeres$pd0.total
pd1.total=pdrangeres$pd1.total

# clinical information
clinical_data_brain = read_csv(paste0(getCurrentFileLocation(),"/clinical_data_brain.csv"))

# smoothing parameters
sig0cv = sig1cv = seq(0.2,5,by=0.2)
sigcv = expand.grid(sig0cv,sig1cv)

# mdw
## CV
ncore = parallel::detectCores()
cl = snow::makeCluster(ncore-1, type="SOCK")
clusterEvalQ(cl,{library(spatstat);library(tidyverse);library(tools);library(survival)})
source(paste0(getCurrentFileLocation(),"/brain_cv_mdw_function.R"))
prescv=snow::parApply(cl,matrix(1:nrow(sigcv),nrow(sigcv),1),1,main,
                      clinical_data_brain,minmax0,minmax1,
                      slideid100,ids,patientid,ntumor,pd0.total,pd1.total) # pvalue results from cv
stopCluster(cl)
save(prescv, file = paste0(getCurrentFileLocation(),
                           "/brain_cv_mdw_result.Rdata", sep=""))
## plot
sigcv.df = data.frame(sigcv,prescv) 
ggplot(data = sigcv.df, aes(x=Var1, y=Var2, fill=prescv)) + 
  geom_tile() +
  theme_minimal() +
  scale_fill_distiller(limits=c(0, 0.35),palette="YlGnBu", direction=-1) +
  labs(x=expression(sigma[0]),y=expression(sigma[1]),fill="p-value")

# linear weight
## CV
cl = snow::makeCluster(ncore-1, type="SOCK")
clusterEvalQ(cl,{library(spatstat);library(tidyverse);library(tools);library(survival)})
source(paste0(getCurrentFileLocation(),"/brain_cv_linear_function.R"))
prescv=snow::parApply(cl,matrix(1:nrow(sigcv),nrow(sigcv),1),1,main,
                      clinical_data_brain,minmax0,minmax1,
                      slideid100,ids,patientid,ntumor,pd0.total,pd1.total) # pvalue results from cv
stopCluster(cl)
save(prescv, file = paste0(getCurrentFileLocation(),
                           "/brain_cv_linear_result.Rdata", sep=""))

## plot
sigcv.df = data.frame(sigcv,prescv) 
ggplot(data = sigcv.df, aes(x=Var1, y=Var2, fill=prescv)) + 
  geom_tile() +
  theme_minimal() +
  scale_fill_distiller(limits=c(0, 0.35),palette="YlGnBu", direction=-1) +
  labs(x=expression(sigma[0]),y=expression(sigma[1]),fill="p-value")

# pwgk
## CV
cl = snow::makeCluster(ncore-1, type="SOCK")
clusterEvalQ(cl,{library(spatstat);library(tidyverse);library(tools);library(survival)})
source(paste0(getCurrentFileLocation(),"/brain_cv_pwgk_function.R"))
prescv=snow::parApply(cl,matrix(1:nrow(sigcv),nrow(sigcv),1),1,main,
                      clinical_data_brain,minmax0,minmax1,
                      slideid100,ids,patientid,ntumor,pd0.total,pd1.total) # pvalue results from cv
stopCluster(cl)
save(prescv, file = paste0(getCurrentFileLocation(),
                           "/brain_cv_pwgk_result.Rdata", sep=""))
## plot
sigcv.df = data.frame(sigcv,prescv) 
ggplot(data = sigcv.df, aes(x=Var1, y=Var2, fill=prescv)) + 
  geom_tile() +
  theme_minimal() +
  scale_fill_distiller(limits=c(0, 0.35),palette="YlGnBu", direction=-1) +
  labs(x=expression(sigma[0]),y=expression(sigma[1]),fill="p-value")

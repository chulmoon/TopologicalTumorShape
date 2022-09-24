library(snow)
library(parallel)

niter = 100 
permres = rep(NA,niter)
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

ncore = parallel::detectCores()
cl = snow::makeCluster(ncore-1, type="SOCK")
snow::clusterEvalQ(cl,{library(spatstat);library(tidyverse);library(tools);library(survival)})
source(paste0(getCurrentFileLocation(),"/lung_simulation_function.R"))
permres=snow::parApply(cl,matrix(1:niter,niter,1),1,main)
snow::stopCluster(cl)

# option to save result
# save(permres, file = paste0(getCurrentFileLocation(),"/lung_simulation_result.Rdata"))

pvals = rep(NA,length(permres))
count = rep(NA,length(permres))

for (ii in 1:length(permres)){
  pvals[ii]=permres[[ii]]$p.val
  count[ii]=permres[[ii]]$coxph_count
}

length(pvals[pvals==pvals[2]]) # cases when the p-values of log rank tests of FCoxPH 
                               # is equal to those of CoxPH  
sum(count)/(230*niter)         # proportion of all models 
                               # that did not select topological features

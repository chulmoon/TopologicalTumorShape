library(snow)
library(parallel)
library(tidyverse)

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

pvals_df=data.frame(pvals)
pvals_df$CoxPH=1-(pvals==pvals[2])*1
pvals_df$CoxPH=pvals_df$CoxPH+(pvals!=pvals[2] & pvals>(pvals[2]+0.0004))
pvals_df = pvals_df %>% mutate(
  cases = case_when(
    CoxPH ==0 ~ "0",
    CoxPH ==1 ~ "<0.0004",
    CoxPH == 2 ~ ">0.0004"))
pvals_df$cases=factor(pvals_df$cases,levels=c("0","<0.0004",">0.0004"))

ggplot(pvals_df)+
  geom_bar(aes(cases,fill=cases)) +
  labs(x="Difference of p-values",y="Frequency",fill="Difference")
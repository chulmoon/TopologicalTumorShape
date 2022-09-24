library(tidyverse)

#############################################################################
# current file location
#############################################################################
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


#############################################################################
# range of persistence function
#############################################################################
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

#############################################################################
# generate persistence function
#############################################################################
# mdw
pfmdw = pfgen(sigma0=1.4,sigma1=0.4,weightfun=mdw,
              filename,minmax0,minmax1,patientid,ids,pd0.total,pd1.total)

# pwgk
pfpwgk = pfgen(sigma0=4.2,sigma1=2.6,weightfun=pwgk,
              filename,minmax0,minmax1,patientid,ids,pd0.total,pd1.total)

# linear
pflinear = pfgen(sigma0=1.2,sigma1=4.6,weightfun=lin,
              filename,minmax0,minmax1,patientid,ids,pd0.total,pd1.total)

#############################################################################
# FCoxPH model
#############################################################################
# clinical information
clinical_data_brain = read_csv(paste0(getCurrentFileLocation(),"/clinical_data_brain.csv"))

## mdw
fcoxph_mdw_model=fcoxph(X0=pfmdw$pfres0,X1=pfmdw$pfres1,ntumor,ids,clinical_data_brain)

summary(fcoxph_mdw_model$model.fcoxph) # FCoxPH
summary(fcoxph_mdw_model$model.coxph)  # CoxPH

## plot estimated coefficient functions
fcoef_mdw=fcoef(model.fcoxph=fcoxph_mdw_model$model.fcoxph,
                Eigvec0=fcoxph_mdw_model$Eigvec0,Eigvec1=fcoxph_mdw_model$Eigvec1,
                minmax0,minmax1,AICselectindex=fcoxph_mdw_model$AICselectindex)

df0=fcoef_mdw[[1]]
df1=fcoef_mdw[[2]]

ggplot(df0,aes(x,y,fill=z))+
  geom_raster()+
  theme_tufte()+
  labs(x="birth",y="death",fill="") +
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  geom_abline(slope=1,intercept=0)+
  scale_fill_gradient2(limits=c( -max(abs(df0$z)),max(abs(df0$z)) ) ) +
  theme(legend.position = c(0.9, 0.25))

ggplot(df1,aes(x,y,fill=z))+
  geom_raster()+
  theme_tufte()+
  labs(x="birth",y="death",fill="") +
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  geom_abline(slope=1,intercept=0)+
  scale_fill_gradient2(limits=c( -max(abs(df1$z)),max(abs(df1$z)) ) ) +
  theme(legend.position = c(0.9, 0.25))

## pwgk
fcoxph_pwgk_model=fcoxph(X0=pfpwgk$pfres0,X1=pfpwgk$pfres1,ntumor,ids,clinical_data_brain)
summary(fcoxph_pwgk_model$model.fcoxph)
summary(fcoxph_pwgk_model$model.coxph)

fcoef_pwgk=fcoef(model.fcoxph=fcoxph_pwgk_model$model.fcoxph,
                Eigvec0=fcoxph_pwgk_model$Eigvec0,Eigvec1=fcoxph_pwgk_model$Eigvec1,
                minmax0,minmax1,AICselectindex=fcoxph_pwgk_model$AICselectindex)

df0=fcoef_pwgk[[1]]
df1=fcoef_pwgk[[2]]

ggplot(df0,aes(x,y,fill=z))+
  geom_raster()+
  theme_tufte()+
  labs(x="birth",y="death",fill="") +
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  geom_abline(slope=1,intercept=0)+
  scale_fill_gradient2(limits=c( -max(abs(df0$z)),max(abs(df0$z)) ) ) +
  theme(legend.position = c(0.9, 0.25))

ggplot(df1,aes(x,y,fill=z))+
  geom_raster()+
  theme_tufte()+
  labs(x="birth",y="death",fill="") +
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  geom_abline(slope=1,intercept=0)+
  scale_fill_gradient2(limits=c( -max(abs(df1$z)),max(abs(df1$z)) ) ) +
  theme(legend.position = c(0.9, 0.25))

## linear
fcoxph_linear_model=fcoxph(X0=pflinear$pfres0,X1=pflinear$pfres1,ntumor,ids,clinical_data_brain)
summary(fcoxph_linear_model$model.fcoxph)
summary(fcoxph_linear_model$model.coxph)

fcoef_linear=fcoef(model.fcoxph=fcoxph_linear_model$model.fcoxph,
                 Eigvec0=fcoxph_linear_model$Eigvec0,Eigvec1=fcoxph_linear_model$Eigvec1,
                 minmax0,minmax1,AICselectindex=fcoxph_linear_model$AICselectindex)

df0=fcoef_linear[[1]]
df1=fcoef_linear[[2]]

ggplot(df0,aes(x,y,fill=z))+
  geom_raster()+
  theme_tufte()+
  labs(x="birth",y="death",fill="") +
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  geom_abline(slope=1,intercept=0)+
  scale_fill_gradient2(limits=c( -max(abs(df0$z)),max(abs(df0$z)) ) ) +
  theme(legend.position = c(0.9, 0.25))

ggplot(df1,aes(x,y,fill=z))+
  geom_raster()+
  theme_tufte()+
  labs(x="birth",y="death",fill="") +
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  geom_abline(slope=1,intercept=0)+
  scale_fill_gradient2(limits=c( -max(abs(df1$z)),max(abs(df1$z)) ) ) +
  theme(legend.position = c(0.9, 0.25))

#################################################################
# validation test using chi-square test
#################################################################
validity_test_mdw=validtest(threshold=0.9,X0=pfmdw$pfres0,X1=pfmdw$pfres1,
                            ntumor,ids,clinical_data_brain)
validity_test_mdw

validity_test_pwgk=validtest(threshold=0.9,X0=pfpwgk$pfres0,X1=pfpwgk$pfres1,
                             ntumor,ids,clinical_data_brain)
validity_test_pwgk

validity_test_linear=validtest(threshold=0.9,X0=pflinear$pfres0,X1=pflinear$pfres1,
                               ntumor,ids,clinical_data_brain)
validity_test_linear

##################################################################
# LOOCV risk prediction
##################################################################

# FCoxPH
cvresult_fcoxph=cvpred_fcoxph(X0=pfmdw$pfres0,X1=pfmdw$pfres1,ntumor,ids,clinical_data_brain)
survdiff(Surv(os,status) ~ group, data = cvresult_fcoxph$risk.mean.fcoxph)

# CoxPH
cvresult_coxph=cvpred_coxph(ntumor,ids,clinical_data_brain)
survdiff(Surv(os,status) ~ group, data = cvresult_coxph$risk.mean.coxph)

## Kaplan-Meier plots
### FCoxPH
ggsurvplot(survfit(Surv(os, status) ~ group, data = cvresult_fcoxph$risk.mean.fcoxph), 
           conf.int = TRUE, pval = T, pval.coord = c(1500, 0.9), 
           xlim = c(0, max(cvresult_fcoxph$risk.mean.fcoxph$os)),
           legend.labs = c("High", "Low"),legend.title="Predicted risk")

### CoxPH
ggsurvplot(survfit(Surv(os, status) ~ group, data = cvresult_coxph$risk.mean.coxph), 
           conf.int = TRUE, pval = T, pval.coord = c(1500, 0.9), 
           xlim = c(0, max(cvresult_coxph$risk.mean.coxph$os)),
           legend.labs = c("High", "Low"),legend.title="Predicted risk")

## Hazard ratios
### FCoxPH
data.survdiff.fcoxph=survdiff(Surv(os, status) ~ group, data = cvresult_fcoxph$risk.mean.fcoxph)
HR.fcoxph = (data.survdiff.fcoxph$obs[1]/data.survdiff.fcoxph$exp[1])/(data.survdiff.fcoxph$obs[2]/data.survdiff.fcoxph$exp[2])
HR.fcoxph

### CoxPH
data.survdiff.coxph=survdiff(Surv(os, status) ~ group, data = cvresult_coxph$risk.mean.coxph)
HR.coxph = (data.survdiff.coxph$obs[1]/data.survdiff.coxph$exp[1])/(data.survdiff.coxph$obs[2]/data.survdiff.coxph$exp[2])
HR.coxph

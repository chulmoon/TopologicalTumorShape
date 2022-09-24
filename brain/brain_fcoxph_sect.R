# Brain FCoxPH with SECT
library(BGLR)
library(doParallel)
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
library(locpol)
library(MASS)
library(survival)
library(survminer)
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

# set the directory path pulled from https://github.com/lorinanthony/SECT
SECTdirectory=".../SECT"

# Load in the C++ BAKR Kernel functions
sourceCpp(paste0(SECTdirectory,"/Software/Kernel_Functions.cpp"))

# Load clinical information
clinical_data_brain = read_csv(paste0(getCurrentFileLocation(),"/clinical_data_brain.csv"))

# Load in the List that holds the Euler Characteristic (EC) Curves for the TCIA Samples
load(paste0(SECTdirectory,"/Data/MRI_ECs.RData"))
nrot = ncol(MRI_list[[1]]$EC)
stepsize = nrow(MRI_list[[1]]$EC)
ECs = matrix(nrow = length(MRI_list),ncol = nrot*stepsize)
rownames(ECs) = 1:nrow(ECs)

# Place the Curves in an nxp Matrix with Patient Tags as the Row Names
for(i in 1:nrow(ECs)){
  ECs[i,] = c(MRI_list[[i]]$EC)
  rownames(ECs)[i] = MRI_list[[i]]$name
}
# Remove the Vectors of Zeros where a New MRI Slice Begins
ECs = ECs[,-seq(101,ncol(ECs),by=101)]

# slide ids
source(paste0(getCurrentFileLocation(),"/brain_fcoxph_functions.R"))
load(paste0(getCurrentFileLocation(),"/slideid100.Rdata"))
# it can be obtained from the code below
#idres=getid(SECTdirectory)
#slideid100=idres$slideid100
#ids=idres$ids

pdrangeres=pdrange(slideid100)
patientid=pdrangeres$patientid
ntumor=pdrangeres$ntumor

# persistence function with the maximum distance weight
tumor_data=data.frame(ntumor=ntumor,patid=ids)

clinical_data_brain = clinical_data_brain %>%
  left_join(tumor_data)

# Find Samples with Clinical Traits 
ECs = ECs[which(rownames(ECs) %in% clinical_data_brain$patid),]

id_new = data.frame(id=rownames(ECs))

# variable selection using AIC
# hh: number of FPCs to examine in AIC
hh=10 # up to 10 FPCs 

# spectral decomposition of X^0
## X^0-mu^0 = U*Sigma*V^T
scaleX0=scale(ECs,center=TRUE,scale=FALSE)
## eigenfunction of X^0-mu^0. phi_j
Eigvec0=svd(scaleX0)$v
## eigenvalue of X^0-mu^0
Singval0=svd(scaleX0)$d
## eigenvalue of (X^0-mu^0)'(X^0-mu^0)
est_Eigval0=Singval0^2

comb=1:hh
loglikelihood=AIC=rep(NA,length(comb))

# select PCs according to AIC
for (ss in 1:length(comb)) {
  # eigenfunction
  est_Eigfun0=Eigvec0[,1:comb[ss]]
  
  # eigenscore
  eigenscore0=scaleX0%*%est_Eigfun0
  
  # designmatrix
  designmatrix=data.frame(eigenscore0,id=id_new)
  colnames(designmatrix)=c(paste0("dim0.",1:comb[ss]),"id")
  
  # combined data
  comdf = designmatrix %>% 
    left_join(clinical_data_brain, c("id"="patid"))

  # FCoxPH model
  # overall survival
  inifor="Surv(os,status)~gender+age+ks+ntumor"
  for (s in 1:comb[ss]){
    inifor=paste0(inifor,paste0("+dim0.",s))
  }
  
  ## FCoxPH model fitting
  model.cox = coxph(as.formula(inifor),data=comdf)
  finalresult=summary(model.cox)
  loglikelihood[ss]=finalresult$loglik[2]
  AIC[ss]=2*(comb[ss])-2*loglikelihood[ss]
}	

# selected number of FPCs by AIC
AICselectindex=which.min(AIC)

# eigenfunction
est_Eigfun0=Eigvec0[,1:comb[AICselectindex]]

# eigenscore
eigenscore0=scaleX0%*%est_Eigfun0

designmatrix=data.frame(eigenscore0,id=id_new)
colnames(designmatrix)=c(paste0("dim0.",1:comb[AICselectindex]),"id")

# combined data
comdf = designmatrix %>% 
  left_join(clinical_data_brain, c("id"="patid"))

# FCoxPH model
## formula for FCoxPH model
inifor="Surv(os,status)~gender+age+ks+ntumor"
for (s in 1:comb[AICselectindex]){
  inifor=paste0(inifor,paste0("+dim0.",s))
}

## FCoxPH model fitting
model.fcoxph = coxph(as.formula(inifor),data=comdf)
summary(model.fcoxph)



##################################################################
# LOOCV risk prediction
##################################################################
X0 = ECs

# predicted risk of FCoxPH for LOOCV
risk.pred.fcoxph = rep(NA,nrow(X0))
# predicted risk of CoxPH for LOOCV
risk.pred.coxph = rep(NA,nrow(X0))

# LOOCV
hh=10
for (ii in 1:nrow(X0)){
  # train and test data
  X0.train = X0[-ii,]
  X0.test = X0[ii,]
  id.train = rownames(X0.train)
  id.test = rownames(data.frame(X0)[ii,])
  
  # spectral decomposition of X^0
  ## X^0-mu^0 = U*Sigma*V^T
  scaleX0.train=scale(X0.train,center=TRUE,scale=FALSE)
  ## eigenfunction of X^0-mu^0. phi_j
  Eigvec0=svd(scaleX0.train)$v
  ## eigenvalue of X^0-mu^0
  #Singval0=svd(scaleX0)$d
  ## eigenvalue of (X^0-mu^0)'(X^0-mu^0)
  #est_Eigval0=Singval0^2
  
  comb=1:hh
  loglikelihood=AIC=rep(NA,length(comb))
  
  # select PCs according to AIC
  for (ss in 1:length(comb)) {
    # eigenfunction
    est_Eigfun0=Eigvec0[,1:comb[ss]]
    
    # eigenscore
    eigenscore0=scaleX0.train%*%est_Eigfun0
    
    # designmatrix
    designmatrix.train=data.frame(eigenscore0,id=id.train)
    colnames(designmatrix.train)=c(paste0("dim0.",1:comb[ss]),"id")
    
    # combined data
    comdf.train = designmatrix.train %>% 
      left_join(clinical_data_brain, c("id"="patid"))
    
    # FCoxPH model
    ## formula for FCoxPH model
    inifor="Surv(os,status)~gender+age+ks+ntumor"
    for (s in 1:comb[ss]){
      inifor=paste0(inifor,paste0("+dim0.",s))
    }
    
    ## FCoxPH model fitting
    model.cox = coxph(as.formula(inifor),data=comdf.train)
    finalresult=summary(model.cox)
    loglikelihood[ss]=finalresult$loglik[2]
    AIC[ss]=2*(comb[ss])-2*loglikelihood[ss]
  }
  
  # selected number of FPCs by AIC
  AICselectindex=which.min(AIC)
  
  # eigenfunction
  est_Eigfun0=as.matrix(Eigvec0[,1:comb[AICselectindex]])
  
  # eigenscore
  eigenscore0=scaleX0.train%*%est_Eigfun0
  
  designmatrix.train=data.frame(eigenscore0,id=id.train)
  colnames(designmatrix.train)=c(paste0("dim0.",1:comb[AICselectindex])
                                 ,"id")
  
  # combined data
  comdf.train = designmatrix.train %>% 
    left_join(clinical_data_brain, c("id"="patid"))
  
  # FCoxPH model
  ## formula for FCoxPH model
  inifor="Surv(os,status)~gender+age+ks+ntumor"
  for (s in 1:comb[AICselectindex]){
    inifor=paste0(inifor,paste0("+dim0.",s))
  }
  
  ## FCoxPH model fitting
  model.fcoxph = coxph(as.formula(inifor),data=comdf.train)
  
  # test data 
  ## substract mean of training data
  scaleX0.test = t(as.matrix(X0.test-colMeans(X0.train)))
  ## estimated eigenscore of test data
  eigenscore0.test=as.matrix(scaleX0.test)%*%est_Eigfun0
  ## design matrix of test data
  designmatrix.test=data.frame(eigenscore0.test,id=id.test)
  colnames(designmatrix.test)=c(paste0("dim0.",1:comb[AICselectindex]),"id")
  ## combined data
  comdf.test = designmatrix.test %>% 
    left_join(clinical_data_brain, c("id"="patid"))
  
  ## risk prediction
  risk.pred.fcoxph[ii]=predict(model.fcoxph,comdf.test,type="risk")
}

##################################################################
# Kaplan-Meier plots
##################################################################
# FcoxPH
## risk data
risk.df.fcoxph=data.frame(id=rownames(X0),risk=risk.pred.fcoxph)
risk.mean.fcoxph = comdf %>% 
  left_join(risk.df.fcoxph, by=c("id"="id"))
risk.mean.fcoxph$group=ifelse(risk.mean.fcoxph$risk<=median(risk.mean.fcoxph$risk),"low","high")

risk.result = survival::survdiff(Surv(os, status) ~ group, data = risk.mean.fcoxph)
rescv = pchisq(risk.result$chisq, length(risk.result$n)-1, lower.tail = FALSE)
rescv
survdiff(Surv(os,status) ~ group, data = risk.mean.fcoxph)

## KM plot
ggsurvplot(survfit(Surv(os, status) ~ group, data = risk.mean.fcoxph), 
           conf.int = TRUE, pval = T, pval.coord = c(1500, 0.9), 
           xlim = c(0, max(risk.mean.fcoxph$os)),
           legend.labs = c("High", "Low"),legend.title="Predicted risk")

##################################################################
# Hazard ratios
##################################################################
# FCoxPH
## difference of curves in KM plot
data.survdiff.fcoxph=survdiff(Surv(os, status) ~ group, data = risk.mean.fcoxph)
## hazard ratio
HR.fcoxph = (data.survdiff.fcoxph$obs[1]/data.survdiff.fcoxph$exp[1])/(data.survdiff.fcoxph$obs[2]/data.survdiff.fcoxph$exp[2])
HR.fcoxph
























# Setup Parameters for Analysis
# Load in the Cross-Validated Bandwidths
load("C:/Users/Chul/Box/r_Papers/ShapeAnalysis/SECT-master/Analysis/Cross_Validation_Results/GaussCV_Results.RData");
# Call the Spectrum of Considered Bandwidths 
theta = seq(from = 0.1, to = 10, by = 0.1)


##################################################################
# LOOCV (scaled) survival prediction for Gaussian Process (GP) model
##################################################################
# use scaled overall survival
y = scale(clinical.info$os)
# predicted survival
surv.pred = rep(NA,nrow(clinical.info))

j=1 # overall survival
for (ii in 1:length(y)){
  print(ii)
  # Create the training and test sets
  ind = 1:length(y)
  ind = ind[-ii]
  
  # Topological Summary Statistics Analysis
  X = ECs
  X = scale(X)
  X = cbind(rep(1,nrow(X)),X)
  theta_hat = theta[cvs[[j]][4]]
  K = GaussKernel(t(X),theta_hat)
  diag(K) = 1
  
  # Center and Scale the Covariance Matrix
  n=nrow(K)
  v=matrix(1, n, 1)
  M=diag(n)-v%*%t(v)/n
  K=M%*%K%*%M
  K=K/mean(diag(K))
  
  # Compute the Posterior Predictive Mean
  Kn = K[ind,ind]
  fhat = K[,ind] %*% solve(Kn + diag(sigma,nrow(Kn)),y[ind]) # predicted survival
  surv.pred[ii]=fhat[ii]
}

##################################################################
# Kaplan-Meier plots
##################################################################
# Gaussian Process
## risk data
risk.fcoxph.sect=data.frame(id=clinical.info$patid,surv.time=surv.pred)
risk.fcoxph.sect.joined = risk.fcoxph.sect %>%
  left_join(clinical.info.brain, c("id"="patid"))
risk.fcoxph.sect.joined$group=
  ifelse(risk.fcoxph.sect.joined$surv.time<=median(risk.fcoxph.sect.joined$surv.time)," Short","Long")

survdiff(Surv(os, event) ~ group, data = risk.fcoxph.sect.joined)

# FCoxPH
## difference of curves in KM plot
data.survdiff.fcoxph.sect=survdiff(Surv(os, status) ~ group, data = risk.fcoxph.sect.joined)
## hazard ratio
HR.fcoxph = (data.survdiff.fcoxph$obs[1]/data.survdiff.fcoxph$exp[1])/(data.survdiff.fcoxph$obs[2]/data.survdiff.fcoxph$exp[2])
HR.fcoxph

## KM plot

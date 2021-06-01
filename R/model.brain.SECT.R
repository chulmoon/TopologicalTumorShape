library(BGLR)
library(doParallel)
library(Rcpp)
library(RcppArmadillo)
library(RcppParallel)
library(locpol)
library(MASS)
library(survival)
library(survcomp)
library(survminer)
library(tidyverse)

# Load in the C++ BAKR Kernel functions from https://github.com/lorinanthony/SECT
sourceCpp("./SECT/Software/Kernel_Functions.cpp")

# Load in the List that holds the Euler Characteristic (EC) Curves for the TCIA Samples
# from https://github.com/lorinanthony/SECT
load("./SECT/Data/MRI_ECs.RData")
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


# Load clinical information 
# from https://github.com/chulmoon/TopologicalTumorShape
load("./TopologicalTumorShape/data/brain/clinical_info_brain.Rdata")

# Find Samples with Clinical Traits 
ECs = ECs[which(rownames(ECs)%in%clinical.info.brain$patid),]

id.info = data.frame(id=rownames(ECs))

clinical.info= id.info %>%
	left_join(clinical.info.brain,c("id"="patid")) %>%
	rename(patid=id)

# Setup Parameters for Analysis
# Load in the Cross-Validated Bandwidths
# from https://github.com/lorinanthony/SECT
load("./SECT/Analysis/Cross_Validation_Results/GaussCV_Results.RData");
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
risk.gp=data.frame(id=clinical.info$patid,surv.time=surv.pred)
risk.gp.joined = risk.gp %>%
	left_join(clinical.info.brain, c("id"="patid"))
risk.gp.joined$group=
	ifelse(risk.gp.joined$surv.time<=median(risk.gp.joined$surv.time)," Short","Long")

survdiff(Surv(os, event) ~ group, data = risk.gp.joined)

## KM plot
ggsurvplot(survfit(Surv(os, event) ~ group, data = risk.gp.joined), 
					 conf.int = TRUE,pval = T, pval.coord = c(1000, 0.9),
					 legend.labs = c("Short", "Long"),legend.title="Predicted survival")

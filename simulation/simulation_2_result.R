library(spatstat)
library(tidyverse)
library(tools)
library(survival)
library(ggthemes)

iter=30

minmax0 = c(-80,80)
minmax1 = c(-20,80)

beta0 = matrix(NA,30,(80+80)^2)
beta1 = matrix(NA,30,(100)^2)

set.seed(40)

validity.pval = wald.fcoxph.pval = wald.coxph.pval = rep(NA,iter)

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

for (tt in 1:iter){
  print(tt)
  # load PD
  filename_control = list.files(path=paste0(getCurrentFileLocation(),"/simulation_2_pd")
                                ,pattern=paste0("control_iter_",tt-1,"_"),full.names = T)
  filename_test = list.files(path=paste0(getCurrentFileLocation(),"/simulation_2_pd")
                             ,pattern=paste0("feature_iter_",tt-1,"_"),full.names = T)
  filename = c(filename_control,filename_test)
  # recover slide id
  pd.total = NULL
  for (ii in 1:length(filename)) {
    pd=read.table(filename[ii],col.names = c("dim","birth","death"))
    pd$group = substr(basename(filename[ii]),1,1)
    pd$id = strtoi(strsplit(file_path_sans_ext(basename(filename[ii])),"_")[[1]][5])
    pd = pd %>% filter(!(death == Inf & dim==1)) # remove dim1 Inf death values
    pd[pd$death==Inf,3]=pd[pd$death==Inf,2] #replace Inf death with birth value
    pd.total = rbind(pd.total,pd)
  }
  
  
  pd0.total = pd.total %>% filter(dim==0)
  pd1.total = pd.total %>% filter(dim==1)
  
  # number of pixels
  pixnum0 = ceiling(minmax0[2])-floor(minmax0[1])
  pixnum1 = ceiling(minmax1[2])-floor(minmax1[1])
  
  gid=unique(pd.total[,4:5]) %>%
    arrange(group,id)
  
  # Persistence function
  # persistence functions with the maximum distance weight
  pfmdw0 = as.data.frame(matrix(NA,length(filename),pixnum0*(pixnum0+1)/2))
  pfmdw1 = as.data.frame(matrix(NA,length(filename),pixnum1*(pixnum1+1)/2))
  
  mdw = function(pd) {
    mdw=apply(cbind(pd$death-pd$birth,abs(pd$birth),abs(pd$death)),1,max)
    return(mdw)
  }
  
  sigma0=2
  sigma1=2
  
  
  for (ii in 1:nrow(gid)){
    # load aggregated persistence diagrams
    pd0=pd0.total %>% filter(group==gid$group[ii],id==gid$id[ii])
    pd1=pd1.total %>% filter(group==gid$group[ii],id==gid$id[ii])
    
    if (nrow(pd1)==0){
      ########## dimension 0 persistence function ########## 
      dim0.x = pd0$birth
      dim0.y = pd0$death
      
      # ppp object
      pd0.ppp = ppp(dim0.x, dim0.y, minmax0, minmax0)
      
      # maximum distance weight for dimension-zero
      mdw0 = mdw(pd0)
      
      # compute persistence function
      pf0=density(pd0.ppp,sigma=sigma0,dimyx=c(pixnum0,pixnum0),
                  weights=mdw0)
      
      # for center of the pixel
      cunit0=(minmax0[2]-minmax0[1])/(2*pixnum0)
      # grids for birth and death
      gridx0=seq(minmax0[1]+cunit0,minmax0[2]-cunit0,length.out = pixnum0)
      gridy0=seq(minmax0[1]+cunit0,minmax0[2]-cunit0,length.out = pixnum0)
      
      # persistence function
      pf0.mdw=pf0$v
      df0=data.frame(x=rep(gridx0,each=pixnum0),
                     y=rep(gridy0,pixnum0),
                     z=as.vector(pf0.mdw))
      # some pixels are negative, so convert them to zero
      dfmw0 = df0 %>% 
        filter(y>=x) %>%
        mutate(newz=case_when(
          z < 0 ~ 0,
          TRUE ~ z)
        )
      pfmdw0[ii,] = dfmw0$newz # persistence function	
      
      ########## dimension 1 persistence function ########## 
      cunit1=(minmax1[2]-minmax1[1])/(2*pixnum1)
      gridx1=seq(minmax1[1]+cunit1,minmax1[2]-cunit1,length.out = pixnum1)
      gridy1=seq(minmax1[1]+cunit1,minmax1[2]-cunit1,length.out = pixnum1)
      
      # maximum distance weight
      df1=data.frame(x=rep(gridx1,each=pixnum1),
                     y=rep(gridy1,pixnum1),
                     z=0)
      dfmdw1 = df1 %>% 
        filter(y>=x) %>%
        mutate(newz=case_when(
          z < 0 ~ 0,
          TRUE ~ z)
        )
      pfmdw1[ii,] = dfmdw1$newz # persistence function
    } else{
      ########## dimension 0 persistence function ########## 
      dim0.x = pd0$birth
      dim0.y = pd0$death
      
      # ppp object
      pd0.ppp = ppp(dim0.x, dim0.y, minmax0, minmax0)
      
      # maximum distance weight for dimension-zero
      mdw0 = mdw(pd0)
      
      # compute persistence function
      pf0=density(pd0.ppp,sigma=sigma0,dimyx=c(pixnum0,pixnum0),
                  weights=mdw0)
      
      # for center of the pixel
      cunit0=(minmax0[2]-minmax0[1])/(2*pixnum0)
      # grids for birth and death
      gridx0=seq(minmax0[1]+cunit0,minmax0[2]-cunit0,length.out = pixnum0)
      gridy0=seq(minmax0[1]+cunit0,minmax0[2]-cunit0,length.out = pixnum0)
      
      # persistence function
      pf0.mdw=pf0$v
      df0=data.frame(x=rep(gridx0,each=pixnum0),
                     y=rep(gridy0,pixnum0),
                     z=as.vector(pf0.mdw))
      # some pixels are negative, so convert them to zero
      dfmw0 = df0 %>% 
        filter(y>=x) %>%
        mutate(newz=case_when(
          z < 0 ~ 0,
          TRUE ~ z)
        )
      pfmdw0[ii,] = dfmw0$newz # persistence function	
      
      ########## dimension 1 persistence function ########## 
      dim1.x = pd1$birth
      dim1.y = pd1$death
      
      # ppp object
      pd1.ppp = ppp(dim1.x, dim1.y, minmax1, minmax1)
      
      # weights
      mdw1 = mdw(pd1)
      
      # compute persistence function
      pf1=density(pd1.ppp,sigma=sigma1,dimyx=c(pixnum1,pixnum1),
                  weights=mdw1)
      
      cunit1=(minmax1[2]-minmax1[1])/(2*pixnum1)
      gridx1=seq(minmax1[1]+cunit1,minmax1[2]-cunit1,length.out = pixnum1)
      gridy1=seq(minmax1[1]+cunit1,minmax1[2]-cunit1,length.out = pixnum1)
      
      # maximum distance weight
      pf1.mdw=pf1$v
      df1=data.frame(x=rep(gridx1,each=pixnum1),
                     y=rep(gridy1,pixnum1),
                     z=as.vector(pf1.mdw))
      dfmdw1 = df1 %>% 
        filter(y>=x) %>%
        mutate(newz=case_when(
          z < 0 ~ 0,
          TRUE ~ z)
        )
      pfmdw1[ii,] = dfmdw1$newz # persistence function
    }
  }
  
  
  # persistence function with the maximum distance weight
  X0=pfmdw0
  X1=pfmdw1
  
  ####################### survival data generation
  
  # baseline hazard: Weibull
  
  # N = sample size    
  # lambda = scale parameter in h0()
  # rho = shape parameter in h0()
  # beta = fixed effect parameter
  # rateC = rate parameter of the exponential distribution of C
  
  simulWeib <- function(N, lambda, rho, beta, rateC){
    # covariate --> N Bernoulli trials
    #age = rpois(N,40)
    x = c(rep(0,N/2),rep(1,N/2))
    
    # Weibull latent event times
    v <- runif(n=N)
    Tlat <- (- log(v) / (lambda * exp(x * beta)))^(1 / rho)
    
    # censoring times
    C <- rexp(n=N, rate=rateC)
    
    # follow-up times and event indicators
    time <- pmin(Tlat, C)
    status <- as.numeric(Tlat <= C)
    
    # data set
    return(data.frame(id=1:N,time=time,status=status,x=x))
  }
  
  
  N=200
  dat <- simulWeib(N=N, lambda=0.01, rho=1, beta=0.8, rateC=0.001)
  dat$age = rpois(N,40)
  dat$sex = sample(x=c(0, 1), size=N, replace=TRUE, prob=c(0.5, 0.5))
  
  # combined data
  comdf = dat
  # CoxPH model
  ## formula 
  inifor="Surv(time, status)~age+sex"
  ## CoxPH model fitting
  model.coxph = coxph(as.formula(inifor),data=comdf)
  
  
  # spectral decomposition of X^0
  ## X^0-mu^0 = U*Sigma*V^T
  scaleX0=scale(X0,center=TRUE,scale=FALSE)
  ## eigenfunction of X^0-mu^0. phi_j
  Eigvec0=svd(scaleX0)$v
  ## eigenvalue of X^0-mu^0
  Singval0=svd(scaleX0)$d
  ## eigenvalue of (X^0-mu^0)'(X^0-mu^0)
  est_Eigval0=Singval0^2
  
  # spectral decomposition of X^1
  ## X^1-mu^1 = U*Sigma*V^T
  scaleX1=scale(X1,center=TRUE,scale=FALSE)
  ## eigenfunction of X^1-mu^1. pi_k
  Eigvec1=svd(scaleX1)$v
  ## eigenvalue of X^1-mu^10.00
  Singval1=svd(scaleX1)$d
  ## eigenvalue of (X^1-mu^1)'(X^1-mu^1)
  est_Eigval1=Singval1^2
  
  hh=5
  comb=expand.grid(0:hh,0:hh)
  comb=comb[-1,]
  loglikelihood=AIC=rep(NA,nrow(comb))
  
  
  # select PCs according to AIC
  for (ss in 1:nrow(comb)) {
    #print(ss)
    
    if(comb[ss,1]==0){
      
      # eigenfunction
      est_Eigfun1=Eigvec1[,1:comb[ss,2]]
      
      # eigenscore
      eigenscore1=scaleX1%*%est_Eigfun1
      
      # combined data
      comdf=cbind(eigenscore1,dat)
      colnames(comdf)=c(paste0("dim1.",1:comb[ss,2]),colnames(dat))
      
      # FCoxPH model
      ## formula for FCoxPH model
      inifor="Surv(time, status)~age+sex"
      for (s in 1:comb[ss,2]){
        inifor=paste0(inifor,paste0("+dim1.",s))
      }
      
      ## FCoxPH model fitting
      model.fcox = coxph(as.formula(inifor),data=comdf)
      finalresult=summary(model.fcox)
      loglikelihood[ss]=finalresult$loglik[2]
      AIC[ss]=2*(comb[ss,1]+comb[ss,2])-2*loglikelihood[ss]
      
    } else if (comb[ss,2]==0){
      
      # eigenfunction
      est_Eigfun0=Eigvec0[,1:comb[ss,1]]
      
      # eigenscore
      eigenscore0=scaleX0%*%est_Eigfun0
      
      # combined data
      comdf=cbind(eigenscore0,dat)
      colnames(comdf)=c(paste0("dim0.",1:comb[ss,1]),colnames(dat))
      
      # FCoxPH model
      ## formula for FCoxPH model
      inifor="Surv(time, status)~age+sex"
      for (s in 1:comb[ss,1]){
        inifor=paste0(inifor,paste0("+dim0.",s))
      }
      
      ## FCoxPH model fitting
      model.fcox = coxph(as.formula(inifor),data=comdf)
      finalresult=summary(model.fcox)
      loglikelihood[ss]=finalresult$loglik[2]
      AIC[ss]=2*(comb[ss,1]+comb[ss,2])-2*loglikelihood[ss]
      
    } else{
      # eigenfunction
      est_Eigfun0=Eigvec0[,1:comb[ss,1]]
      est_Eigfun1=Eigvec1[,1:comb[ss,2]]
      
      # eigenscore
      eigenscore0=scaleX0%*%est_Eigfun0
      eigenscore1=scaleX1%*%est_Eigfun1
      
      # combined data
      comdf=cbind(eigenscore0,eigenscore1,dat)
      colnames(comdf)=c(paste0("dim0.",1:comb[ss,1]),paste0("dim1.",1:comb[ss,2]),colnames(dat))
      
      # FCoxPH model
      ## formula for FCoxPH model
      inifor="Surv(time, status)~age+sex"
      for (s in 1:comb[ss,1]){
        inifor=paste0(inifor,paste0("+dim0.",s))
      }
      for (s in 1:comb[ss,2]){
        inifor=paste0(inifor,paste0("+dim1.",s))
      }
      
      ## FCoxPH model fitting
      model.fcox = coxph(as.formula(inifor),data=comdf)
      finalresult=summary(model.fcox)
      loglikelihood[ss]=finalresult$loglik[2]
      AIC[ss]=2*(comb[ss,1]+comb[ss,2])-2*loglikelihood[ss]
    }
  }
  
  # selected number of FPCs by AIC
  AICselectindex=which.min(AIC)
  
  if(comb[AICselectindex,1]==0){
    
    # eigenfunction
    est_Eigfun1=Eigvec1[,1:comb[AICselectindex,2]]
    
    # eigenscore
    eigenscore1=scaleX1%*%est_Eigfun1
    
    # combined data
    comdf=cbind(eigenscore1,dat)
    colnames(comdf)=c(paste0("dim1.",1:comb[AICselectindex,2]),colnames(dat))
    
    # FCoxPH model
    ## formula for FCoxPH model
    inifor="Surv(time, status)~age+sex"
    for (s in 1:comb[AICselectindex,2]){
      inifor=paste0(inifor,paste0("+dim1.",s))
    }
    
    ## FCoxPH model fitting
    model.fcox = coxph(as.formula(inifor),data=comdf)
    summary(model.fcox)
    
  } else if (comb[AICselectindex,2]==0){
    
    # eigenfunction
    est_Eigfun0=Eigvec0[,1:comb[AICselectindex,1]]
    
    # eigenscore
    eigenscore0=scaleX0%*%est_Eigfun0
    
    # combined data
    comdf=cbind(eigenscore0,dat)
    colnames(comdf)=c(paste0("dim0.",1:comb[AICselectindex,1]),colnames(dat))
    
    # FCoxPH model
    ## formula for FCoxPH model
    inifor="Surv(time, status)~age+sex"
    for (s in 1:comb[AICselectindex,1]){
      inifor=paste0(inifor,paste0("+dim0.",s))
    }
    
    ## FCoxPH model fitting
    model.fcox = coxph(as.formula(inifor),data=comdf)
    summary(model.fcox)
    
  } else{
    # eigenfunction
    est_Eigfun0=Eigvec0[,1:comb[AICselectindex,1]]
    est_Eigfun1=Eigvec1[,1:comb[AICselectindex,2]]
    
    # eigenscore
    eigenscore0=scaleX0%*%est_Eigfun0
    eigenscore1=scaleX1%*%est_Eigfun1
    
    # combined data
    comdf=cbind(eigenscore0,eigenscore1,dat)
    colnames(comdf)=c(paste0("dim0.",1:comb[AICselectindex,1]),paste0("dim1.",1:comb[AICselectindex,2]),colnames(dat))
    
    # FCoxPH model
    ## formula for FCoxPH model
    inifor="Surv(time, status)~age+sex"
    for (s in 1:comb[AICselectindex,1]){
      inifor=paste0(inifor,paste0("+dim0.",s))
    }
    for (s in 1:comb[AICselectindex,2]){
      inifor=paste0(inifor,paste0("+dim1.",s))
    }
    
    ## FCoxPH model fitting
    model.fcox = coxph(as.formula(inifor),data=comdf)
    summary(model.fcox)
    
  }
  
  num.clinical = 2 # number of clinical variables
  
  if (comb[AICselectindex,1]==0) {
    beta.est1 = as.matrix(Eigvec1[,1:comb[AICselectindex,2]]) %*% as.matrix(model.fcox$coefficients[(num.clinical+comb[AICselectindex,1]+1):(num.clinical+comb[AICselectindex,1]+comb[AICselectindex,2])] )
    
    minmax1[2]=ceiling(minmax1[2])
    minmax1[1]=floor(minmax1[1])
    
    res1 = minmax1[2]-minmax1[1]
    
    cunit1=(minmax1[2]-minmax1[1])/(2*res1)
    gridx1=seq(minmax1[1]+cunit1,minmax1[2]-cunit1,length.out = res1)
    gridy1=seq(minmax1[1]+cunit1,minmax1[2]-cunit1,length.out = res1)
    
    # data for dimension zero
    beta0[tt,]=NA
    
    # data for dimension one
    df1=data.frame(x=rep(gridx1,each=res1),y=rep(gridy1,res1))
    df1_part = data.frame(x=rep(gridx1,each=res1),y=rep(gridy1,res1)) %>%
      filter(y>=x)
    df1_part$z=beta.est1[1:nrow(df1_part),]
    df1_join = left_join(df1,df1_part,by=c("x"="x","y"="y"))
    df1_join$z = ifelse(is.na(df1_join$z),0,df1_join$z)
    beta1[tt,]=df1_join$z
    
  } else if (comb[AICselectindex,2]==0) {
    beta.est0 = as.matrix(Eigvec0[,1:comb[AICselectindex,1]]) %*% as.matrix(model.fcox$coefficients[(num.clinical+1):(num.clinical+comb[AICselectindex,1])] )
    
    minmax0[2]=ceiling(minmax0[2])
    minmax0[1]=floor(minmax0[1])
    
    
    res0 = minmax0[2]-minmax0[1]
    
    cunit0=(minmax0[2]-minmax0[1])/(2*res0)
    gridx0=seq(minmax0[1]+cunit0,minmax0[2]-cunit0,length.out = res0)
    gridy0=seq(minmax0[1]+cunit0,minmax0[2]-cunit0,length.out = res0)
    
    # data for dimension zero
    df0=data.frame(x=rep(gridx0,each=res0),y=rep(gridy0,res0))  
    df0_part = data.frame(x=rep(gridx0,each=res0),y=rep(gridy0,res0)) %>%
      filter(y>=x)
    df0_part$z=beta.est0[1:nrow(df0_part),]
    df0_join = left_join(df0,df0_part,by=c("x"="x","y"="y"))
    df0_join$z = ifelse(is.na(df0_join$z),0,df0_join$z)
    beta0[tt,]=df0_join$z
    
    beta1[tt,]=NA
  } else{
    beta.est0 = as.matrix(Eigvec0[,1:comb[AICselectindex,1]]) %*% as.matrix(model.fcox$coefficients[(num.clinical+1):(num.clinical+comb[AICselectindex,1])] )
    beta.est1 = as.matrix(Eigvec1[,1:comb[AICselectindex,2]]) %*% as.matrix(model.fcox$coefficients[(num.clinical+comb[AICselectindex,1]+1):(num.clinical+comb[AICselectindex,1]+comb[AICselectindex,2])] )
    
    minmax0[2]=ceiling(minmax0[2])
    minmax0[1]=floor(minmax0[1])
    
    minmax1[2]=ceiling(minmax1[2])
    minmax1[1]=floor(minmax1[1])
    
    res0 = minmax0[2]-minmax0[1]
    res1 = minmax1[2]-minmax1[1]
    
    cunit0=(minmax0[2]-minmax0[1])/(2*res0)
    gridx0=seq(minmax0[1]+cunit0,minmax0[2]-cunit0,length.out = res0)
    gridy0=seq(minmax0[1]+cunit0,minmax0[2]-cunit0,length.out = res0)
    
    cunit1=(minmax1[2]-minmax1[1])/(2*res1)
    gridx1=seq(minmax1[1]+cunit1,minmax1[2]-cunit1,length.out = res1)
    gridy1=seq(minmax1[1]+cunit1,minmax1[2]-cunit1,length.out = res1)
    
    # data for dimension zero
    df0=data.frame(x=rep(gridx0,each=res0),y=rep(gridy0,res0))  
    df0_part = data.frame(x=rep(gridx0,each=res0),y=rep(gridy0,res0)) %>%
      filter(y>=x)
    df0_part$z=beta.est0[1:nrow(df0_part),]
    df0_join = left_join(df0,df0_part,by=c("x"="x","y"="y"))
    df0_join$z = ifelse(is.na(df0_join$z),0,df0_join$z)
    beta0[tt,]=df0_join$z
    
    # data for dimension one
    df1=data.frame(x=rep(gridx1,each=res1),y=rep(gridy1,res1))
    df1_part = data.frame(x=rep(gridx1,each=res1),y=rep(gridy1,res1)) %>%
      filter(y>=x)
    df1_part$z=beta.est1[1:nrow(df1_part),]
    df1_join = left_join(df1,df1_part,by=c("x"="x","y"="y"))
    df1_join$z = ifelse(is.na(df1_join$z),0,df1_join$z)
    beta1[tt,]=df1_join$z
  }
  
  #################################################################
  # validation test using chi-square test
  #################################################################
  # threshold values
  thresholdvalues = c(0.80, 0.85, 0.90, 0.95, 0.99)
  j=4 # use 90%
  threshold=thresholdvalues[j]
  
  # select dimension one FPCs according to the variability threshold
  for (ss in 1:length(est_Eigval0)) {
    if (sum(est_Eigval0[1:ss])/sum(est_Eigval0)>threshold)
      break
  }
  selectindex0=ss
  # select dimension one FPCs according to the variability threshold
  for (qq in 1:length(est_Eigval1)) {
    if (sum(est_Eigval1[1:qq])/sum(est_Eigval1)>threshold)
      break
  }
  selectindex1=qq
  
  est_Eigfun0_var=Eigvec0[,1:selectindex0]
  est_Eigfun1_var=Eigvec1[,1:selectindex1]
  
  eigenscore0_var=scaleX0%*%est_Eigfun0_var
  eigenscore1_var=scaleX1%*%est_Eigfun1_var
  
  comdf_var=cbind(eigenscore0_var,eigenscore1_var,dat)
  colnames(comdf_var)=c(paste0("dim0.", 1:selectindex0),paste0("dim1.", 1:selectindex1),colnames(dat))
  
  # FCoxPH model variable selected by variance threshold
  ## formula
  inifor="Surv(time, status)~age+sex"
  for (s in 1:selectindex0){
    inifor=paste0(inifor,paste0("+dim0.",s))
  }
  for (s in 1:selectindex1){
    inifor=paste0(inifor,paste0("+dim1.",s))
  }
  
  ## FCoxPH model fitting
  model.fcoxph.var = coxph(as.formula(inifor),data=comdf_var)
  finalresult=summary(model.fcoxph.var)
  
  ## initial value for chi-sq test
  initial.value=c(model.coxph$coef,rep(0,selectindex0),rep(0,selectindex1) )
  
  ## coefficient for full model
  testscore=model.fcoxph.var$score
  
  # chi-sq test
  validity.pval[tt] =1-pchisq(testscore, selectindex0+selectindex1, lower.tail = TRUE, log.p = FALSE)
  
  coxphsum = summary(model.coxph)
  fcoxphsum = summary(model.fcox)
  
  wald.coxph.pval[tt]=coxphsum$waldtest[3]
  wald.fcoxph.pval[tt]=fcoxphsum$waldtest[3]
}

mean(validity.pval)
mean(wald.coxph.pval)
mean(wald.fcoxph.pval)

minmax0[2]=ceiling(minmax0[2])
minmax0[1]=floor(minmax0[1])

minmax1[2]=ceiling(minmax1[2])
minmax1[1]=floor(minmax1[1])

res0 = minmax0[2]-minmax0[1]
res1 = minmax1[2]-minmax1[1]

cunit0=(minmax0[2]-minmax0[1])/(2*res0)
gridx0=seq(minmax0[1]+cunit0,minmax0[2]-cunit0,length.out = res0)
gridy0=seq(minmax0[1]+cunit0,minmax0[2]-cunit0,length.out = res0)

cunit1=(minmax1[2]-minmax1[1])/(2*res1)
gridx1=seq(minmax1[1]+cunit1,minmax1[2]-cunit1,length.out = res1)
gridy1=seq(minmax1[1]+cunit1,minmax1[2]-cunit1,length.out = res1)

# data for dimension zero
df0=data.frame(x=rep(gridx0,each=res0),y=rep(gridy0,res0))  
df0_part = data.frame(x=rep(gridx0,each=res0),y=rep(gridy0,res0)) %>%
  filter(y>=x)
df0_part$z=beta.est0[1:nrow(df0_part),]
df0_join = left_join(df0,df0_part,by=c("x"="x","y"="y"))
df0_join$z = colMeans(beta0,na.rm = T)

# data for dimension one
df1=data.frame(x=rep(gridx1,each=res1),y=rep(gridy1,res1))
df1_part = data.frame(x=rep(gridx1,each=res1),y=rep(gridy1,res1)) %>%
  filter(y>=x)
df1_part$z=beta.est1[1:nrow(df1_part),]
df1_join = left_join(df1,df1_part,by=c("x"="x","y"="y"))
df1_join$z = colMeans(beta1,na.rm = T)


minmax0[2]=ceiling(minmax0[2])
minmax0[1]=floor(minmax0[1])

minmax1[2]=ceiling(minmax1[2])
minmax1[1]=floor(minmax1[1])

res0 = minmax0[2]-minmax0[1]
res1 = minmax1[2]-minmax1[1]

cunit0=(minmax0[2]-minmax0[1])/(2*res0)
gridx0=seq(minmax0[1]+cunit0,minmax0[2]-cunit0,length.out = res0)
gridy0=seq(minmax0[1]+cunit0,minmax0[2]-cunit0,length.out = res0)

cunit1=(minmax1[2]-minmax1[1])/(2*res1)
gridx1=seq(minmax1[1]+cunit1,minmax1[2]-cunit1,length.out = res1)
gridy1=seq(minmax1[1]+cunit1,minmax1[2]-cunit1,length.out = res1)

# data for dimension zero
df0=data.frame(x=rep(gridx0,each=res0),y=rep(gridy0,res0))  
df0_part = data.frame(x=rep(gridx0,each=res0),y=rep(gridy0,res0)) %>%
  filter(y>=x)
df0_part$z=beta.est0[1:nrow(df0_part),]
df0_join = left_join(df0,df0_part,by=c("x"="x","y"="y"))
df0_join$z = colMeans(beta0,na.rm = T)

# data for dimension one
df1=data.frame(x=rep(gridx1,each=res1),y=rep(gridy1,res1))
df1_part = data.frame(x=rep(gridx1,each=res1),y=rep(gridy1,res1)) %>%
  filter(y>=x)
df1_part$z=beta.est1[1:nrow(df1_part),]
df1_join = left_join(df1,df1_part,by=c("x"="x","y"="y"))
df1_join$z = colMeans(beta1,na.rm = T)


# Plot estimated coefficient function for dimension zero
ggplot(df0_join,aes(x,y,fill=z))+
  geom_raster()+
  theme_tufte()+
  labs(x="birth",y="death",fill="") +
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  geom_abline(slope=1,intercept=0)+
  scale_fill_gradient2(limits=c( -max(abs(df0_join$z)),max(abs(df0_join$z)) ) ) +
  theme(legend.position = c(0.9, 0.25))

# Plot estimated coefficient function for dimension one
ggplot(df1_join,aes(x,y,fill=z))+
  geom_raster()+
  theme_tufte()+
  labs(x="birth",y="death",fill="") +
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  geom_abline(slope=1,intercept=0)+
  scale_fill_gradient2(limits=c( -max(abs(df1_join$z)),max(abs(df1_join$z)) ) )+
  theme(legend.position = c(0.9, 0.25))


library(spatstat)
library(tidyverse)
library(tools)
library(survival)
library(survminer)
library(ggthemes)

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

#############################################################################
# range of all persistence diagrams
#############################################################################

pdrange = function(){
  filename = list.files(path=paste0(getCurrentFileLocation(),"/persistencediagram/"),
                        pattern=".txt$",full.names = T)
  
  # recover slide id
  pd.total = NULL
  for (ii in 1:length(filename)) {
    pd=read.table(filename[ii],col.names = c("dim","birth","death"))
    ## slide id
    pd$slide=strtoi(strsplit(file_path_sans_ext(basename(filename[ii])),"_pd")[[1]][1])
    pd = pd %>% filter(!(death == Inf & dim==1)) # remove dim1 Inf death values
    pd[pd$death==Inf,3]=pd[pd$death==Inf,2] #replace Inf death with birth value
    pd.total = rbind(pd.total,pd)
  }
  pd0.total = pd.total %>% filter(dim==0)
  pd1.total = pd.total %>% filter(dim==1)
  
  # range of all persistence diagrams
  minmax0=c(floor(min(pd0.total$birth)),ceiling(max(pd0.total$death)))
  minmax1=c(floor(min(pd1.total$birth)),ceiling(max(pd1.total$death)))
  
  slideid=unique(pd.total$slide)
  
  return(list(filename=filename,minmax0=minmax0,minmax1=minmax1,
              slideid=slideid,pd0.total=pd0.total,pd1.total=pd1.total))
}

#############################################################################
# weights of persistence function
#############################################################################

# maximum distance weight (maximum of distance, absolute birth, and absolute death)
mdw=function(pd) {
  mdw=apply(cbind(pd$death-pd$birth,abs(pd$birth),abs(pd$death)),1,max)
  return(mdw)
}

# persistence weighted kernel
pwgk=function(pd,C=1,d=1) {
  pwgk=atan(C*(pd$death-pd$birth)^d)
  return(pwgk)
}

lin=function(pd) {
  lin=pd$death-pd$birth
  return(lin)
}

#############################################################################
# generate persistence function
#############################################################################

pfgen = function(sigma0,sigma1,weightfun,
                 filename,minmax0,minmax1,
                 slideid,pd0.total,pd1.total){
  
  # number of pixels in persistence function
  pixnum0 = minmax0[2]-minmax0[1]
  pixnum1 = minmax1[2]-minmax1[1]
  
  # persistence functions with the maximum distance weight
  pfres0 = as.data.frame(matrix(NA,length(filename),pixnum0*(pixnum0+1)/2))
  pfres1 = as.data.frame(matrix(NA,length(filename),pixnum1*(pixnum1+1)/2))
  
  for (ii in 1:length(filename)){
    # load aggregated persistence diagrams
    pd0=pd0.total %>% filter(slide==slideid[ii])
    pd1=pd1.total %>% filter(slide==slideid[ii])
    
    ########## dimension 0 persistence function ########## 
    dim0.x = pd0$birth
    dim0.y = pd0$death
    
    # ppp object
    pd0.ppp = ppp(dim0.x, dim0.y, minmax0, minmax0)
    
    # maximum distance weight for dimension-zero
    weight0 = weightfun(pd0)
    
    # compute persistence function
    pf0=density(pd0.ppp,sigma=sigma0,dimyx=c(pixnum0,pixnum0),
                weights=weight0)
    
    # for center of the pixel
    cunit0=(minmax0[2]-minmax0[1])/(2*pixnum0)
    # grids for birth and death
    gridx0=seq(minmax0[1]+cunit0,minmax0[2]-cunit0,length.out = pixnum0)
    gridy0=seq(minmax0[1]+cunit0,minmax0[2]-cunit0,length.out = pixnum0)
    
    # persistence function
    pf0.res=pf0$v
    df0=data.frame(x=rep(gridx0,each=pixnum0),
                   y=rep(gridy0,pixnum0),
                   z=as.vector(pf0.res))
    # some pixels may be negative, so convert them to zero
    dfres0 = df0 %>% 
      filter(y>=x) %>%
      mutate(newz=case_when(
        z < 0 ~ 0,
        TRUE ~ z)
      )
    pfres0[ii,] = dfres0$newz	
    
    ########## dimension 1 persistence function ########## 
    dim1.x = pd1$birth
    dim1.y = pd1$death
    
    # ppp object
    pd1.ppp = ppp(dim1.x, dim1.y, minmax1, minmax1)
    
    # weights
    weight1 = weightfun(pd1)
    
    # compute persistence function
    pf1=density(pd1.ppp,sigma=sigma1,dimyx=c(pixnum1,pixnum1),
                weights=weight1)
    
    cunit1=(minmax1[2]-minmax1[1])/(2*pixnum1)
    gridx1=seq(minmax1[1]+cunit1,minmax1[2]-cunit1,length.out = pixnum1)
    gridy1=seq(minmax1[1]+cunit1,minmax1[2]-cunit1,length.out = pixnum1)
    
    # persistence function
    pf1.res=pf1$v
    df1=data.frame(x=rep(gridx1,each=pixnum1),
                   y=rep(gridy1,pixnum1),
                   z=as.vector(pf1.res))
    dfres1 = df1 %>% 
      filter(y>=x) %>%
      mutate(newz=case_when(
        z < 0 ~ 0,
        TRUE ~ z)
      )
    pfres1[ii,] = dfres1$newz
  }
  
  rownames(pfres0) =rownames(pfres1) = slideid

  return(list(pfres0=pfres0,pfres1=pfres1))
}


#################################################################
# FCoxPH model estimation
#################################################################

fcoxph = function(X0,X1,slideid,clinical_info){
  hh=6 # up to 6 FPCs for each dimension
  
  # select data that have clinical information
  ## X^0 dimension zero
  X0=data.frame(X0,id=as.character(slideid))
  X0 = X0 %>% 
    semi_join(clinical_info, by=c("id"="slide_id"))	%>% 
    select(-id)

  ## X^1 dimension one
  X1=data.frame(X1,id=as.character(slideid))
  X1 = X1 %>% 
    semi_join(clinical_info, by=c("id"="slide_id"))	
  id_new = X1$id
  X1=X1 %>% 
    select(-id)
  
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
  
  comb=expand.grid(1:hh,1:hh)
  loglikelihood=AIC=rep(NA,nrow(comb))
  
  # select PCs according to AIC
  for (ss in 1:nrow(comb)) {
    # eigenfunction
    est_Eigfun0=Eigvec0[,1:comb[ss,1]]
    est_Eigfun1=Eigvec1[,1:comb[ss,2]]
    
    # eigenscore
    eigenscore0=scaleX0%*%est_Eigfun0
    eigenscore1=scaleX1%*%est_Eigfun1
    
    # designmatrix
    designmatrix=data.frame(eigenscore0,eigenscore1,id=id_new)
    colnames(designmatrix)=c(paste0("dim0.",1:comb[ss,1]),paste0("dim1.",1:comb[ss,2]),"id")
    
    # combined data
    comdf = clinical_info %>% left_join(designmatrix, by=c("slide_id"="id"))
    
    # FCoxPH model
    ## formula for FCoxPH model
    inifor="Surv(survival_time_new, dead)~age+tobacco+female+stage+tumorsize+grade+grade"
    for (s in 1:comb[ss,1]){
      inifor=paste0(inifor,paste0("+dim0.",s))
    }
    for (s in 1:comb[ss,2]){
      inifor=paste0(inifor,paste0("+dim1.",s))
    }
    inifor=paste0(inifor,paste0("+cluster(patient_id)"))
    
    ## FCoxPH model fitting
    model.cox = coxph(as.formula(inifor),data=comdf)
    finalresult=summary(model.cox)
    loglikelihood[ss]=finalresult$loglik[2]
    AIC[ss]=2*(comb[ss,1]+comb[ss,2])-2*loglikelihood[ss]
  }
  
  # selected number of FPCs by AIC
  AICselectindex=which.min(AIC)
  
  # eigenfunction
  est_Eigfun0=Eigvec0[,1:comb[AICselectindex,1]]
  est_Eigfun1=Eigvec1[,1:comb[AICselectindex,2]]
  
  # eigenscore
  eigenscore0=scaleX0%*%est_Eigfun0
  eigenscore1=scaleX1%*%est_Eigfun1
  
  designmatrix=data.frame(eigenscore0,eigenscore1,id=id_new)
  colnames(designmatrix)=c(paste0("dim0.",1:comb[AICselectindex,1]),paste0("dim1.",1:comb[AICselectindex,2]),"id")
  
  # combined data
  comdf = clinical_info %>% left_join(designmatrix, by=c("slide_id"="id"))
  
  # selected FCoxPH model by AIC
  ## formula for selected FCoxPH model
  inifor="Surv(survival_time_new, dead)~age+tobacco+female+stage+tumorsize+grade"
  for (s in 1:comb[AICselectindex,1]){
    inifor=paste0(inifor,paste0("+dim0.",s))
  }
  for (s in 1:comb[AICselectindex,2]){
    inifor=paste0(inifor,paste0("+dim1.",s))
  }
  inifor=paste0(inifor,paste0("+cluster(patient_id)"))
  
  ## FCoxPH model fitting
  model.fcoxph = coxph(as.formula(inifor),data=comdf)
  model.coxph = coxph(Surv(survival_time_new, dead)~age+tobacco+female+stage+tumorsize+grade+cluster(patient_id),data=comdf)
  return(list(model.fcoxph=model.fcoxph,model.coxph=model.coxph,
              Eigvec0=Eigvec0,Eigvec1=Eigvec1,id=id_new,
              est_Eigval0=est_Eigval0,est_Eigval1=est_Eigval1,
              AICselectindex=AICselectindex))
}

fcoef = function(model.fcoxph,Eigvec0,Eigvec1,minmax0,minmax1,AICselectindex){
  hh=6 # up to 6 FPCs for each dimension
  comb=expand.grid(1:hh,1:hh)
  num.clinical = 11 # number of clinical variables
  beta.est0 = as.matrix(Eigvec0[,1:comb[AICselectindex,1]]) %*% as.matrix(model.fcoxph$coefficients[(num.clinical+1):(num.clinical+comb[AICselectindex,1])] )
  beta.est1 = as.matrix(Eigvec1[,1:comb[AICselectindex,2]]) %*% as.matrix(model.fcoxph$coefficients[(num.clinical+comb[AICselectindex,1]+1):(num.clinical+comb[AICselectindex,1]+comb[AICselectindex,2])] )
  
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
  
  # data for dimension one
  df1=data.frame(x=rep(gridx1,each=res1),y=rep(gridy1,res1))  %>%
    filter(y>=x)
  df1_part = data.frame(x=rep(gridx1,each=res1),y=rep(gridy1,res1)) %>%
    filter(y>=x)
  df1_part$z=beta.est1[1:nrow(df1_part),]
  df1_join = left_join(df1,df1_part,by=c("x"="x","y"="y"))
  df1_join$z = ifelse(is.na(df1_join$z),0,df1_join$z)
  return(list(df0_join,df1_join))
}


validtest = function(threshold=0.9,X0,X1,slideid,clinical_info){

  # select data that have clinical information
  ## X^0 dimension zero
  X0=data.frame(X0,id=as.character(slideid))
  X0 = X0 %>% 
    semi_join(clinical_info, by=c("id"="slide_id"))	%>% 
    select(-id)
  
  ## X^1 dimension one
  X1=data.frame(X1,id=as.character(slideid))
  X1 = X1 %>% 
    semi_join(clinical_info, by=c("id"="slide_id"))	
  id_new = X1$id
  X1=X1 %>% 
    select(-id)
  
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
  
  # select dimension zero FPCs according to the variability threshold
  for (ss in 1:length(est_Eigval0)) {
    if (sum(est_Eigval0[1:ss])/sum(est_Eigval0)>threshold)
      break
  }
  selectindex0=ss
  est_Eigfun0_var=Eigvec0[,1:selectindex0]
  eigenscore0_var=data.frame(scaleX0%*%est_Eigfun0_var)
  
  # select dimension one FPCs according to the variability threshold
  for (ss in 1:length(est_Eigval1)) {
    if (sum(est_Eigval1[1:ss])/sum(est_Eigval1)>threshold)
      break
  }
  selectindex1=ss
  
  est_Eigfun0_var=Eigvec0[,1:selectindex0]
  est_Eigfun1_var=Eigvec1[,1:selectindex1]
  
  eigenscore0_var=scaleX0%*%est_Eigfun0_var
  eigenscore1_var=scaleX1%*%est_Eigfun1_var
  
  designmatrix_var=data.frame(eigenscore0_var,eigenscore1_var,id=id_new)
  colnames(designmatrix_var)=c(paste0("dim0.",1:selectindex0),paste0("dim1.",1:selectindex1),"id")
  
  # combined data
  comdf_var = clinical_info %>% left_join(designmatrix_var, by=c("slide_id"="id"))
  
  # FCoxPH model variable selected by variance threshold
  ## formula
  inifor="Surv(survival_time_new, dead)~age+tobacco+female+stage+tumorsize+grade"
  for (s in 1:selectindex0){
    inifor=paste0(inifor,paste0("+dim0.",s))
  }
  for (s in 1:selectindex1){
    inifor=paste0(inifor,paste0("+dim1.",s))
  }
  inifor=paste0(inifor,paste0("+cluster(patient_id)"))
  
  ## model fitting
  model.fcoxph_var = coxph(as.formula(inifor),data=comdf_var)
  finalresult=summary(model.fcoxph_var)
  
  ## initial value for chi-sq test
  model.coxph = coxph(Surv(survival_time_new, dead)~age+tobacco+female+stage+tumorsize+grade+cluster(patient_id),data=clinical_info)
  initial.value=c(model.coxph$coef,rep(0,selectindex0+selectindex1 ) )
  
  ## coefficient for full model
  testscore=model.fcoxph_var$score
  
  # chi-sq test
  criticalvalue=1-pchisq(testscore, selectindex0+selectindex1, lower.tail = TRUE, log.p = FALSE)
  
  return(list(selectindex0,selectindex1,criticalvalue))
}

##################################################################
# LOOCV risk prediction for FCoxPH model
##################################################################
cvpred = function(X0,X1,slideid,clinical_info){
  # select data that have clinical information
  ## X^0 dimension zero
  X0=data.frame(X0,id=as.character(slideid))
  X0 = X0 %>% 
    semi_join(clinical_info, by=c("id"="slide_id"))	%>% 
    select(-id)
  
  ## X^1 dimension one
  X1=data.frame(X1,id=as.character(slideid))
  X1 = X1 %>% 
    semi_join(clinical_info, by=c("id"="slide_id"))	
  id_new = X1$id
  X1=X1 %>% 
    select(-id)
  
  # predicted risk of FCoxPH for LOOCV
  risk.pred.fcoxph = rep(NA,nrow(X0))
  # predicted risk of CoxPH for LOOCV
  risk.pred.coxph = rep(NA,nrow(X0))
  
  hh=6
  
  # LOOCV
  for (ii in 1:nrow(X0)){
    print(ii)
    rownames(X0) = id_new
    rownames(X1) = id_new
    X0.train = X0[-ii,]
    X0.test = X0[ii,]
    X1.train = X1[-ii,]
    X1.test = X1[ii,]
    id.train = rownames(X0.train)
    id.test = rownames(X0.test)
    
    # spectral decomposition of X^0
    ## X^0-mu^0 = U*Sigma*V^T
    scaleX0.train=scale(X0.train,center=TRUE,scale=FALSE)
    ## eigenfunction of X^0-mu^0. phi_j
    Eigvec0=svd(scaleX0.train)$v

    # spectral decomposition of X^1
    ## X^1-mu^1 = U*Sigma*V^T
    scaleX1.train=scale(X1.train,center=TRUE,scale=FALSE)
    ## eigenfunction of X^1-mu^1. pi_k
    Eigvec1=svd(scaleX1.train)$v

    
    comb=expand.grid(1:hh,1:hh)
    loglikelihood=AIC=rep(NA,nrow(comb))
    
    # select PCs according to AIC
    for (ss in 1:nrow(comb)) {
      # eigenfunction
      est_Eigfun0=Eigvec0[,1:comb[ss,1]]
      est_Eigfun1=Eigvec1[,1:comb[ss,2]]
      
      # eigenscore
      eigenscore0=scaleX0.train%*%est_Eigfun0
      eigenscore1=scaleX1.train%*%est_Eigfun1
      
      # designmatrix
      designmatrix.train=data.frame(eigenscore0,eigenscore1,id=id.train)
      colnames(designmatrix.train)=c(paste0("dim0.",1:comb[ss,1]),paste0("dim1.",1:comb[ss,2]),"id")
      
      # combined data
      comdf.train = designmatrix.train %>% left_join(clinical_info, by=c("id"="slide_id"))
      
      # FCoxPH model
      ## formula for FCoxPH model
      inifor="Surv(survival_time_new, dead)~age+tobacco+female+stage+tumorsize+grade"
      for (s in 1:comb[ss,1]){
        inifor=paste0(inifor,paste0("+dim0.",s))
      }
      for (s in 1:comb[ss,2]){
        inifor=paste0(inifor,paste0("+dim1.",s))
      }
      inifor=paste0(inifor,paste0("+cluster(patient_id)"))
      
      ## FCoxPH model fitting
      model.cox = coxph(as.formula(inifor),data=comdf.train)
      finalresult=summary(model.cox)
      loglikelihood[ss]=finalresult$loglik[2]
      AIC[ss]=2*(comb[ss,1]+comb[ss,2])-2*loglikelihood[ss]
    }
    
    # selected number of FPCs by AIC
    AICselectindex=which.min(AIC)
    
    # eigenfunction
    est_Eigfun0=Eigvec0[,1:comb[AICselectindex,1]]
    est_Eigfun1=Eigvec1[,1:comb[AICselectindex,2]]
    
    # eigenscore
    eigenscore0=scaleX0.train%*%est_Eigfun0
    eigenscore1=scaleX1.train%*%est_Eigfun1
    
    designmatrix.train=data.frame(eigenscore0,eigenscore1,id=id.train)
    colnames(designmatrix.train)=c(paste0("dim0.",1:comb[AICselectindex,1]),paste0("dim1.",1:comb[AICselectindex,2]),"id")
    
    # combined data
    comdf.train = designmatrix.train %>% left_join(clinical_info, by=c("id"="slide_id"))
    
    # FCoxPH model by AIC
    ## formula for FCoxPH model
    inifor="Surv(survival_time_new, dead)~age+tobacco+female+stage+tumorsize+grade"
    for (s in 1:comb[AICselectindex,1]){
      inifor=paste0(inifor,paste0("+dim0.",s))
    }
    for (s in 1:comb[AICselectindex,2]){
      inifor=paste0(inifor,paste0("+dim1.",s))
    }
    inifor=paste0(inifor,paste0("+cluster(patient_id)"))
    
    ## FCoxPH model fitting
    model.fcoxph = coxph(as.formula(inifor),data=comdf.train)
    
    # test set 
    ## substract mean of training data
    scaleX0.test = X0.test-colMeans(X0.train)
    scaleX1.test = X1.test-colMeans(X1.train)
    ## estimated eigenscore of test data
    eigenscore0.test=as.matrix(scaleX0.test)%*%est_Eigfun0
    eigenscore1.test=as.matrix(scaleX1.test)%*%est_Eigfun1
    ## design matrix of test data
    designmatrix.test=data.frame(eigenscore0.test,eigenscore1.test,id=id.test)
    colnames(designmatrix.test)=c(paste0("dim0.",1:comb[AICselectindex,1]),paste0("dim1.",1:comb[AICselectindex,2]),"id")
    # combined test data
    comdf.test = designmatrix.test %>% left_join(clinical_info, by=c("id"="slide_id"))
    ## risk prediction
    risk.pred.fcoxph[ii]=predict(model.fcoxph,comdf.test,type="risk")
    
    ## CoxPH model fitting
    ## formula for CoxPH model
    inifor.coxph="Surv(survival_time_new, dead)~age+tobacco+female+stage+tumorsize+grade+cluster(patient_id)"
    model.coxph = coxph(as.formula(inifor.coxph),data=comdf.train)
    ## risk prediction
    risk.pred.coxph[ii]=predict(model.coxph,comdf.test,type="risk")
  }
  
  ## risk data
  ### FCoxPH
  risk.df.fcoxph=data.frame(id=rownames(X0),risk=risk.pred.fcoxph)
  risk.mean.fcoxph = clinical_info %>% 
    left_join(risk.df.fcoxph, by=c("slide_id"="id")) %>% 
    group_by(patient_id) %>% # mean predicted risk per patient
    summarize(meanrisk=mean(risk),
              survival_time_new=mean(survival_time_new),
              dead=mean(dead))
  risk.mean.fcoxph$group=ifelse(risk.mean.fcoxph$meanrisk<=median(risk.mean.fcoxph$meanrisk),"low","high")
  
  ### CoxPH
  risk.df.coxph=data.frame(id=rownames(X0),risk=risk.pred.coxph)
  risk.mean.coxph = clinical_info %>% 
    left_join(risk.df.coxph, by=c("slide_id"="id")) %>% 
    group_by(patient_id) %>% # mean predicted risk per patient
    summarize(meanrisk=mean(risk),
              survival_time_new=mean(survival_time_new),
              dead=mean(dead))
  risk.mean.coxph$group=ifelse(risk.mean.coxph$meanrisk<=median(risk.mean.coxph$meanrisk),"low","high")
  
  return(list(risk.pred.fcoxph=risk.pred.fcoxph, risk.pred.coxph=risk.pred.coxph,
              risk.mean.fcoxph=risk.mean.fcoxph, risk.mean.coxph=risk.mean.coxph))
}
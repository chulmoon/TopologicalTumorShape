main=function(iter){
  # smoothing parameters
  sigma0=0.3
  sigma1=0.5
  
  # range of persistence functions
  minmax0=c(-22,26)
  minmax1=c(-6,32)
  pixnum0=48
  pixnum1=38
  
  filename = list.files(path="./TopologicalTumorShape/data/brain/PersistenceDiagram-rearranged/"
                        ,pattern=paste0("_",iter-1,"_pd.txt$"),full.names = T)
  
  # aggregated pd
  pd.total=NULL
  for (ii in 1:length(filename)) {
    print(ii)
    pd=read.table(filename[ii],col.names = c("dim","birth","death"))
    
    strsplit(strsplit(file_path_sans_ext(basename(filename[ii])),"-")[[1]][3],"_")
    
    id1=strsplit(file_path_sans_ext(basename(filename[ii])),"-")[[1]][2]
    id2=strsplit(strsplit(file_path_sans_ext(basename(filename[ii])),"-")[[1]][3],"_")[[1]][[1]]
    slideid=strsplit(strsplit(file_path_sans_ext(basename(filename[ii])),"-")[[1]][3],"_")[[1]][[2]]
    
    pd$id1 = id1
    pd$id2 = id2
    pd$slideid = slideid
    
    pd.total = rbind(pd.total,pd)
  }
  
  
  # replace Inf death values with its birth values
  pd.total$death[pd.total$death==Inf]=pd.total$birth[pd.total$death==Inf]
  
  # rescale the birth and death values
  pd.total$death = pd.total$death/2
  pd.total$birth = pd.total$birth/2
  
  pd0.total = pd.total %>% 
    filter(dim==0)
  pd1.total = pd.total %>% 
    filter(dim==1)
  
  # patient id
  patientid = pd.total %>%
    select(id1,id2) %>%
    distinct()
  
  # maximum distance weight (maximum of distance, absolute birth, and absolute death)
  mdw=function(pd) {
    mdw=apply(cbind(pd$death-pd$birth,abs(pd$birth),abs(pd$death)),1,max)
    return(mdw)
  }
  
  ##############################################################
  # generate persistence surface function
  ##############################################################

  # persistence functions with the maximum distance weight
  pfmdw0 = as.data.frame(matrix(NA,nrow(patientid),pixnum0*(pixnum0+1)/2))
  pfmdw1 =	as.data.frame(matrix(NA,nrow(patientid),pixnum1*(pixnum1+1)/2))
  ntumor  =rep(NA,nrow(patientid))
  
  for (ii in 1:nrow(patientid)){
    pd0 = pd0.total %>% 
      filter(id1 == patientid$id1[ii], id2 == patientid$id2[ii])
    pd1 = pd1.total %>% 
      filter(id1 == patientid$id1[ii], id2 == patientid$id2[ii])
    sizedat = idslice100mat %>% 
      filter(id1 == patientid$id1[ii], id2 == patientid$id2[ii])
    
    ntumor[ii] = mean(sizedat$ntumor)
    nslide=length(unique(pd0$slideid))
    
    ### dimension 0
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
    pfmdw0[ii,] = dfmw0$newz/nslide # average density function	
    
    #################### dimension 1
    dim1.x = pd1$birth
    dim1.y = pd1$death
    nslide=length(unique(pd1$slideid))
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
    pfmdw1[ii,] = dfmdw1$newz/nslide # average density function	
  }
  
  #################################################################
  # FCoxPH model estimation
  #################################################################
  # persistence function with the maximum distance weight
  rownames(pfmdw0) = ids
  rownames(pfmdw1) = ids
  X0=pfmdw0
  X1=pfmdw1
  
  # select non overlapping data
  X0 = X0[which(rownames(X0)%in%clinical.info.brain$patid),]
  X1 = X1[which(rownames(X1)%in%clinical.info.brain$patid),]
  
  ntumor=data.frame(ntumor)
  rownames(ntumor) = ids
  ntumor = ntumor[which(rownames(ntumor)%in%clinical.info.brain$patid),]
  ntumor = as.vector(ntumor)
  
  # same id for X0 and X1
  id_new = rownames(X0)
  
  
  # predicted risk of FCoxPH for LOOCV
  risk.pred.fcoxph = rep(NA,nrow(X0))
  # predicted risk of CoxPH for LOOCV
  risk.pred.coxph = rep(NA,nrow(X0))
  
  # LOOCV
  hh=4
  for (ii in 1:nrow(X0)){
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
    ## eigenvalue of X^0-mu^0
    #Singval0=svd(scaleX0)$d
    ## eigenvalue of (X^0-mu^0)'(X^0-mu^0)
    #est_Eigval0=Singval0^2
    
    # spectral decomposition of X^1
    ## X^1-mu^1 = U*Sigma*V^T
    scaleX1.train=scale(X1.train,center=TRUE,scale=FALSE)
    ## eigenfunction of X^1-mu^1. pi_k
    Eigvec1=svd(scaleX1.train)$v
    ## eigenvalue of X^1-mu^1
    #Singval1=svd(scaleX1)$d
    ## eigenvalue of (X^1-mu^1)'(X^1-mu^1)
    #est_Eigval1=Singval1^2
    
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
      designmatrix.train=data.frame(eigenscore0,eigenscore1,ntumor[-ii],id=id.train)
      colnames(designmatrix.train)=c(paste0("dim0.",1:comb[ss,1]),paste0("dim1.",1:comb[ss,2]),"ntumor", "id")
      
      # combined data
      comdf.train = designmatrix.train %>% 
        left_join(clinical.info.brain, c("id"="patid"))
      
      # FCoxPH model
      ## formula for FCoxPH model
      inifor="Surv(os,event)~age+ks+gender+ntumor"
      for (s in 1:comb[ss,1]){
        inifor=paste0(inifor,paste0("+dim0.",s))
      }
      for (s in 1:comb[ss,2]){
        inifor=paste0(inifor,paste0("+dim1.",s))
      }
      
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
    
    designmatrix.train=data.frame(eigenscore0,eigenscore1,ntumor[-ii],id=id.train)
    colnames(designmatrix.train)=c(paste0("dim0.",1:comb[AICselectindex,1]),paste0("dim1.",1:comb[AICselectindex,2]),"ntumor", "id")
    
    # combined data
    comdf.train = designmatrix.train %>% 
      left_join(clinical.info.brain, c("id"="patid"))
    
    # FCoxPH model by AIC
    ## formula for FCoxPH model
    inifor="Surv(os,event)~age+ks+gender+ntumor"
    for (s in 1:comb[AICselectindex,1]){
      inifor=paste0(inifor,paste0("+dim0.",s))
    }
    for (s in 1:comb[AICselectindex,2]){
      inifor=paste0(inifor,paste0("+dim1.",s))
    }
    ## FCoxPH model fitting
    model.fcoxph = coxph(as.formula(inifor),data=comdf.train)
    
    # test data 
    ## subtract mean of training data
    scaleX0.test = X0.test-colMeans(X0.train)
    scaleX1.test = X1.test-colMeans(X1.train)
    ## estimated eigenscore of test data
    eigenscore0.test=as.matrix(scaleX0.test)%*%est_Eigfun0
    eigenscore1.test=as.matrix(scaleX1.test)%*%est_Eigfun1
    ## design matrix of test data
    designmatrix.test=data.frame(eigenscore0.test,eigenscore1.test,ntumor[ii],id=id.test)
    colnames(designmatrix.test)=c(paste0("dim0.",1:comb[AICselectindex,1]),paste0("dim1.",1:comb[AICselectindex,2]),"ntumor","id")
    # combined test data
    comdf.test = designmatrix.test %>% 
      left_join(clinical.info.brain, c("id"="patid"))
    ## risk prediction
    risk.pred.fcoxph[ii]=predict(model.fcoxph,comdf.test,type="risk")
    
  }
  
  ##################################################################
  # Kaplan-Meier plots
  ##################################################################
  # FcoxPH
  ## risk data
  risk.df.fcoxph=data.frame(id=rownames(X0),risk=risk.pred.fcoxph)
  risk.mean.fcoxph = data.frame(id=id_new) %>% 
    left_join(clinical.info.brain, by=c("id"="patid")) %>%
    left_join(risk.df.fcoxph, by=c("id"))
  risk.mean.fcoxph$group=ifelse(risk.mean.fcoxph$risk<=median(risk.mean.fcoxph$risk),"low","high")
  
  risk.result = survival::survdiff(Surv(os, event) ~ group, data = risk.mean.fcoxph)
  rescv = pchisq(risk.result$chisq, length(risk.result$n)-1, lower.tail = FALSE)
  return(rescv)
}

library(snow)

# simulation setting
niter = 50
permres = rep(NA,niter)

# number of cores for parallel computation
numCores = parallel::detectCores()
cl=snow::makeCluster(numCores-1, type="SOCK")
clusterEvalQ(cl,{library(spatstat);library(tidyverse);library(tools);library(survival)})
permres=snow::parApply(cl,matrix(1:niter,niter,1),1,main)
stopCluster(cl)

# save p-values
save(permres, file="./TopologicalTumorShape/data/brain/brain_loocv_result.Rdata")

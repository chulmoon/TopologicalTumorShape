main=function(iter){
  # smoothing parameters
  sigma0=2
  sigma1=0.9
  
  # range of persistence functions
  minmax0=c(-41,11)
  minmax1=c(-19,26)
  pixnum0=52
  pixnum1=45
  
  filename = list.files(path="./TopologicalTumorShape/data/lung/PersistenceDiagram-rearranged/"
                        ,pattern=paste0("_",iter-1,"_pd.txt$"),full.names = T)
  # recover slide id
  pd.total = NULL
  for (ii in 1:length(filename)) {
    pd=read.table(filename[ii],col.names = c("dim","birth","death"))
    ## slide id
    pd$slide=strtoi(strsplit(file_path_sans_ext(basename(filename[ii])),paste0("_",iter-1,"_pd"))[[1]][1])
    pd = pd %>% filter(!(death == Inf & dim==1)) # remove dim1 Inf death values
    pd[pd$death==Inf,3]=pd[pd$death==Inf,2] #replace Inf death with birth value
    pd.total = rbind(pd.total,pd)
  }
  pd0.total = pd.total %>% filter(dim==0)
  pd1.total = pd.total %>% filter(dim==1)
  
  slideid=unique(pd.total$slide)
  
  #############################################################################
  # generate persistence functions
  #############################################################################
  
  # number of topological features in each slide
  pd0.int = pd0.total %>%
    group_by(slide) %>%
    summarise(nint0=n())
  pd1.int = pd1.total %>%
    group_by(slide) %>%
    summarise(nint1=n())
  
  # persistence functions with the maximum distance weight
  pfmdw0 = as.data.frame(matrix(NA,length(filename),pixnum0*(pixnum0+1)/2))
  pfmdw1 = as.data.frame(matrix(NA,length(filename),pixnum1*(pixnum1+1)/2))
  
  
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
  
  # persistence function with linear weight
  rownames(pfmdw0) = slideid
  rownames(pfmdw1) = slideid
  X0=pfmdw0
  X1=pfmdw1
  
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
  
  ##################################################################
  # LOOCV risk prediction for FCoxPH model
  ##################################################################
  
  # predicted risk of FCoxPH for LOOCV
  risk.pred.fcoxph = rep(NA,nrow(X0))
  # predicted risk of CoxPH for LOOCV
  risk.pred.coxph = rep(NA,nrow(X0))
  
  hh=5
  
  # LOOCV
  for (ii in 1:nrow(X0)){
    #print(ii)
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
      designmatrix.train=data.frame(eigenscore0,eigenscore1,id=id.train)
      colnames(designmatrix.train)=c(paste0("dim0.",1:comb[ss,1]),paste0("dim1.",1:comb[ss,2]),"id")
      
      # combined data
      comdf.train = designmatrix.train %>% left_join(clinical_info, by=c("id"="slide_id"))
      
      # FCoxPH model
      ## formula for FCoxPH model
      inifor="Surv(survival_time_new, dead)~age+tobacco+female+stage+tumorsize"
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
    inifor="Surv(survival_time_new, dead)~age+tobacco+female+stage+tumorsize"
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
    ## subtract mean of training data
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
    inifor.coxph="Surv(survival_time_new, dead)~age+tobacco+female+stage+tumorsize+cluster(patient_id)"
    model.coxph = coxph(as.formula(inifor.coxph),data=comdf.train)
    ## risk prediction
    risk.pred.coxph[ii]=predict(model.coxph,comdf.test,type="risk")
  }
  
  # FcoxPH
  ## risk data
  risk.df.fcoxph=data.frame(id=rownames(X0),risk=risk.pred.fcoxph)
  risk.mean.fcoxph = clinical_info %>% 
    left_join(risk.df.fcoxph, by=c("slide_id"="id")) %>% 
    group_by(patient_id) %>% # mean predicted risk per patient
    summarize(meanrisk=mean(risk),
              survival_time_new=mean(survival_time_new),
              dead=mean(dead))
  risk.mean.fcoxph$group=ifelse(risk.mean.fcoxph$meanrisk<=median(risk.mean.fcoxph$meanrisk),"low","high")
  risk.result=survdiff(Surv(survival_time_new, dead) ~ group, data = risk.mean.fcoxph)
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
save(permres, file="./TopologicalTumorShape/data/lung/lung_loocv_result.Rdata")

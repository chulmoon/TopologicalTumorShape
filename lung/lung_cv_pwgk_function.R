main=function(cc, clinical_info){

  sig0cv = sig1cv = seq(0.2,5,by=0.2)
  sigcv = expand.grid(sig0cv,sig1cv)
  
  # sigma for smoothing function
  sigma0 = sigcv[cc,1]
  sigma1 = sigcv[cc,2]
  
  # persistence functions with the maximum distance weight
  pfpwgk0 = as.data.frame(matrix(NA,length(filename),pixnum0*(pixnum0+1)/2))
  pfpwgk1 = as.data.frame(matrix(NA,length(filename),pixnum1*(pixnum1+1)/2))
  

  pwgk=function(pd,C=1,d=1) {
    pwgk=atan(C*(pd$death-pd$birth)^d)
    return(pwgk)
  }

  
  for (ii in 1:length(filename)){
    # load aggregated persistence diagrams
    pd0=pd0.total %>% filter(slide==slideid[ii])
    pd1=pd1.total %>% filter(slide==slideid[ii])
    
    ########## dimension 0 persistence function ########## 
    dim0.x = pd0$birth
    dim0.y = pd0$death
    
    # ppp object
    pd0.ppp = ppp(dim0.x, dim0.y, minmax0, minmax0)
    
    # pwgk for dimension-zero
    pwgk0 = pwgk(pd0)
    
    # compute persistence function
    pf0=density.ppp(pd0.ppp,sigma=sigma0,dimyx=c(pixnum0,pixnum0),
                weights=pwgk0)
    
    # for center of the pixel
    cunit0=(minmax0[2]-minmax0[1])/(2*pixnum0)
    # grids for birth and death
    gridx0=seq(minmax0[1]+cunit0,minmax0[2]-cunit0,length.out = pixnum0)
    gridy0=seq(minmax0[1]+cunit0,minmax0[2]-cunit0,length.out = pixnum0)
    
    # persistence function
    pf0.pwgk=pf0$v
    df0=data.frame(x=rep(gridx0,each=pixnum0),
                   y=rep(gridy0,pixnum0),
                   z=as.vector(pf0.pwgk))
    # some pixels are negative, so convert them to zero
    dfmw0 = df0 %>% 
      filter(y>=x) %>%
      mutate(newz=case_when(
        z < 0 ~ 0,
        TRUE ~ z)
      )
    pfpwgk0[ii,] = dfmw0$newz # persistence function	
    
    ########## dimension 1 persistence function ########## 
    dim1.x = pd1$birth
    dim1.y = pd1$death
    
    # ppp object
    pd1.ppp = ppp(dim1.x, dim1.y, minmax1, minmax1)
    
    # weights
    pwgk1 = pwgk(pd1)
    
    # compute persistence function
    pf1=density(pd1.ppp,sigma=sigma1,dimyx=c(pixnum1,pixnum1),
                weights=pwgk1)
    
    cunit1=(minmax1[2]-minmax1[1])/(2*pixnum1)
    gridx1=seq(minmax1[1]+cunit1,minmax1[2]-cunit1,length.out = pixnum1)
    gridy1=seq(minmax1[1]+cunit1,minmax1[2]-cunit1,length.out = pixnum1)
    
    # maximum distance weight
    pf1.pwgk=pf1$v
    df1=data.frame(x=rep(gridx1,each=pixnum1),
                   y=rep(gridy1,pixnum1),
                   z=as.vector(pf1.pwgk))
    dfpwgk1 = df1 %>% 
      filter(y>=x) %>%
      mutate(newz=case_when(
        z < 0 ~ 0,
        TRUE ~ z)
      )
    pfpwgk1[ii,] = dfpwgk1$newz # persistence function
  }
  
  
  # persistence function with the maximum distance weight
  rownames(pfpwgk0) = slideid
  rownames(pfpwgk1) = slideid
  X0=pfpwgk0
  X1=pfpwgk1
  
  
  hh=5
  # select data that have clinical information
  ## X^0 dimension zero
  X0=data.frame(X0,id=as.character(slideid))
  X0 = X0 %>% 
    dplyr::semi_join(clinical_info, by=c("id"="slide_id"))	%>% 
    dplyr::select(-id)
  
  ## X^1 dimension one
  X1=data.frame(X1,id=as.character(slideid))
  X1 = X1 %>% 
    dplyr::semi_join(clinical_info, by=c("id"="slide_id"))	
  id_new = X1$id
  X1=X1 %>% 
    dplyr::select(-id)
  
  # predicted risk of FCoxPH for LOOCV
  risk.pred.fcoxph = rep(NA,nrow(X0))
  # predicted risk of CoxPH for LOOCV
  risk.pred.coxph = rep(NA,nrow(X0))
  
  # LOOCV
  for (ii in 1:nrow(X0)){
    
    
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
        eigenscore1=scaleX1.train%*%est_Eigfun1
        
        # designmatrix
        designmatrix.train=data.frame(eigenscore1,id=id.train)
        colnames(designmatrix.train)=c(paste0("dim1.",1:comb[ss,2]),"id")
        
        # combined data
        comdf.train = designmatrix.train %>% dplyr::left_join(clinical_info, by=c("id"="slide_id"))
        
        # FCoxPH model
        ## formula for FCoxPH model
        inifor="Surv(survival_time_new, dead)~age+tobacco+female+stage+tumorsize+grade"
        for (s in 1:comb[ss,2]){
          inifor=paste0(inifor,paste0("+dim1.",s))
        }
        inifor=paste0(inifor,paste0("+cluster(patient_id)"))
        
        ## FCoxPH model fitting
        model.cox = coxph(as.formula(inifor),data=comdf.train)
        finalresult=summary(model.cox)
        loglikelihood[ss]=finalresult$loglik[2]
        AIC[ss]=2*(comb[ss,1]+comb[ss,2])-2*loglikelihood[ss]
        
      } else if (comb[ss,2]==0){
        
        # eigenfunction
        est_Eigfun0=Eigvec0[,1:comb[ss,1]]
        
        # eigenscore
        eigenscore0=scaleX0.train%*%est_Eigfun0
        
        # designmatrix
        designmatrix.train=data.frame(eigenscore0,id=id.train)
        colnames(designmatrix.train)=c(paste0("dim0.",1:comb[ss,1]),"id")
        
        # combined data
        comdf.train = designmatrix.train %>% dplyr::left_join(clinical_info, by=c("id"="slide_id"))
        
        # FCoxPH model
        ## formula for FCoxPH model
        inifor="Surv(survival_time_new, dead)~age+tobacco+female+stage+tumorsize+grade"
        for (s in 1:comb[ss,1]){
          inifor=paste0(inifor,paste0("+dim0.",s))
        }
        inifor=paste0(inifor,paste0("+cluster(patient_id)"))
        
        ## FCoxPH model fitting
        model.cox = coxph(as.formula(inifor),data=comdf.train)
        finalresult=summary(model.cox)
        loglikelihood[ss]=finalresult$loglik[2]
        AIC[ss]=2*(comb[ss,1]+comb[ss,2])-2*loglikelihood[ss]
        
      } else{
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
        comdf.train = designmatrix.train %>% dplyr::left_join(clinical_info, by=c("id"="slide_id"))
        
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
    }
    
    # selected number of FPCs by AIC
    AICselectindex=which.min(AIC)
    
    if(comb[AICselectindex,1]==0 & comb[AICselectindex,2]==0){

      designmatrix.train=data.frame(id=id.train)
      
      # combined data
      comdf.train = designmatrix.train %>% dplyr::left_join(clinical_info, by=c("id"="slide_id"))
      
      # FCoxPH model by AIC
      ## formula for FCoxPH model
      inifor="Surv(survival_time_new, dead)~age+tobacco+female+stage+tumorsize+grade"
      inifor=paste0(inifor,paste0("+cluster(patient_id)"))
      
      ## FCoxPH model fitting
      model.fcoxph = coxph(as.formula(inifor),data=comdf.train)
      
      # test set 
      designmatrix.test=data.frame(id=id.test)
      # combined test data
      comdf.test = designmatrix.test %>% dplyr::left_join(clinical_info, by=c("id"="slide_id"))
      ## risk prediction
      risk.pred.fcoxph[ii]=predict(model.fcoxph,comdf.test,type="risk")
      
    } else if (comb[AICselectindex,1]==0){
      
      # eigenfunction
      est_Eigfun1=Eigvec1[,1:comb[AICselectindex,2]]
      
      # eigenscore
      eigenscore1=scaleX1.train%*%est_Eigfun1
      
      designmatrix.train=data.frame(eigenscore1,id=id.train)
      colnames(designmatrix.train)=c(paste0("dim1.",1:comb[AICselectindex,2]),"id")
      
      # combined data
      comdf.train = designmatrix.train %>% dplyr::left_join(clinical_info, by=c("id"="slide_id"))
      
      # FCoxPH model by AIC
      ## formula for FCoxPH model
      inifor="Surv(survival_time_new, dead)~age+tobacco+female+stage+tumorsize+grade"
      for (s in 1:comb[AICselectindex,2]){
        inifor=paste0(inifor,paste0("+dim1.",s))
      }
      inifor=paste0(inifor,paste0("+cluster(patient_id)"))
      
      ## FCoxPH model fitting
      model.fcoxph = coxph(as.formula(inifor),data=comdf.train)
      
      # test set 
      ## substract mean of training data
      scaleX1.test = X1.test-colMeans(X1.train)
      ## estimated eigenscore of test data
      eigenscore1.test=as.matrix(scaleX1.test)%*%est_Eigfun1
      ## design matrix of test data
      designmatrix.test=data.frame(eigenscore1.test,id=id.test)
      colnames(designmatrix.test)=c(paste0("dim1.",1:comb[AICselectindex,2]),"id")
      # combined test data
      comdf.test = designmatrix.test %>% dplyr::left_join(clinical_info, by=c("id"="slide_id"))
      ## risk prediction
      risk.pred.fcoxph[ii]=predict(model.fcoxph,comdf.test,type="risk")
      
    } else if (comb[AICselectindex,2]==0){
      
      # eigenfunction
      est_Eigfun0=Eigvec0[,1:comb[AICselectindex,1]]
      
      # eigenscore
      eigenscore0=scaleX0.train%*%est_Eigfun0
      
      designmatrix.train=data.frame(eigenscore0,id=id.train)
      colnames(designmatrix.train)=c(paste0("dim0.",1:comb[AICselectindex,1]),"id")
      
      # combined data
      comdf.train = designmatrix.train %>% dplyr::left_join(clinical_info, by=c("id"="slide_id"))
      
      # FCoxPH model by AIC
      ## formula for FCoxPH model
      inifor="Surv(survival_time_new, dead)~age+tobacco+female+stage+tumorsize+grade"
      for (s in 1:comb[AICselectindex,1]){
        inifor=paste0(inifor,paste0("+dim0.",s))
      }
      inifor=paste0(inifor,paste0("+cluster(patient_id)"))
      
      ## FCoxPH model fitting
      model.fcoxph = coxph(as.formula(inifor),data=comdf.train)
      
      # test set 
      ## substract mean of training data
      scaleX0.test = X0.test-colMeans(X0.train)
      ## estimated eigenscore of test data
      eigenscore0.test=as.matrix(scaleX0.test)%*%est_Eigfun0
      ## design matrix of test data
      designmatrix.test=data.frame(eigenscore0.test,id=id.test)
      colnames(designmatrix.test)=c(paste0("dim0.",1:comb[AICselectindex,1]),"id")
      # combined test data
      comdf.test = designmatrix.test %>% dplyr::left_join(clinical_info, by=c("id"="slide_id"))
      ## risk prediction
      risk.pred.fcoxph[ii]=predict(model.fcoxph,comdf.test,type="risk")
    } else{
      # eigenfunction
      est_Eigfun0=Eigvec0[,1:comb[AICselectindex,1]]
      est_Eigfun1=Eigvec1[,1:comb[AICselectindex,2]]
      
      # eigenscore
      eigenscore0=scaleX0.train%*%est_Eigfun0
      eigenscore1=scaleX1.train%*%est_Eigfun1
      
      designmatrix.train=data.frame(eigenscore0,eigenscore1,id=id.train)
      colnames(designmatrix.train)=c(paste0("dim0.",1:comb[AICselectindex,1]),paste0("dim1.",1:comb[AICselectindex,2]),"id")
      
      # combined data
      comdf.train = designmatrix.train %>% dplyr::left_join(clinical_info, by=c("id"="slide_id"))
      
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
      comdf.test = designmatrix.test %>% dplyr::left_join(clinical_info, by=c("id"="slide_id"))
      ## risk prediction
      risk.pred.fcoxph[ii]=predict(model.fcoxph,comdf.test,type="risk")
    }
    
  }
  
  # FcoxPH
  ## risk data
  risk.df.fcoxph=data.frame(id=rownames(X0),risk=risk.pred.fcoxph)
  risk.mean.fcoxph = clinical_info %>% 
    dplyr::left_join(risk.df.fcoxph, by=c("slide_id"="id")) %>% 
    dplyr::group_by(patient_id) %>% # mean predicted risk per patient
    dplyr::summarize(meanrisk=mean(risk),
                     survival_time_new=mean(survival_time_new),
                     dead=mean(dead))
  risk.mean.fcoxph$group=ifelse(risk.mean.fcoxph$meanrisk<=median(risk.mean.fcoxph$meanrisk),"low","high")
  
  if (length(unique(risk.mean.fcoxph$group))==2) {
    risk.result = survival::survdiff(Surv(survival_time_new, dead) ~ group, data = risk.mean.fcoxph)
    p.val = pchisq(risk.result$chisq, length(risk.result$n)-1, lower.tail = FALSE)
  } else p.val=0
  
  return(p.val)
}

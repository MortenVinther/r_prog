pick_species<-function(ps=c(1), inp,smsConf=0L) {
 #test  ps<-c(7L,8L); ps<-c(1L,8L); inp=inp_all; smsConf=0L
  
  ps<-sort(unique(ps))
  nps<-length(ps)
  parameters<-inp[["parameters"]]
  data<-inp[["data"]]
  
  p<-parameters
  d<-data
  
  if (smsConf==0) d$info<-data$info[ps,,drop=FALSE]
  d$info[,'s']<-1L:length(ps)
  d$nSpecies<-length(ps)
  
  d$nAges<-max(d$info[,'last-age'])-data$minAge+1 
  ages<-1:d$nAges
  d$spNames=rownames(d$info)                
  d$fbarRange<-data$fbarRange[ps,,drop=FALSE]     
  d$useRho<-data$useRho[ps]
  d$stockRecruitmentModelCode<-data$stockRecruitmentModelCode[ps]
  d$zeroCatchYearExistsSp<-data$zeroCatchYearExistsSp[ps]
  d$zeroCatchYearExists<-any(d$zeroCatchYearExistsSp==1)  
  d$zeroCatchYear<-data$zeroCatchYear[ps]
  d$seasonalCatches<-data$seasonalCatches[ps]
  cut_tab<-function(tab,reNumber=FALSE,surv=FALSE) {
   if (surv) tab<-tab[data$keySurvey.overview[,"s"] %in% ps,ages,drop=FALSE]  else tab<-tab[ps,ages,drop=FALSE] 
   if (reNumber) {
     k<-sort(unique(tab[tab>0]))
     kk<-1:length(k)
     names(kk)<-k
     tab[tab>0]<-kk[as.character(tab[tab>0])]
   }
   tab
  }
  
  cutFromTo<-function(tab){
    tab<-tab[ps,,drop=FALSE]
    n<-tab[,2]-tab[,1]+1L
    tab[,2]<-cumsum(n)
    tab[,1]<-tab[,2]-n+1L
    tab
  }
  
  d$keyLogFsta<-cut_tab(data$keyLogFsta,reNumber=FALSE)   
  d$seasonalCatches<-data$seasonalCatches[ps]
  d$keyLogFstaSd<-cut_tab(data$keyLogFstaSd,reNumber=TRUE)

  d$nlogFfromTo<-cutFromTo(data$nlogFfromTo)
  d$nlogNfromTo<-cutFromTo(data$nlogNfromTo)
  d$nlogF <- d$nlogFfromTo[,2]- d$nlogFfromTo[,1]+ 1L
  d$nlogN <- d$nlogNfromTo[,2]- d$nlogNfromTo[,1]+ 1L
  
  cut_tab_3<-function(tab,reNumber=FALSE) {
    tab<-tab[ps] 
    if (reNumber) {
      x<-  do.call(c,lapply(tab,function(x) unique(as.vector(x))))
      k<-x[x>0]  
      kk<-1:length(k)
      names(kk)<-k
      tab<-lapply(tab,function(x)  {
        dims<-dim(x); 
        x[x>0]<-kk[as.character(x[x>0])]
        matrix(x,ncol=dims[[2]])
      })
    }
    tab
  }
  d$keyLogSeparF<-cut_tab_3(tab=data$keyLogSeparF,reNumber=TRUE) 
  

                             
  d$keyVarLogN<-cut_tab(data$keyVarLogN,reNumber=TRUE)
  d$keyVarObsCatch<-cut_tab(data$keyVarObsCatch,reNumber=TRUE)
  d$keyVarObsSurvey<-cut_tab(tab=data$keyVarObsSurvey,reNumber=TRUE,surv=TRUE)
  d$keyCatchability<-cut_tab(data$keyCatchability,reNumber=TRUE,surv=TRUE)
 
  d$propMat<-data$propMat[ps]; names(d$propMat)<-1:nps
  d$stockMeanWeight<-data$stockMeanWeight[ps]; names(d$stockMeanWeight)<-1:nps
  d$catchMeanWeight<-data$catchMeanWeight[ps]; names(d$catchMeanWeight)<-1:nps
  d$catchNumber<-data$catchNumber[ps]; names(d$catchNumber)<-1:nps
  
  d$seasFprop<-data$seasFprop[ps]; names(d$seasFprop)<-1:nps
  d$natMor<-data$natMor[ps]; names(d$natMor)<-1:nps
  d$propF<-data$propF[ps]; names(d$propF)<-1:nps
  d$propM<-data$propM[ps]; names(d$propM)<-1:nps
  
  d$keyCatch<-data$keyCatch[data$keyCatch[,'s'] %in% ps,]
    k<-sort(unique(d$keyCatch[,'keyVarObsCatch']))
    kk<-1:length(k)
    names(kk)<-k
    d$keyCatch[,'keyVarObsCatch']<-kk[as.character(d$keyCatch[,'keyVarObsCatch'])]
    
    k<-ps
    kk<-1:length(ps)
    names(kk)<-k
    d$keyCatch[,'s']<-kk[as.character(d$keyCatch[,'s'])]
    
  d$logCatchObs<-data$logCatchObs[d$keyCatch[,'obs.no']]
  d$keyCatch[,'obs.no']<-1:dim(d$keyCatch)[[1]]
  
 
  # survey
  d$keySurvey.overview<-data$keySurvey.overview[data$keySurvey.overview[,'s'] %in% ps,]
  d$fleetNames<-rownames(d$keySurvey.overview)
  foundTC<-d$keySurvey.overview[,'techCreep']>0
  d$keySurvey.overview[foundTC,'techCreep']<-1:sum(foundTC)
  if (sum(data$keySurvey.overview[,'s'] %in% ps) ==1) {
    d$keySurvey.overview<-matrix(d$keySurvey.overview,nrow=1)
    colnames(d$keySurvey.overview)<- colnames(data$keySurvey.overview)
  }
  d$sampleTimeWithinSurvey<-data$sampleTimeWithinSurvey[d$keySurvey.overview[,'f']]
  
  d$keySurvey<-data$keySurvey[data$keySurvey[,'s'] %in% ps,,drop=FALSE]
  k<-sort(unique(d$keySurvey[,'keyVarObsSurvey']))
  kk<-1:length(k)
  names(kk)<-k
  d$keySurvey[,'keyVarObsSurvey']<-kk[as.character(d$keySurvey[,'keyVarObsSurvey'])]
  
  k<-sort(unique(d$keySurvey[,'keyCatchability']))
  kk<-1:length(k)
  names(kk)<-k
  d$keySurvey[,'keyCatchability']<-kk[as.character(d$keySurvey[,'keyCatchability'])]
  d$logSurveyObs<-data$logSurveyObs[d$keySurvey[,'obs.no']]             
 
  k<-d$keySurvey[,'obs.no']
  d$keySurvey[,'obs.no']<-1:length(d$logSurveyObs)
  kk<-d$keySurvey[,'obs.no']
  names(kk)<-k
  d$keySurvey.overview[,'first']<-kk[as.character(d$keySurvey.overview[,'first'])]
  d$keySurvey.overview[,'last']<-kk[as.character(d$keySurvey.overview[,'last'])]
  
   
  k<-ps
  kk<-1:length(ps)
  names(kk)<-k
  d$keySurvey[,'s']<-kk[as.character(d$keySurvey[,'s'])]
  d$keySurvey.overview[,'s']<-kk[as.character(d$keySurvey.overview[,'s'])]
  
  k<-sort(unique(d$keySurvey[,'f']))
  kk<-1:length(k)
  names(kk)<-k
  d$keySurvey[,'f']<-kk[as.character(d$keySurvey[,'f'])]
  d$keySurvey.overview[,'f']<-kk[as.character(d$keySurvey.overview[,'f'])]
  
  d$recruitYears<-data$recruitYears[ps,,drop=FALSE]
 
  if (smsConf==0) {
    d$consum<-NULL
    d$meanL<-NULL
    d$propM2<-NULL
    d$natMor1<-NULL
    d$otherN<-NULL
    d$stom<-NULL
    d$otherFood<-NULL
    d$predPreySize<-NULL
    d$vulneraIdx<-NULL
  }
  
  
  #parameters
  if (smsConf==0) {  # single species
    p$vulnera<-rep(0.0,0)
  }
    
  p$logSdLogObsCatch<-parameters$logSdLogObsCatch[1:max(d$keyVarObsCatch)]
  
  p$logCatchability<-parameters$logCatchability[1:max(d$keyCatchability)] 
  p$logSdLogObsSurvey<-parameters$logSdLogObsSurvey[1:max(d$keyVarObsSurvey)] 
  p$logSdLogFsta<-parameters$logSdLogFsta[1:max(d$keyLogFstaSd)]  
  p$logSdLogN<-parameters$logSdLogN[1:max(d$keyVarLogN)]  
  
  x<-do.call(c,lapply(d$keyLogSeparF,function(x) unique(as.vector(x)))); x<-x[x>0]
  
  if (length(x)>0) p$logSeparF<-parameters$logSeparF[1:length(x)] else  p$logSeparF=numeric(0)
  p$rho<-parameters$rho[ps] 
  p$rec_loga<-parameters$rec_loga[ps]         
  p$rec_logb<-parameters$rec_logb[ps] 
  if (sum(foundTC)>0) p$logTechCreep<-parameters$logTechCreep[1:sum(foundTC)] else p$logTechCreep<-numeric(0)
  p$Un<-parameters$Un[1:sum(d$nlogN),,drop=FALSE]   
  p$Uf<-parameters$Uf[1:sum(d$nlogF),,drop=FALSE]  
  
  list(data=d,parameters=p)
}

pick_species<-function(ps=c(1), inp) {
 #test  ps<-c(2L,9L); 
  
  ps<-sort(unique(ps))
  nps<-length(ps)
  parameters<-inp[["parameters"]]
  data<-inp[["data"]]
  
  p<-parameters
  d<-data
  
  d$info<-data$info[ps,,drop=FALSE]
  d$nSpecies<-length(ps)
  d$doSpecies <-ps
  
  d$nAges<-max(d$info[,'last-age'])-data$minAge+1 
  ages<-1:d$nAges
  d$spNames=rownames(d$info)                
  d$fbarRange<-data$fbarRange[ps,,drop=FALSE]     
  
  d$stockRecruitmentModelCode<-data$stockRecruitmentModelCode[ps]
  d$zeroCatchYearExistsSp<-data$zeroCatchYearExistsSp[ps]   
  d$zeroCatchYearExists<-any(d$zeroCatchYearExistsSp==1)  
  d$zeroCatchYear<-data$zeroCatchYear[ps]
  
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
  d$keyLogFstaSd<-cut_tab(data$keyLogFstaSd,reNumber=TRUE)
  d$nlogF <-max(d$keyLogFstaSd)
  d$nlogFfromTo<-cutFromTo(data$nlogFfromTo)
  d$nlogNfromTo<-cutFromTo(data$nlogNfromTo)
  d$nlogN<-max(d$nlogNfromTo) 
  d$keyVarLogN<-cut_tab(data$keyVarLogN,reNumber=TRUE)
  d$keyVarObsCatch<-cut_tab(data$keyVarObsCatch,reNumber=TRUE)
  d$keyVarObsSurvey<-cut_tab(data$keyVarObsSurvey,reNumber=TRUE,surv=TRUE)
  d$keyCatchability<-cut_tab(data$keyCatchability,reNumber=TRUE,surv=TRUE)
 
  d$propMat<-data$propMat[ps]
  d$stockMeanWeight<-data$stockMeanWeight[ps]
  d$catchMeanWeight<-data$catchMeanWeight[ps]
  d$seasFprop<-data$seasFprop[ps]
  d$natMor<-data$natMor[ps]
  d$propF<-data$propF[ps]
  d$propM<-data$propM[ps]
  
  d$keyCatch<-data$keyCatch[data$keyCatch[,'s'] %in% ps,]
    k<-sort(unique(d$keyCatch[,'keyVarObsCatch']))
    kk<-1:length(k)
    names(kk)<-k
    d$keyCatch[,'keyVarObsCatch']<-kk[as.character(d$keyCatch[,'keyVarObsCatch'])]
    
    k<-ps
    kk<-1:length(ps)
    names(kk)<-k
    d$keyCatch[,'s']<-kk[as.character(d$keyCatch[,'s'])]
    
#  summary(d$keyCatch)
  d$logCatchObs<-data$logCatchObs[d$keyCatch[,'obs.no']]
  d$keyCatch[,'obs.no']<-1:dim(d$keyCatch)[[1]]
  
  d$keySurvey<-data$keySurvey[data$keySurvey[,'s'] %in% ps,]
  k<-sort(unique(d$keySurvey[,'keyVarObsSurvey']))
  kk<-1:length(k)
  names(kk)<-k
  d$keySurvey[,'keyVarObsSurvey']<-kk[as.character(d$keySurvey[,'keyVarObsSurvey'])]
  
  k<-sort(unique(d$keySurvey[,'keyCatchability']))
  kk<-1:length(k)
  names(kk)<-k
  d$keySurvey[,'keyCatchability']<-kk[as.character(d$keySurvey[,'keyCatchability'])]
  d$logSurveyObs<-data$logSurveyObs[d$keySurvey[,'obs.no']]             
  d$keySurvey[,'obs.no']<-1:length(d$logSurveyObs)
  k<-ps
  kk<-1:length(ps)
  names(kk)<-k
  d$keySurvey[,'s']<-kk[as.character(d$keySurvey[,'s'])]

  k<-sort(unique(d$keySurvey[,'f']))
  kk<-1:length(k)
  names(kk)<-k
  d$keySurvey[,'f']<-kk[as.character(d$keySurvey[,'f'])]
  
  d$keySurvey.overview<-data$keySurvey.overview[data$keySurvey.overview[,'s'] %in% ps,]  

  p$logSdLogObsCatch<-parameters$logSdLogObsCatch[1:max(d$keyVarObsCatch)]
  
  p$logCatchability<-parameters$logCatchability[1:max(d$keyCatchability)] 
  p$logSdLogObsSurvey<-parameters$logSdLogObsSurvey[1:max(d$keyVarObsSurvey)] 
  p$logSdLogFsta<-parameters$logSdLogFsta[1:max(d$keyLogFstaSd)]  
  p$logSdLogN<-parameters$logSdLogN[1:max(d$keyVarLogN)]      
  p$rho<-parameters$rho[ps]
  p$rec_loga<-parameters$rec_loga[ps]         
  p$rec_logb<-parameters$rec_logb[ps]         
  p$Un<-parameters$Un[1:d$nlogN,,drop=FALSE]   
  p$Uf<-parameters$Uf[1:d$nlogF,,drop=FALSE]  
  
  list(data=d,parameters=p)
}

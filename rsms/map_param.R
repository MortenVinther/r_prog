# adjust if the are species/year combination with zero catches (or assumed F is very low and highly uncertain)
if (data$zeroCatchYearExists==1 |length(parameters$logSeparF)>0) {
  UfMap<-matrix(1L:(dim(parameters$Uf)[[1]]*dim(parameters$Uf)[[2]]),nrow=sum(data$nlogF),byrow=TRUE)
  for (s in 1:data$nSpecies) if (length(data$zeroCatchYear[[s]]) >0 ) {
    zy<-data$zeroCatchYear[[s]]
    fromTo<-data$nlogFfromTo[s,]
    UfMap[fromTo[1]:fromTo[2],zy]<-NA
    parameters$Uf[fromTo[1]:fromTo[2],zy]<-log(0.001)
  }
  if (length(parameters$logSeparF)>0) {
    ii<-data$info[data$info[,'fSepar']<99,]
    fromTo<-data$nlogFfromTo[ii[,'s'],][,1]
    
    parameters$Uf[fromTo,1]<-0.0  
    UfMap[fromTo,1]<-NA
  }
  
  UfMap<-factor(UfMap)
}

if (data$zeroCatchYearExists==1 | length(parameters$logSeparF)>0) my.map<-list(Uf=UfMap) else my.map=list()

if (any(data$stockRecruitmentModelCode %in% c(0,3))) { #random walk recruitment, no need for recruitment parameters
  aMap<-1L:length(parameters$rec_loga)
  bMap<-1L:length(parameters$rec_logb)
  aMap[data$stockRecruitmentModelCode==0]<-NA
  bMap[data$stockRecruitmentModelCode==0]<-NA
  bMap[data$stockRecruitmentModelCode==3]<-NA
  aMap<-factor(aMap)
  bMap<-factor(bMap)
  my.map<-c(my.map,list(rec_loga=aMap),list(rec_logb=bMap))
}

if (any(!data$useRho)) {
  rhoMap<-1L:length(data$useRho)
  rhoMap[data$useRho==0]<-NA
  rhoMap<-factor(rhoMap)
  my.map<-c(my.map,list(rho=rhoMap))
}

if (FALSE & length(parameters$logSeparF)>0) {
  sepMap<-1L:length(parameters$logSeparF)
  for (s in 1:data$nSpecies) {
    if (data$info[s,'fModel']==2) {
      fix<-max(unique(as.vector(data$keyLogSeparF[[s]])))
      sepMap[fix]<-NA
    }
  } 
  sepMap<-factor(sepMap)
  my.map<-c(my.map,list(logSeparF=sepMap))
}

#grep('logSeparF',nl)

# #sandeel test
# i<-grep('logSeparF',nl)[2]
# lower[i]<-0
# upper[i]<-0

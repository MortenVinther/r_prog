lock_param<-function(data,parameters,my.map,lockP){
  for (p in lockP) {
    if (p %in% names(my.map)) {
      my.map[[p]][]<-NA
      my.map[[p]]<- factor(my.map[[p]])
    } else {
      m<-1L:length(parameters[[p]])
      m[]<-NA
      my.map[[p]]<-factor(m)
    }
  }
  return(my.map)
}

map_param<-function(data,parameters) {
  
  my.map<-lapply(parameters,function(x) {if (length(x)>0) a<-factor(1L:length(x)) else a<-factor(numeric(0)); a})
  
  # adjust if the are species/year combination with zero catches (or assumed F is very low and highly uncertain)
  if (data$zeroCatchYearExists==1 |length(parameters$logSeparF)>0) {
    UfMap<-matrix(1L:(dim(parameters$Uf)[[1]]*dim(parameters$Uf)[[2]]),nrow=sum(data$nlogF),byrow=TRUE)
    for (s in 1:data$nSpecies) if (length(data$zeroCatchYear[[s]]) >0 ) {
      zy<-data$zeroCatchYear[[s]]
      fromTo<-data$nlogFfromTo[s,]
      UfMap[fromTo[1]:fromTo[2],zy]<-NA
      parameters$Uf[fromTo[1]:fromTo[2],zy]<-log(0.001)
    }
    if (length(parameters$logSeparF)>0) for (s in 1:data$nSpecies){
      fSepar<-data$info[s,'fSepar']
      if (fSepar<99) {
        ii<- which.max(data$nlogFfromTo==fSepar)
        #cat('ii:',ii,'\n')
        parameters$Uf[ii,1]<-0.0  
        UfMap[ii,1]<-NA
      }
    }
    UfMap<-factor(UfMap)
  }
 
  doNarrowProcessN<-FALSE  #VIRKER IKKE !
  if (doNarrowProcessN) {
    #set process noise (except for recruits , to a small number)
    UnMap<-matrix(1L:(dim(parameters$Un)[[1]]*dim(parameters$Un)[[2]]),nrow=sum(data$nlogN),byrow=TRUE)
    row<-0
    for (s in 1:data$nSpecies) {
      for (a in 1:data$nlogN[s] ) {
        row<-row+1
        if (a>1) {
         # cat('s:',s,' a:',a ,'row:',row,'\n')
          parameters$Un[row,4:49]<-log(0.1)
          UnMap[row,]<-NA
        }
      }
    }
    UnMap<-factor(UnMap)
  }
  
  if (data$zeroCatchYearExists==1 | length(parameters$logSeparF)>0) my.map<-purrr::assign_in(my.map,'Uf',UfMap)
  if (doNarrowProcessN) my.map<-purrr::assign_in(my.map,'Un',UnMap)
  
  # my.map<-list(Uf=UfMap) else my.map=list()
  
  if (any(data$stockRecruitmentModelCode %in% c(0,3))) { #random walk recruitment, no need for recruitment parameters
    aMap<-1L:length(parameters$rec_loga)
    bMap<-1L:length(parameters$rec_logb)
    aMap[data$stockRecruitmentModelCode==0]<-NA
    bMap[data$stockRecruitmentModelCode==0]<-NA
    bMap[data$stockRecruitmentModelCode==3]<-NA
    aMap<-factor(aMap)
    bMap<-factor(bMap)
    #my.map<-c(my.map,list(rec_loga=aMap),list(rec_logb=bMap))
    my.map<-purrr::assign_in(my.map,'rec_loga',aMap)
    my.map<-purrr::assign_in(my.map,'rec_logb',bMap)
  }
  
  if (any(data$useRho<=0)) {
    rhoMap<-1L:length(data$useRho)
    rhoMap[data$useRho==0]<-NA
    rhoMap<-factor(rhoMap)
    my.map<-purrr::assign_in(my.map,'rho',rhoMap)
  }

  if (data$sms.mode>0) {
    if( any(!data$usestomObsVar) ) {
      sOVMap<-1L:length(data$usestomObsVar)
      sOVMap[!data$usestomObsVar]<-NA
      sOVMap<-factor(sOVMap)
      my.map<-purrr::assign_in(my.map,'logStomObsVar',sOVMap)
    }
  
    if (length(my.map$vulnera) >0) {
      vulMap<-1:length(parameters$vulnera)
      vulNo<-data$vulneraIdx[data$otherFoodName,]
      vulNo<-vulNo[vulNo>0]
      vulMap[vulNo]<-NA
      vulMap<-factor(vulMap)
      my.map<-purrr::assign_in(my.map,'vulnera',vulMap)
    }
  }  
  # 
  if (length(parameters$logSeparF)>0) {
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
  return(my.map)
}

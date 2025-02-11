
func <- function(parameters) {
  square <- function(x){x*x}
  ### Parameter transform
  # f <- function(x) {2/(1 + exp(-2 * x)) - 1}
 
  getAll(data, parameters,warn=TRUE)

  
  ## Optional (enables extra RTMB features) 
  #obs<-OBS(obs) #statement tells RTMB that obs is the response. This is needed to enable automatic simulation and residual calculations from the model object.
  #logCatchObs<-OBS(logCatchObs)
  #logSurveyObs<-OBS(logSurveyObs)
 
  
  # first we have all the N states
  logN<-list()
  for (s in seq_len(nSpecies)) logN<-c(logN, list(Un[(nlogNfromTo[s,1]:nlogNfromTo[s,2]),,drop=FALSE]))

  # and then the F states
  # with seasonal data, logF is redefined as the sum of seasonal F values (which is not the same as the annual F)
  
  logF<-list()
  for (s in seq_len(nSpecies)) logF<-c(logF, list(Uf[(nlogFfromTo[s,1]:nlogFfromTo[s,2]),,drop=FALSE]))
  
  timeSteps <- nYears
  stateDimF <- nlogF
  stateDimN <- nlogN
  sdLogFsta <- exp(logSdLogFsta)
  varLogN <- exp(logSdLogN*2)
  varLogObsCatch <- exp(logSdLogObsCatch*2)
  maxAgePlusGroup <-info[,'+group']
  varLogObsSurvey = exp(logSdLogObsSurvey*2)
  
  
  makeVar2<-function() {
    d<-lapply(1:nSpecies,function(x) {
      matrix(0,nrow=info[x,'la'],ncol=timeSteps,dimnames=list(a=1:info[x,'la'],y=years))
    })
    names(d)<-spNames
    d
  }
  makeVar2sp<-function() {
    d<-lapply(1:max(yqIdx),function(x) {
      matrix(0,nrow=nSpecies,ncol=nAges,dimnames=list(sp=spNames,a=minAge:(nAges-1)))
    })
    names(d)<- paste(rep(1:nYears,each=nSeasons),rep(1:nSeasons,times=nYears),sep='_')
    d
  }
  makeVar2spAll<-function() {
    d<-lapply(1:max(yqIdx),function(x) {
      matrix(0,nrow=nSpecies+nOthSpecies,ncol=nAges,dimnames=list(sp=allSpNames,a=minAge:(nAges-1)))
    })
    names(d)<- paste(rep(1:nYears,each=nSeasons),rep(1:nSeasons,times=nYears),sep='_')
    d
  }
  
  makeVar3<-function() {
    d<-lapply(1:nSpecies,function(x) {
      array(0,dim=c(timeSteps,nSeasons,info[x,'la']),dimnames=list(y=years,q=1:nSeasons,a=1:info[x,'la']))
    })
    names(d)<-spNames
    d
  }
  
  makeVar3all<-function() {
    d<-lapply(1:(nSpecies+nOthSpecies),function(x) {
      array(0,dim=c(timeSteps,nSeasons,info[x,'la']),dimnames=list(y=years,q=1:nSeasons,a=1:info[x,'la']))
    })
    if (nOthSpecies>0) names(d)<-c(spNames,othspNames) else names(d)<-spNames
    d
  }
  
  
  makeVar3allLogical<-function() {
    d<-lapply(1:(nSpecies+nOthSpecies),function(x) {
      array(FALSE,dim=c(timeSteps,nSeasons,info[x,'la']),dimnames=list(y=years,q=1:nSeasons,a=1:info[x,'la']))
    })
    if (nOthSpecies>0) names(d)<-c(spNames,othspNames) else names(d)<-spNames
    d
  }
  
  logNq<-makeVar3all()
  logNbarq<-makeVar3()
  Zq<-makeVar3()
  
  MM<-makeVar3()  # natural mortality
  for (i in seq_along(MM)) MM[[i]]<-natMor[[i]]
 
  #Cq<-makeVar3()
  Chat<-makeVar2()
  predN<-makeVar2()

  
  ssb<-matrix(0,nrow=nSpecies,ncol=nYears,dimnames=list(species=spNames,year=years))
  Fbar<-matrix(0,nrow=nSpecies,ncol=nYears,dimnames=list(species=spNames,year=years))
  recruit<-matrix(0,nrow=nSpecies,ncol=nYears,dimnames=list(species=spNames,year=years))
  
  fq<-1L;     # first season (e.g. quarter)
  lq<-nSeasons #  use 1 for annual data
  la<-info[,'la']

  nlls<-matrix(0,ncol=4,nrow=nSpecies,dimnames=list(species=spNames,nll=c("catch","F","N","survey")))

  
  # survey initialization
  predSurveyObs<-numeric(length(logSurveyObs))
  surveyType<-keySurvey.overview[,'type']
  mina<-keySurvey.overview[,'mina']
  maxa<-keySurvey.overview[,'maxa']
  techCreepIdx<-keySurvey.overview[,'techCreep']
  
  
  ## multispecies 
  if (sms.mode>0) {
    #M2<-makeVar3()  # predation mortality (M2)
    availFood<-makeVar2spAll()
    M2<-makeVar2sp()
    
    nSuit<-dim(suitIdx)[[1]]
    #suit<-rep(0.0,nSuit)
    otherNpositiv<-makeVar3allLogical()
    stomVarDist<-formatC(info[,'stomachVariance'],width=2,flag='0')
    if (nOthSpecies >0) for (ii in ((nSpecies+1):(nSpecies+nOthSpecies))) {otherNpositiv[[ii]]<-otherN[[ii-nSpecies]] >0}
     
    if (nOthSpecies >0) for (ii in ((nSpecies+1):(nSpecies+nOthSpecies))) {
     for (iy in (1:timeSteps)) {
       for (iq in (1:nSeasons)) {
         for (ia in (info[ii,'faf']:info[ii,'la'])) {
           if (nOthSpecies >0) if (otherNpositiv[[ii]][iy,iq,ia]) logNq[[ii]][iy,iq,ia]<-log(otherN[[ii-nSpecies]][iy,iq,ia])   
    }}}}
    logOtherFood<-log(otherFood+1)
    otherFoodIdx<-nSpecies+1L
    nlls<-matrix(0,ncol=5,nrow=length(allSpNames),dimnames=list(species=allSpNames,nll=c("catch","F","N","survey",'Stomach')))
    
  } else M2<-0
  
  ## end multispecies
  
  
  ## Initialize joint negative log likelihood
  ans <- 0
  
  # needed for simulation  and consistencyCheck ? does not work
  # i<-1;
  # for (s in 1:nSpecies) {faf=info[s,'faf']; for (a in nlogFfromTo[s,1]:nlogFfromTo[s,2]) ans<- ans -dnorm(logF[[s]][a,i], sdLogFsta[keyLogFstaSd[s,a+faf-1L]], log=TRUE)}
  # for (s in 1:nSpecies) {for (a in nlogNfromTo[s,1]:nlogNfromTo[s,2]) ans<- ans -dnorm(logN[[s]][a,i]-predN[[s]][a,i],  log=TRUE)}
  

  ##########################################################################################
  
  SSB_R<-function(s,y,a=1) {
    if (stockRecruitmentModelCode[s]==0 | !recruitYears[s,y]){    ## straight RW
      rec = logN[[s]][a, y-1]
    } else {
      if (stockRecruitmentModelCode[s]==1){ ## Ricker
        rec<-rec_loga[s]+log(ssb[s,y-recAge])-exp(rec_logb[s])*ssb[s,y-recAge]
      } else {
        if(stockRecruitmentModelCode[s]==2){  ## B&H
          rec<-rec_loga[s]+log(ssb[s,y-recAge])-log(1+exp(rec_logb[s])*ssb[s,y-recAge])
        } else {
          if(stockRecruitmentModelCode[s]==3){  ## GM
            rec<-rec_loga[s]
          } else {
            stop(paste0("SR model code ",stockRecruitmentModelCode[s]," not recognized"))
          }
        }
      }
    }
    return(rec)
   }

  
  #function for Fishing mortality
  fiMo<-function(s,y,a,q=0L,expIt=TRUE) {
    if (y %in% zeroCatchYear[[s]]) {
      if (expIt) return(0.0) else return(-10)
    } else {
        ff<-logF[[s]][keyLogFsta[s,a],y]  # annual F
        if (fModel==2 & a >=fSepar)  ff<-ff+logSeparF[keyLogSeparF[[s]][y,a]]
        if  (q==0L) { 
           if (expIt) return(exp(ff)) else return(ff)
        } else {
          if (seasFprop[[s]][y,q,a] > 0) ff<-ff+logSeasFprop[[s]][y,q,a] 
          if (expIt) return(exp(ff)) else return(ff)
        }
   }}

 
  ###################  now we begin 
  for (s in seq_len(nSpecies)) {
  
    cat('Species:',s,spNames[s],'\n')
    fModel <- as.integer(info[s,'fModel']) 
    fSepar <- as.integer(info[s,'fSepar'])
    
    ## First take care of F
    fsd <- sdLogFsta[keyLogFstaSd[s,keyLogFstaSd[s,]>0]]

    if (useRho[s]==0) {  # no correlation
         fvar <- diag(fsd[1:stateDimF[s]]*fsd[1:stateDimF[s]],ncol=stateDimF[s], nrow=stateDimF[s])
    } else if (useRho[s]==1) {  # compound_symmetry
         fcor <- outer(1:stateDimF[s],
                       1:stateDimF[s],
                       function(i,j)(i!=j)*rho[s] + (i==j))
         fvar <- outer(1:stateDimF[s],
                       1:stateDimF[s],
                       function(i,j)fcor[cbind(i,j)]*fsd[i]*fsd[j])  
    } else if (useRho[s]==2) {  # AR(1)
         fcor<-rho^(abs(matrix(1:stateDimF[s] - 1, nrow = stateDimF[s], ncol = stateDimF[s], byrow = TRUE) - 1:stateDimF[s] - 1))
         fvar <- outer(1:stateDimF[s],
                       1:stateDimF[s],
                       function(i,j)fcor[cbind(i,j)]*fsd[i]*fsd[j])  
     } else  stop(paste0("No handler for useRho", useRho[s],'\n'))
    
    if (zeroCatchYearExistsSp[s]==1)  llh<- sum(dmvnorm( diff( t(logF[[s]][,-zeroCatchYear[[s]]])),mu=0, Sigma=fvar, log=TRUE)) else {
                                      llh<- sum(dmvnorm( diff( t(logF[[s]])),mu=0, Sigma=fvar, log=TRUE))
    }
    ans<-ans-llh
    nlls[s,"F"]<- nlls[s,"F"] - llh
    
   # Spawning Stock Biomass

  # # SSB 1. January !!   
  # for (s in 1:nSpecies) { 
  #   for (y in (1:nYears)) {
  #     ssb[s,y]<-sum(exp(logN[[s]][,y]) *propMat[[s]][y,q,]*stockMeanWeight[[s]][y,q,])  #HUSK AT SÆTTE PropMat for ikke rekruterede aldre til 0
  #   }
  # }  
  
  # SSB 1. January !!
  #  Class 'simref' does not allow operation 'sum' (argument must be fully simulated)
    for (y in seq_len(nYears)) {
      ssb[s,y]<- sum(exp(logN[[s]][,y])*propMat[[s]][y,spawnSeason,]*stockMeanWeight[[s]][y,spawnSeason,])
    
    #  for (a in seq_len(stateDimN[s])) {
    #    ssb[s,y]<- ssb[s,y] + exp(logN[[s]][a,y])*propMat[[s]][y,spawnSeason,a]*stockMeanWeight[[s]][y,spawnSeason,a]
    #  }
    }

  # Natural mortalities
    if (sms.mode>0) {
      yq<-0L
      for (y in seq_len(nYears)) {
        for (q in seq_len(nSeasons)) {
          yq<-yq+1L
          for (a in seq_len(stateDimN[s])) { 
            MM[[s]][y,q,q]<-natMor1[[s]][y,q,q]+M2[[yq]][s,a]
    }}}}
   
    # ##fishing mortality
    # for (y in seq_len(timeSteps)) {
    #   for (a in (info[s,'faf']:info[s,'la'])) {
    #     for (q in 1:lq) if (a >1 | q>=recSeason ) {
    #         if (y %in% zeroCatchYear[[s]]) {
    #           logFiq[[s]][y,q,a]<- -10  # lidt snyd, skulle være log(0)
    #         } else {
    #            logFiq[[s]][y,q,a]<-logF[[s]][keyLogFsta[s,a],y]+logSeasFprop[[s]][y,q,a] 
    #            if (fModel==2) if (a >=fSepar)  logFiq[[s]][y,q,a]<-logFiq[[s]][y,q,a]+logSeparF[keyLogSeparF[[s]][y,a]]
    #          }
    # }}}     
    
     
    
    ## Now take care of N
    nvar <- outer(1:stateDimN[s], 1:stateDimN[s],
                  function(i,j) (i==j)*varLogN[ keyVarLogN[s,i]])
   
    #recruits first year given recruitment at age 0 (later in the same year)
    if (stockRecruitmentModelCode[s] >=1 & recAge==0 & recruitYears[s,1]){
      predN[[s]][1,1]<-SSB_R(s,y=1);
      llh <- dmvnorm(logN[[s]][1,1], predN[[s]][1,1], nvar[1,1], log=TRUE) ## N-Process likelihood
      ans <- ans -llh
      nlls[s,"N"]<-  nlls[s,"N"] -llh 
    }
    
   for(i in seq.int(2L,timeSteps)) { 
      predN[[s]][1,i]<-SSB_R(s,y=i);
  
      for(j in seq.int(2L,stateDimN[s])) {
          if (keyLogFsta[s,j-1]>0) {  # to take account for ages with no F
             predN[[s]][j,i]=logN[[s]][j-1,i-1]-fiMo(s,i-1,j-1)-sum(MM[[s]][i-1,,j-1]) 
          } else { 
            predN[[s]][j,i]=logN[[s]][j-1,i-1]-sum(MM[[s]][i-1,,j-1]) 
          }
       }
       if(maxAgePlusGroup[s]==1){
          predN[[s]][stateDimN[s],i] = log( exp(logN[[s]][stateDimN[s]-1,i-1]- fiMo(s,i-1,stateDimN[s]-1)-sum(MM[[s]][i-1,,stateDimN[s]-1])) +
                                            exp(logN[[s]][stateDimN[s],i-1]  - fiMo(s,i-1,stateDimN[s])  -sum(MM[[s]][i-1,,stateDimN[s]]))  )
       }
      llh <- dmvnorm(logN[[s]][,i], predN[[s]][,i], nvar, log=TRUE) ## N-Process likelihood
      ans <- ans -llh
      nlls[s,"N"]<-  nlls[s,"N"] -llh 
    }

  #We have now calculated predN (1 quarter)
    for(i in seq_len(timeSteps)) {
      for(j in seq_len(stateDimN[s])) {
       
        #first season (q=fq)
        if (j >1 | recSeason==1 ) {  # first age might not have enter the model yet
          logNq[[s]][i,fq,j]<-logN[[s]][j,i]
          if((keyLogFsta[s,j])>0) {
            Zq[[s]][i,fq,j] <- fiMo(s,i,j,q=fq) + MM[[s]][i,fq,j]
          } else {
            Zq[[s]][i,fq,j] <- MM[[s]][i,fq,j]
          }
          logNbarq[[s]][i,fq,j] <- logNq[[s]][i,fq,j]-log(Zq[[s]][i,fq,j]) + log(1.0 -exp(-Zq[[s]][i,fq,j]))
          if (keyLogFsta[s,j]>0) Chat[[s]][j,i] <- exp(logNbarq[[s]][i,fq,j]+fiMo(s,i,j,q=fq,expIt=FALSE))
        } else Chat[[s]][j,i] = 0.1; # SNYD for at undgå log(0), værdien bruges ikke, så OK
        
        # recruits within the year, need simpler code
        if (j ==1 & recSeason>1) {
           logNq[[s]][i,recSeason,j]<-logN[[s]][j,i] 
           if (keyLogFsta[s,j]>0) {
             Zq[[s]][i,recSeason,j] <- fiMo(s,i,j,q=recSeason) + MM[[s]][i,recSeason,j]
           } else {
             Zq[[s]][i,recSeason,j] <- MM[[s]][i,recSeason,j]
           }
           logNbarq[[s]][i,recSeason,j] <- logNq[[s]][i,recSeason,j]-log(Zq[[s]][i,recSeason,j]) +log(1.0 -exp(-Zq[[s]][i,recSeason,j]))
           if (keyLogFsta[s,j]>0)  {Chat[[s]][j,i] <-  exp(logNbarq[[s]][i,recSeason,j]+fiMo(s,i,j,q=recSeason,expIt=FALSE))
           } else { Chat[[s]][j,i]<-0.1  }
        }
        
        #remaining seasons
        if (nSeasons>1) for (q in 2:lq) if (j >1 | q>recSeason ) {
          logNq[[s]][i,q,j]<- logNq[[s]][i,q-1,j]-Zq[[s]][i,q-1,j]
          if (keyLogFsta[s,j]>0) {
            Zq[[s]][i,q,j] <- fiMo(s,i,j,q) + MM[[s]][i,q,j]
          } else {
            Zq[[s]][i,q,j] <- MM[[s]][i,q,j]
          }
          logNbarq[[s]][i,q,j] <- logNq[[s]][i,q,j]-log(Zq[[s]][i,q,j]) +log(1.0 -exp(-Zq[[s]][i,q,j]))
          if (keyLogFsta[s,j]>0)  {Chat[[s]][j,i] <- Chat[[s]][j,i]+ exp(logNbarq[[s]][i,q,j]+fiMo(s,i,j,q,expIt=FALSE))}  else { Chat[[s]][j,i]<-0.1  }
        }
      }
    }
  
  Chat[[s]]<-log(Chat[[s]])
  
  ########  match to observations
  
  # first catches
    #cat('Species: ',s,'\n')
    idx<-keyCatch[,"s"]==s
    key.v<-keyCatch[idx,"keyVarObsCatch"]
    yy<-keyCatch[idx,"y"]
    aa<-keyCatch[idx,'a']
    ay<-cbind(aa,yy)
    obs<-logCatchObs[idx]
    predObs<-Chat[[s]][ay]
    var <- varLogObsCatch[key.v]
    llh<- sum(dnorm(obs,predObs,sqrt(var),log=TRUE))
    ans <- ans - llh
    nlls[s,'catch']<- nlls[s,'catch'] - llh
      

  # and now surveys
  fleets<-keySurvey.overview[keySurvey.overview[,'s']==s,'f']
    for (fl in  fleets)  {
    #cat("survey: ",fl,'\n')
    keys<-keySurvey[keySurvey[,"f"]==fl ,]
    q<-keys[1,"q"]
    if (surveyType[fl]==1) {
      for (a in mina[fl]:maxa[fl]) {
         keysA<-keys[keys[,"a"]==a  ,]
         flYears<-keysA[,'y']
         keyPowerQ<-keysA[1,"keyPowerQ"]
         keyCatchability<-keysA[1,"keyCatchability"]
         keyVarObsSurvey<-keysA[1,"keyVarObsSurvey"]
         obs.no<-keysA[,'obs.no']
         predSurveyObs[obs.no]<- logNq[[s]][flYears,q,a] - Zq[[s]][flYears,q,a]*sampleTimeWithinSurvey[fl]
         if(keyPowerQ>0) predSurveyObs[obs.no] <- predSurveyObs[obs.no]*exp(logQpow[keyPowerQ])
         predSurveyObs[obs.no]<- predSurveyObs[obs.no]+logCatchability[keyCatchability]
         var <- varLogObsSurvey[keyVarObsSurvey]
         llh<-sum(dnorm(logSurveyObs[obs.no], predSurveyObs[obs.no],sqrt(var),log=TRUE))
         ans <- ans - llh
         nlls[s,'survey']<- nlls[s,'survey']- llh
       } 
    } else if (surveyType[fl]==2) {  # exploitable biomass (assumed all age with F>0)
      flYears<-keys[,'y']
      faf<-info[s,'faf']; laf<-info[s,'last-age']
      obs.no<-keys[,'obs.no']
      keyCatchability<-keys[1,"keyCatchability"]
      keyVarObsSurvey<-keys[1,"keyVarObsSurvey"]
      # does not work     predObs<-sapply(flYears,function(y) log(sum(exp(logNq[[s]][y,q,faf:laf] - Zq[[s]][y,q,faf:laf]*sampleTimeWithinSurvey[fl])*stockMeanWeight[[s]][y,q,faf:laf]))+ logCatchability[keyCatchability] )
      for (f in (1:length(obs.no)))  {
        y<-flYears[f]
        predSurveyObs[obs.no[f]]<- log(sum(exp(logNq[[s]][y,q,faf:laf] - Zq[[s]][y,q,faf:laf]*sampleTimeWithinSurvey[fl])*stockMeanWeight[[s]][y,q,faf:laf]))+ logCatchability[keyCatchability] 
      }
      var <- varLogObsSurvey[keyVarObsSurvey]
      llh <- sum(dnorm(logSurveyObs[obs.no],predSurveyObs[obs.no],sqrt(var),log=TRUE))
      ans <- ans - llh
      nlls[s,'survey']<- nlls[s,'survey'] - llh
    } else if (surveyType[fl]==3) {  # SSB index
      flYears<-keys[,'y']
      faf<-1L; laf<-info[s,'la']
      obs.no<-keys[,'obs.no']
      keyCatchability<-keys[1,"keyCatchability"]
      keyVarObsSurvey<-keys[1,"keyVarObsSurvey"]
      for (f in 1:length(obs.no))  {
        y<-flYears[f]
        predSurveyObs[obs.no[f]]<- log(sum(exp(logNq[[s]][y,q,faf:laf] - Zq[[s]][y,q,faf:laf]*sampleTimeWithinSurvey[fl])*stockMeanWeight[[s]][y,q,faf:laf]*propMat[[s]][y,q,faf:laf]))+ logCatchability[keyCatchability] 
      }
      var <- varLogObsSurvey[keyVarObsSurvey]
      llh<- sum(dnorm(logSurveyObs[obs.no],predSurveyObs[obs.no],sqrt(var),log=TRUE))
      ans <- ans -llh
      nlls[s,'survey']<- nlls[s,'survey']-llh
    } else if (surveyType[fl]==4) {  # effort index, one "catchability" by age group
      for (a in mina[fl]:maxa[fl]) {
        keysA<-keys[keys[,"a"]==a, ]
        flYears<-keysA[,'y']
        keyPowerQ<-keysA[1,"keyPowerQ"]
        keyCatchability<-keysA[1,"keyCatchability"]
        keyVarObsSurvey<-keysA[1,"keyVarObsSurvey"]
        techCreepNo<-techCreepIdx[fl]
        obs.no<-keysA[,'obs.no']
        predSurveyObs[obs.no]<- fiMo(s,flYearss,a,q,expIt=FALSE) + logCatchability[keyCatchability] 
        if (techCreepNo>0) for (n in seq_len((obs.no)))  {
            predSurveyObs[obs.no[n]]<- predSurveyObs[obs.no[n]]+ n*logTechCreep[techCreepNo] 
        }
        var <- varLogObsSurvey[keyVarObsSurvey]
        llh<-sum(dnorm(logSurveyObs[obs.no],predSurveyObs[obs.no],sqrt(var),log=TRUE))
        ans <- ans -llh
        nlls[s,'survey']<- nlls[s,'survey'] -llh
      }
    } else if (surveyType[fl]==5) {  # effort used as index for Fbar
      flYears<-keys[,'y']
      faf<-fbarRange[s,1]; laf<-fbarRange[s,2]; naf<-laf-faf+1
      obs.no<-keys[,'obs.no']
      keyCatchability<-keys[1,"keyCatchability"]
      keyVarObsSurvey<-keys[1,"keyVarObsSurvey"]
      techCreepNo<-techCreepIdx[fl]
      for (n in 1L:length(obs.no))  {
        y<-flYears[n]
        predSurveyObs[obs.no[n]]<-0.0;
        for(j in seq.int(faf,laf)) {
           predSurveyObs[obs.no[n]]<- predSurveyObs[obs.no[n]] + fiMo(s,y,j,q) #sum seasonal F within Fbar range
        }
        predSurveyObs[obs.no[n]]<- log(predSurveyObs[obs.no[n]]/naf) + logCatchability[keyCatchability] 
        if (techCreepNo>0) predSurveyObs[obs.no[n]]<- predSurveyObs[obs.no[n]]+ n*logTechCreep[techCreepNo] 
      }
      var <- varLogObsSurvey[keyVarObsSurvey]
      llh<-sum(dnorm(logSurveyObs[obs.no],predSurveyObs[obs.no],sqrt(var),log=TRUE))
      ans <- ans -llh
      nlls[s,'survey']<- nlls[s,'survey']  -llh
    }  else stop(paste("s:",s,"  fleet:",fl,'  fleet type:',surveyType[fl],' is not known\n'))
    
  } # end fleet loop
  } #end species loop
 
  # end single species mode
 
 
 # food suitability
 suitability <-function(q,pred,prey,predSize,preySize,ratio,vulIdx,logIt=TRUE) {
   opt<-info[pred,'sizeSlct']
   overl<-overlap[q,pred,prey]
   if (opt %in% c(0L,4L)) {         # (0=uniform) no size selection, or 4=confined uniform) no size selection, but within limit
                                    # the prey/pred size ratios have already been checked in input, so no check here
     suit<-overl+vulnera[vulIdx]
     if (logIt) return(suit) else return(exp(suit))
   }
 } #end function
 
  # else {           //  size selection
 #   if (ratio >= min_pred_prey_size_ratio(pred,prey) &&  ratio <= max_pred_prey_size_ratio(pred,prey)){    
 #     
 #     if (size_selection(pred)==1  || size_selection(pred)==2 ) { //normal distribution or asymmetric normal distribution
 #       tmp=log(ratio)-(pref_size_ratio(pred)*prey_size_adjustment(prey)+pref_size_ratio_correction(pred)*log(pred_size));     
 #       return vul*exp(-square(tmp)/(2.0* var_size_ratio(pred)));
 #     }
 #     else if (size_selection(pred)==3) {      //Gamma     1.0/(b^a*gamma(a))*x^(a-1)*exp(-x/b)
 #       // use gammln function: exp(-(log(b)*a+log(gamma(a)))+log(x)*(a-1)-x/b)
 #       tmp=log(ratio);
 #       // cout<<"gamma var:"<<var_size_ratio(pred)<<" pref: "<<pref_size_ratio(pred)<<" tmp:"<<tmp<<endl;                                                   
 #       size_sel= exp(-(log(var_size_ratio(pred))*pref_size_ratio(pred)+gammln(var_size_ratio(pred)))+log(tmp)*(pref_size_ratio(pred)-1.0)-tmp/pref_size_ratio(pred));
 #       return vul*size_sel;
 #     }
 #     else if (size_selection(pred)==5 || size_selection(pred)==6) {      //beta     gamma(a+b)/gamma(a)/gamma(b)*x^(a-1)*(1-x)^(b-1)
 #       // use gammln function: exp(lgamma(a+b)-lgamma(a)-lgamma(b)+log(x)*(a-1)+log(1-x)*(b-1))                                             
 #       //rescale ratio to [0;1]                                             
 #       ratio=(log(ratio)-all_min_pred_prey_size_ratio(pred))/all_range_pred_prey_size_ratio(pred); 
 #       size_sel=exp(gammln(pref_size_ratio(pred)+var_size_ratio(pred))-gammln(pref_size_ratio(pred))-gammln(var_size_ratio(pred))+log(ratio)*(pref_size_ratio(pred)-1.0)+log(1.0-ratio)*(var_size_ratio(pred)-1));
 #       //cout<<"pred:"<<pred<<" prey:"<<prey<<setprecision(3)<<" ratio:"<<ratio<<" size_sel:"<<size_sel<<" suit:"<<vul*size_sel<<endl;
 #       return vul*size_sel;
 #     }
 #     
 #     else return -1000.0;  // error
 #   }
 #   else return 0.0;
 # }  
 # 
 
   ### Multispecies mode 
   if (sms.mode>0) {
     #update predator prey overlap when estimated
     if (dim(overlapIdx)[[1]]>0)  for (ii in seq_len(dim(overlapIdx)[[1]])) {overlap[overlapIdx[ii,"q"],overlapIdx[ii,"predNo"],overlapIdx[ii,"preyNo"]] <- overlapP[overlapIdx[ii,'overlap']]}
     
     #calculate M2
     cat('M2 calculation\n')
     for (yqa in (1:dim(suitIdx)[[1]])) {
       x1<- suitIdx[yqa,]
       y<-x1$y
       q<-x1$q
       yq<-yqIdx[y,q]
       #area<-x1$area
  #     localNq<-rbind(suppressWarnings(do.call(rbind,lapply(1:nSpecies,function(x)logNq[[x]][y,q,]))),rep(0,nAges)) # reformat prey abundance  to allow vectorization
       localNq<-rbind(do.call(rbind,lapply(1:nSpecies,function(x) c(logNq[[x]][y,q,],rep(-9,nAges-la[x])))),rep(0,nAges)) # reformat prey abundance  to allow vectorization
       #if (q==1) cat('M2 y:',y,'\n')
       M2[[yq]][,]<-0
       
       for (yqapa in seq_len(dim(x1$data[[1]])[[1]])) {
         x2<-x1$data[[1]][yqapa,]
         predNo<-x2$predNo
         predAge<-x2$predAge
         predW<-x2$predW
         localNq[otherFoodIdx,]<-logOtherFood[predNo] #predator dependent other food
         predAbun<-logNq[[predNo]][y,q,predAge]
         predCons<-consum[[predNo]][y,q,predAge]

         xx<-x2$data[[1]] %>%
           transmute(preyNo,preyAge,suit=suitability(q,pred=predNo,prey=preyNo, predSize=predW, preySize=preyW, ratio=logRatio, vulIdx=vulneraIdx),
                     availFood=exp(localNq[cbind(preyNo,preyAge)]+preyW+suit),
                     M2=exp(predAbun+suit)*predCons)
         availFood[[yq]][predNo,predAge]<-sum(xx$availFood)
         xx<-subset(xx,preyNo<=nSpecies) # exclude other food from M2 calc
         
         xx$M2<- xx$M2 / availFood[[yq]][predNo,predAge]
         partM2Idx[yqa,]$data[[1]][yqapa,]$data[[1]]<-xx
         
         pa<-cbind(xx$preyNo,xx$preyAge)
         M2[[yq]][pa]<- M2[[yq]][pa]+xx$M2
        }
     }

     DirichletLogLikelihood<-function(alfa0,EstomTot){
       alfa<-EstomTot*alfa0
       lgamma(alfa0)-sum(lgamma(alfa)) + sum((alfa-1)*log(stomconTot))
     }
     
  # Stomach contents observations
     oldY<-0
     cat('Stomach observations\n')
     for (yqa in seq_len(dim(stom)[[1]])) {
       x1<- stom[yqa,]
       y<-x1$y
       q<-x1$q
       yq<-x1$yqALK
       #area<-x1$area
       # N at size class
       nAtL<-matrix(0,nrow=nSpecies+1,ncol=maxSizeCl)
       alk1<-alk[yq,'data'][[1]][[1]]
       for (sa in  seq_len(dim(alk1)[[1]])) {
         alk2<-alk1[sa,]
         s<- alk2$s
         a<- alk2$a
         nAtL[s,alk2$minSize:alk2$maxSize]<- nAtL[s,alk2$minSize:alk2$maxSize]+exp(logNq[[s]][y,q,a])*alk2$data[[1]]$alk
       }
       for (yqapa in seq_len(dim(x1$data[[1]])[[1]])) {
         x2<-x1$data[[1]][yqapa,]
         #stom[yqa,]$data[[1]][yqapa,]
         pred<-x2$pred
         predSize<-x2$predSizeClass
         predW<-x2$predSizeW
         wDiri<-x2$noSampl
         nPreyGroups<-x2$nPreyGroups
         stomVar <- exp(logStomObsVar[x2$stomObsVarIdx])
         nAtL[otherFoodIdx,1]<-otherFood[pred] # predator dependent other food
         suit <-suitability(q,pred,prey=x2$data[[1]]$prey, predSize=predW, preySize=x2$data[[1]]$preySizeW, ratio=x2$data[[1]]$logPPsize, vulIdx=x2$data[[1]]$vulneraIdx,logIt=FALSE)
         availFoodTmp<-nAtL[cbind(x2$data[[1]]$prey,x2$data[[1]]$preySizeClass)]*x2$data[[1]]$preySizeW*suit
         #if (is.na(sum(availFoodTmp))) cat('Some ting is wrong with AvailFood', y,q,' pred',pred,'\n')
         Estom<-availFoodTmp/sum(availFoodTmp)
       
        # data.frame(Estom,Suit=suit,stocon=x2$data[[1]]$stomcon)  
         stom[yqa,]$data[[1]][yqapa,]$data[[1]]$Estom<-Estom
         stom[yqa,]$data[[1]][yqapa,]$data[[1]]$Suit<-suit
         ##stom[yqa,]$data[[1]][yqapa,]$data[[1]]$vulnera<-vulnera[x2$data[[1]]$vulneraIdx]
         stom[yqa,]$data[[1]][yqapa,]$data[[1]]$availFood<-availFoodTmp
       
         
         # stomach contents by prey (summing over prey sizes)
         if (info[pred,'sumStomLike']==1) {  # sum over size
            stomconTot<-unlist(x2$data[[1]][x2$data[[1]]$firstL,'stomconTot'],use.names=FALSE)
            EstomTot<-rep(0,nPreyGroups)
            for (i in seq_len(dim(x2$data[[1]])[[1]])) {
              ii<-unlist(x2$data[[1]][i,'preyIdx'],use.names = FALSE)
             EstomTot[ii]<-EstomTot[ii]+ Estom[i]
            }
         } else { # use observations as they are
           EstomTot<-Estom
           stomconTot<-x2$data[[1]]$stomcon
         }
 
         switch(stomVarDist[pred],
                "01" = {  # log normal distribution
                  llh<-sum(dnorm(log(stomconTot),log(EstomTot),sd=sqrt(stomVar),log=TRUE))
                },
                "02" = {  # normal distribution
                
                },
                "03" = {   # Dirichlet distribution (with estimation of 'concentration' parameter)
                  alfa0<-wDiri*stomVar - 1.0
                  llh<-DirichletLogLikelihood(alfa0,EstomTot)
                  #alfa<-EstomTot*alfa0
                  #llh<-lgamma(alfa0)-sum(lgamma(alfa)) + sum((alfa-1)*log(stomconTot))
                },
                "04" = {   # Dirichlet distribution (with input of known 'concentration' parameter)
                  alfa0<-x2$phi
                  #alfa<-EstomTot*alfa0
                  llh<-DirichletLogLikelihood(alfa0,EstomTot)
                  #llh<-lgamma(alfa0)-sum(lgamma(alfa)) + sum((alfa-1)*log(stomconTot))
                },
                "05" =  {  # Dirichlet distribution (with input of maximum 'concentration' parameter)
                  alfa0<-x2$phi*stomVar
                  #alfa<-EstomTot*alfa0
                  #llh<-lgamma(alfa0)-sum(lgamma(alfa)) + sum((alfa-1)*log(stomconTot))
                  llh<-DirichletLogLikelihood(alfa0,EstomTot)
                },
                
                stop(paste0("No handler for stomach variance", info[pred,'stomachVariance']))
         ) #end switch
         ans <- ans - llh
         nlls[pred,'Stomach']<- nlls[pred,'Stomach']  - llh
        }
     }
   }
 
 ### output  ####
 
 # make results as data.frame
 toDF1<-function(x,val='N',expIt=FALSE){
   a<-do.call(rbind,lapply(names(x),function(xx) {
     b<-array2DF(x[[xx]])
     if (expIt) b<- b %>% mutate(Value=exp(Value))
     colnames(b)<-c('year','quarter','age',val)
     b %>% mutate(age=as.integer(age)-off.age) %>%   mutate_if(is.character,as.integer) %>% mutate(species=xx)
   }))
 }

 dfNq<-toDF1(x=logNq,val='N',expIt=TRUE) 
 dfZq<-toDF1(x=Zq,val='Z',expIt=FALSE) 
 res<-left_join(dfNq,dfZq,by = join_by(year, quarter, age, species)) %>% as_tibble()

 dfFq<-left_join(dfZq %>% mutate(y=year+off.year, a=age+off.age), 
                 data.frame(species=spNames,s=1L:length(spNames)),by = join_by(species)) %>% 
                 mutate(Z=NULL)
 

 dfFq<-do.call(rbind,lapply( 1L:nSpecies,function(sp) {
    fModel <<- as.integer(info[sp,'fModel']) 
    filter(dfFq,s==sp & a>=info[sp,'faf']) %>% filter( !(y %in% zeroCatchYear[[sp]])) %>% rowwise() %>% 
    mutate(FF=fiMo(s=unlist(s),y=unlist(y),a=unlist(a),q=unlist(quarter),expIt=TRUE),a=NULL,y=NULL,s=NULL)
 }))
 
 res<-full_join(res,dfFq,by = join_by(year, quarter, age, species)) %>% as_tibble()

 
 toDF2<-function(x,val='M2',expIt=FALSE){
   a<-do.call(rbind,lapply(names(x),function(xx) {
     b<-array2DF(x[[xx]])
     if (expIt) b<- b %>% mutate(Value=exp(Value))
     colnames(b)<-c('species','age',val)
     b<- b %>% mutate(age=as.integer(age)-off.age) %>%  
       mutate(year= as.integer(strsplit(xx,'_')[[1]][1]) -off.year,quarter= as.integer(strsplit(xx,'_')[[1]][2]))

   }))
 }
 if (sms.mode>0) {
    dfM2q<-toDF2(x=M2,val='M2',expIt=FALSE) 
    res<-left_join(res,dfM2q,by = join_by(year, quarter, age, species))
    res<-left_join(res,toDF1(x=natMor1,val='M1',expIt=FALSE) %>% mutate(year=year-off.year),by = join_by(year, quarter, age, species))
 } else res$M2<- -9;

 res<-left_join(res, data.frame(species=allSpNames,s=1L:length(allSpNames)),by = join_by(species)) %>%
         filter(!(age==recAge & quarter<recSeason)) 
 
 
 toDF3<-function(x,val='predN',expIt=FALSE){
   a<-do.call(rbind,lapply(names(x),function(xx) {
     b<-array2DF(x[[xx]])
     if (expIt) b<- b %>% mutate(Value=exp(Value))
     colnames(b)<-c('age','year',val)
     b %>% mutate(age=as.integer(age)-off.age) %>%   mutate_if(is.character,as.integer) %>% mutate(species=xx)
   }))
 }
 
 toDF4<-function(x,val='SSB',expIt=FALSE){
     b<-array2DF(x)
     if (expIt) b<- b %>% mutate(Value=exp(Value))
     colnames(b)<-c('species','year',val)
     mutate(b,year=as.integer(year)) 
 }
 
 dfpredN<-toDF3(x=predN,val='predN',expIt=TRUE) 
 dfChat<-toDF3(x=Chat,val='Chat',expIt=TRUE) 
 
 resAnnual<-left_join(dfpredN,dfChat,by = join_by(year, age, species)) %>% as_tibble()
 yieldQ<-left_join(toDF1(x=catchNumber,val='Catch',expIt=FALSE) ,
           toDF1(x=catchMeanWeight,val='wCatch',expIt=FALSE),
           by = join_by(year, quarter, age, species) )  %>% mutate(year=year-off.year)
 yield<-yieldQ %>% group_by(species,year,age) %>%summarize(yield=sum(wCatch*Catch),.groups = "drop") 
 resAnnual<- left_join(resAnnual,yield,by = join_by(age, year, species))
 
 res<-left_join(res, yieldQ,by = join_by(year, quarter, age, species)) 
 

 Fbar[,]<-0
 for (s in seq_len(nSpecies)) {
   fModel <- as.integer(info[s,'fModel']) 
   recruit[s,]<-logNq[[s]][,recSeason,1]
   nFa<-fbarRange[s,2]-fbarRange[s,1]+1
   for (y in seq_len(nYears)) {
      for (a in fbarRange[s,1]:fbarRange[s,2]) {
         Fbar[s,y]<- Fbar[s,y] + fiMo(s,y,a,expIt=TRUE) 
      } 
     Fbar[s,y]<- Fbar[s,y] /nFa
   }
 }
 FBAR<-log(Fbar)
 SSB<-log(ssb)
 
 resSummary<-left_join(toDF4(x=ssb,val='SSB',expIt=FALSE),toDF4(x=recruit,val='recruit',expIt=TRUE),by = join_by(species, year)) 
 resSummary<-left_join(resSummary,toDF4(x=Fbar,val='Fbar',expIt=FALSE),by = join_by(species, year))  %>% mutate(s=match(resSummary$species,spNames))
 resSummary<-left_join(resSummary,yield%>% group_by(species,year) %>% 
              summarize(yield=sum(yield),.groups = "drop"),by = join_by(species, year))
 
 REPORT(res) 
 REPORT(resAnnual) 
 REPORT(resSummary) 
 
  ADREPORT(SSB)
  ADREPORT(recruit)
  ADREPORT(FBAR)
  
  #REPORT(predN)
  #REPORT(Zq)
  #REPORT(Chat)
  REPORT(predSurveyObs)
  if (sms.mode>0) {
    REPORT(availFood)
    REPORT(stom)
    REPORT(partM2Idx)
   # REPORT(M2)
  }
  REPORT(nlls)
  
  ans
}


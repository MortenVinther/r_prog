
func <- function(parameters) {
  square <- function(x){x*x}
  ### Parameter transform
  # f <- function(x) {2/(1 + exp(-2 * x)) - 1}
 
  getAll(data, parameters,warn=FALSE)
  
  ## Optional (enables extra RTMB features) 
  #obs<-OBS(obs) #statement tells RTMB that obs is the response. This is needed to enable automatic simulation and residual calculations from the model object.
  logCatchObs<-OBS(logCatchObs)
  logSurveyObs<-OBS(logSurveyObs)
 
  # first we have all the N states
  logN<-list()
  for (s in 1:nSpecies) logN<-c(logN, list(Un[(nlogNfromTo[s,1]:nlogNfromTo[s,2]),,drop=FALSE]))

  # and then the F states
  # with seasonal data, logF is redefined as the sum of seasonal F values (which is not the same as the annual F)
  
  logF<-list()
  for (s in 1:nSpecies) logF<-c(logF, list(Uf[(nlogFfromTo[s,1]:nlogFfromTo[s,2]),,drop=FALSE]))
  
  timeSteps <- nYears
  stateDimF <- nlogF
  stateDimN <- nlogN
  sdLogFsta = exp(logSdLogFsta)
  varLogN = exp(logSdLogN*2)
  varLogObsCatch = exp(logSdLogObsCatch*2)
  maxAgePlusGroup <-info[,'+group']
  
  varLogObsSurvey = exp(logSdLogObsSurvey*2)
  
  makeVar2<-function() {
    d<-lapply(1:nSpecies,function(x) {
      matrix(0,nrow=info[x,'la'],ncol=timeSteps,dimnames=list(a=paste0("a",1:info[x,'la']),y=years))
    })
    names(d)<-spNames
    d
  }
  
  makeVar3<-function() {
    d<-lapply(1:nSpecies,function(x) {
      array(0,dim=c(timeSteps,nSeasons,info[x,'la']),dimnames=list(y=years,q=1:nSeasons,a=paste("a",1:info[x,'la'])))
    })
    names(d)<-spNames
    d
  }
  
  logNq<-makeVar3()
  logNbarq<-makeVar3()
  Zq<-makeVar3()
  Cq<-makeVar3()
  Chat<-makeVar2()
  predN<-makeVar2()
  ssb<-matrix(0,nrow=nSpecies,ncol=nYears)
  fq<-1L;     # first season (e.g. quarter)
  lq<-nSeasons
  nlls<-matrix(0,ncol=4,nrow=nSpecies,dimnames=list(species=spNames,nll=c("catch","F","N","survey")))
  
  
  ## Initialize joint negative log likelihood
  ans <- 0
  
  # needed for simulation  and consistencyCheck ? does not work
  # i<-1;
  # for (s in 1:nSpecies) {faf=info[s,'faf']; for (a in nlogFfromTo[s,1]:nlogFfromTo[s,2]) ans<- ans -dnorm(logF[[s]][a,i], sdLogFsta[keyLogFstaSd[s,a+faf-1L]], log=TRUE)}
  # for (s in 1:nSpecies) {for (a in nlogNfromTo[s,1]:nlogNfromTo[s,2]) ans<- ans -dnorm(logN[[s]][a,i]-predN[[s]][a,i],  log=TRUE)}
  

  ##########################################################################################
  ###################  now we begin 
  
  for (s in 1:nSpecies) {
    cat('Species: ',s,'\n')
    ## First take care of F
    fcor <- outer(1:stateDimF[s],
                  1:stateDimF[s],
                  function(i,j)(i!=j)*rho[s] + (i==j))
    fsd <- sdLogFsta[keyLogFstaSd[s,keyLogFstaSd[s,]>0]]
    fvar <- outer(1:stateDimF[s],
                  1:stateDimF[s],
                  function(i,j)fcor[cbind(i,j)]*fsd[i]*fsd[j])
     
    if (zeroCatchYearExistsSp[s]==1)  ans <- ans - sum(dmvnorm( diff( t(logF[[s]][,-zeroCatchYear[[s]]])),mu=0, Sigma=fvar, log=TRUE)) else {
      ans <- ans - sum(dmvnorm( diff( t(logF[[s]])),mu=0, Sigma=fvar, log=TRUE))
    }
    
    if (Debug) nlls[s,"F"]<-  -sum(dmvnorm( diff( t(logF[[s]])) , 0, Sigma=fvar, log=TRUE))
    
   # Spawning Stock Biomass

  # # SSB 1. January !!   
  # for (s in 1:nSpecies) { 
  #   for (y in (1:nYears)) {
  #     ssb[s,y]<-sum(exp(logN[[s]][,y]) *propMat[[s]][y,q,]*stockMeanWeight[[s]][y,q,])  #HUSK AT SÆTTE PropMat for ikke rekruterede aldre til 0
  #   }
  # }  
  
  # SSB 1. January !!
  #  Class 'simref' does not allow operation 'sum' (argument must be fully simulated)
    for (y in 1:nYears) {
      for (a in 1:stateDimN[s]) {
        ssb[s,y]<- ssb[s,y] + exp(logN[[s]][a,y])*propMat[[s]][y,spawnSeason,a]*stockMeanWeight[[s]][y,spawnSeason,a]
      }
    }

    SSB_R<-function(s,y,a=1) {
      if(stockRecruitmentModelCode[s]==0){    ## straight RW
        rec = logN[[s]][a, y-1]
      } 
      else {
        if (stockRecruitmentModelCode[s]==1){ ## Ricker
          rec= rec_loga[s]+log(ssb[s,y-recAge])-exp(rec_logb[s])*ssb[s,y-recAge]
        }else{
          if(stockRecruitmentModelCode[s]==2){  ## B&H
            rec=rec_loga[s]+log(ssb[s,y-recAge])-log(1+exp(rec_logb[s])*ssb[s,y-recAge])
          }else{
            stop("SR model code not recognized");
          }
        }
      }
    }
    
  
    
    ## Now take care of N
    nvar <- outer(1:stateDimN[s], 1:stateDimN[s],
                  function(i,j) (i==j)*varLogN[ keyVarLogN[s,i]])
   
    # kan måske undværes
    #recruits first year given recruitment at age 0 (later in the same year)
    if(stockRecruitmentModelCode[s] >=1 & recAge==0){
      predN[[s]][1,1]<-SSB_R(s,y=1);
      ans <- ans - dmvnorm(logN[[s]][1,1], predN[[s]][1,1], nvar[1,1], log=TRUE) ## N-Process likelihood
    }
    
   for(i in 2:timeSteps) { 
      predN[[s]][1,i]<-SSB_R(s,y=i);
  
      for(j in 2:stateDimN[s]) {
          if (keyLogFsta[s,j-1]>0) {  # to take account for ages with no F
            predN[[s]][j,i]=logN[[s]][j-1,i-1]-exp(logF[[s]][(keyLogFsta[s,j-1]),i-1])-sum(natMor[[s]][i-1,,j-1]) 
          } else { 
            predN[[s]][j,i]=logN[[s]][j-1,i-1]-sum(natMor[[s]][i-1,,j-1]) 
          }
       }
       if(maxAgePlusGroup[s]==1){
          predN[[s]][stateDimN[s],i] = log( exp(logN[[s]][stateDimN[s]-1,i-1]-exp(logF[[s]][keyLogFsta[s,stateDimN[s]-1],i-1])-sum(natMor[[s]][i-1,,stateDimN[s]-1])) +
                                       exp(logN[[s]][stateDimN[s],i-1]  -exp(logF[[s]][keyLogFsta[s,stateDimN[s]],i-1])  -sum(natMor[[s]][i-1,,stateDimN[s]]))  )
       }
      ans <- ans - dmvnorm(logN[[s]][,i], predN[[s]][,i], nvar, log=TRUE) ## N-Process likelihood
      if (Debug) nlls[s,"N"]<- nlls[s,"N"]  - dmvnorm(logN[[s]][,i], predN[[s]][,i], nvar, log=TRUE) 
    }

  #We have now calculated logN (1 quarter)
  

    for(i in 1:timeSteps) {
      for(j in 1:stateDimN[s]) {
        
        #first season (q=fq)
        if (j >1 | recSeason==1 ) {  # first age might not have enter the model yet
          logNq[[s]][i,fq,j]<-logN[[s]][j,i]
          if((keyLogFsta[s,j])>0) {
            Zq[[s]][i,fq,j] <- exp(logF[[s]][(keyLogFsta[s,j]),i])*seasFprop[[s]][i,fq,j]+ natMor[[s]][i,fq,j]
          } else {
            Zq[[s]][i,fq,j] <- natMor[[s]][i,fq,j]
          }
          logNbarq[[s]][i,fq,j] <- logNq[[s]][i,fq,j]-log(Zq[[s]][i,fq,j]) + log(1.0 -exp(-Zq[[s]][i,fq,j]))
          if (keyLogFsta[s,j]>0) Chat[[s]][j,i] <- exp(logNbarq[[s]][i,fq,j]+logF[[s]][keyLogFsta[s,j],i]+log(seasFprop[[s]][i,fq,j]))
        } else Chat[[s]][j,i] = 0.1; #SNYD for at undgå log(0), værdien bruges ikke, så OK
        
        # recruits within the year, need simpler code
        if (j ==1 & recSeason>1) {
           logNq[[s]][i,recSeason,j]<-logN[[s]][j,i] 
           if (keyLogFsta[s,j]>0) {
             Zq[[s]][i,recSeason,j] <- exp(logF[[s]][(keyLogFsta[s,j]),i])*seasFprop[[s]][i,recSeason,j] + natMor[[s]][i,recSeason,j]
           } else {
             Zq[[s]][i,recSeason,j] <- natMor[[s]][i,recSeason,j]
           }
           logNbarq[[s]][i,recSeason,j] <- logNq[[s]][i,recSeason,j]-log(Zq[[s]][i,recSeason,j]) +log(1.0 -exp(-Zq[[s]][i,recSeason,j]))
           if (keyLogFsta[s,j]>0)  {Chat[[s]][j,i] <-  exp(logNbarq[[s]][i,recSeason,j]+logF[[s]][keyLogFsta[s,j],i]+log(seasFprop[[s]][i,recSeason,j]))
           } else { Chat[[s]][j,i]<-0.1  }
        }
        
        #remaining seasons
        if (nSeasons>1) for (q in 2:lq) if (j >1 | q>recSeason ) {
          logNq[[s]][i,q,j]<- logNq[[s]][i,q-1,j]-Zq[[s]][i,q-1,j]
          if (keyLogFsta[s,j]>0) {
            Zq[[s]][i,q,j] <- exp(logF[[s]][(keyLogFsta[s,j]),i])*seasFprop[[s]][i,q,j] + natMor[[s]][i,q,j]
          } else {
            Zq[[s]][i,q,j] <- natMor[[s]][i,q,j]
          }
          logNbarq[[s]][i,q,j] <- logNq[[s]][i,q,j]-log(Zq[[s]][i,q,j]) +log(1.0 -exp(-Zq[[s]][i,q,j]))
          if (keyLogFsta[s,j]>0)  {Chat[[s]][j,i] <- Chat[[s]][j,i]+ exp(logNbarq[[s]][i,q,j]+logF[[s]][keyLogFsta[s,j],i]+log(seasFprop[[s]][i,q,j]))}  else { Chat[[s]][j,i]<-0.1  }
        }
      }
    }

  
  Chat[[s]]<-log(Chat[[s]])
  
  
  ##  match to observations
  
  # first catches
    #cat('Species: ',s,'\n')
    idx<-keyCatch[,"s"]==s
    key.v<-keyCatch[idx,"keyVarObsCatch"]
    yy<-keyCatch[idx,"y"]
    aa<-keyCatch[idx,'a']
    ay<-cbind(aa,yy)
    obs<-logCatchObs[idx]
    predObs<-Chat[[s]][ay]
    for (i in 1:length(predObs)) {
      var <- varLogObsCatch[key.v[i]]
      ans <- ans - dnorm(obs[i],predObs[i],sqrt(var),log=TRUE)
      
      if (Debug==1)  nlls[s,'catch']<- nlls[s,'catch']   - dnorm(obs[i],predObs[i],sqrt(var),log=TRUE)
    }
    # if (Debug==1)  {  # the same and more efficient than above, but problems with convergence ??
    #    var <- varLogObsCatch[key.v]
    #    nlls[s,'catch']<- nlls[s,'catch']   - sum(dnorm(obs,predObs,sqrt(var),log=TRUE))
    #  }

  # and now surveys
  fleets<-keySurvey.overview[keySurvey.overview[,'s']==s,'f']
  surveyType<-keySurvey.overview[,'type']

  for (fl in  fleets)  {
    #cat("survey: ",fl,'\n')
    keys<-keySurvey[keySurvey[,"f"]==fl,]
    q<-keys[1,"q"]
    if (surveyType[fl]==1) {
      for (i in 1:dim(keys)[[1]]) {
        y<-keys[i,'y']
        a<-keys[i,'a']
        keyPowerQ<-keys[i,"keyPowerQ"]
        keyCatchability<-keys[i,"keyCatchability"]
        keyVarObsSurvey<-keys[i,"keyVarObsSurvey"]
        obs.no<-keys[i,'obs.no']
        predObs<-logNq[[s]][y,q,a] - Zq[[s]][y,q,a]*sampleTimeWithinSurvey[fl]
        if(keyPowerQ>0) predObs <- predObs*exp(logQpow[keyPowerQ])
        predObs <- predObs+logCatchability[keyCatchability]
        var <- varLogObsSurvey[keyVarObsSurvey]
        ans <- ans - dnorm(logSurveyObs[obs.no],predObs,sqrt(var),log=TRUE)
        if (Debug==1)  {  
          nlls[s,'survey']<- nlls[s,'survey']   - dnorm(logSurveyObs[obs.no],predObs,sqrt(var),log=TRUE)
        }
       } 
    } else if (surveyType[fl]==2) {  # exploitable biomass (assumed all age with F>0)
      flYears<-keys[,'y']
      faf<-info[s,'faf']; laf<-info[s,'lalike']
      obs.no<-keys[,'obs.no']
      keyCatchability<-keys[1,"keyCatchability"]
      keyVarObsSurvey<-keys[1,"keyVarObsSurvey"]
      # does not work     predObs<-sapply(flYears,function(y) log(sum(exp(logNq[[s]][y,q,faf:laf] - Zq[[s]][y,q,faf:laf]*sampleTimeWithinSurvey[fl])*stockMeanWeight[[s]][y,q,faf:laf]))+ logCatchability[keyCatchability] )
      predObs<-numeric(length(obs.no))
      n<-0L
      for (y in flYears)  {
        n<-n+1L
       predObs[n]<- log(sum(exp(logNq[[s]][y,q,faf:laf] - Zq[[s]][y,q,faf:laf]*sampleTimeWithinSurvey[fl])*stockMeanWeight[[s]][y,q,faf:laf]))+ logCatchability[keyCatchability] 
      }
      
      var <- varLogObsSurvey[keyVarObsSurvey]
      ans <- ans - sum(dnorm(logSurveyObs[obs.no],predObs,sqrt(var),log=TRUE))
      if (Debug==1)  {  
        nlls[s,'survey']<- nlls[s,'survey']   - sum(dnorm(logSurveyObs[obs.no],predObs,sqrt(var),log=TRUE))
      }
    } else if (surveyType[fl]==3) {  # SSB index
      flYears<-keys[,'y']
      faf<-1L; laf<-info[s,'la']
      obs.no<-keys[,'obs.no']
      keyCatchability<-keys[1,"keyCatchability"]
      keyVarObsSurvey<-keys[1,"keyVarObsSurvey"]
      # does not work     predObs<-sapply(flYears,function(y) log(sum(exp(logNq[[s]][y,q,faf:laf] - Zq[[s]][y,q,faf:laf]*sampleTimeWithinSurvey[fl])*stockMeanWeight[[s]][y,q,faf:laf]))+ logCatchability[keyCatchability] )
      predObs<-numeric(length(obs.no))
      n<-0L
      for (y in flYears)  {
        n<-n+1L
        predObs[n]<- log(sum(exp(logNq[[s]][y,q,faf:laf] - Zq[[s]][y,q,faf:laf]*sampleTimeWithinSurvey[fl])*stockMeanWeight[[s]][y,q,faf:laf]*propMat[[s]][y,q,faf:laf]))+ logCatchability[keyCatchability] 
      }
      
      var <- varLogObsSurvey[keyVarObsSurvey]
      ans <- ans - sum(dnorm(logSurveyObs[obs.no],predObs,sqrt(var),log=TRUE))
      if (Debug==1)  {  
        nlls[s,'survey']<- nlls[s,'survey']   - sum(dnorm(logSurveyObs[obs.no],predObs,sqrt(var),log=TRUE))
      } 
    }  else stop(paste("s:",s,"  fleet:",fl,'  fleet type:',surveyType[fl],' is not known\n'))
  } # end fleet loop
  } #end species loop

  ADREPORT(ssb)
  REPORT(logNq)
  REPORT(predN)
  REPORT(Zq)
  REPORT(Chat)
  REPORT(nlls)
  
  ans
}

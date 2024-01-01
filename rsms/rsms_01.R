# rm(list=ls())
library(RTMB)
library(tidyverse)
sam.root<-file.path("~","cod");
rsms.root<-file.path("~","cod","RSMS");

load(file=file.path(rsms.root,"rsms_input.Rdata"),verbose=TRUE)
#str(parameters)
#str(data)

data$stockRecruitmentModelCode[]<-2L
data$Debug<-1L
setwd(rsms.root)

### Parameter transform
# f <- function(x) {2/(1 + exp(-2 * x)) - 1}
square <- function(x){x*x}

func <- function(parameters) {
    getAll(data, parameters,warn=FALSE)
  
  ## Optional (enables extra RTMB features)  
    logCatchObs<-OBS(logCatchObs)
    logSurveyObs<-OBS(logSurveyObs)
    
    ## Optional (enables extra RTMB features)
    #obs<-OBS(obs) #statement tells RTMB that that obs is the response. This is needed to enable automatic simulation and residual calculations from the model object.

    ## Initialize joint negative log likelihood
    ans <- 0

    # first we have all the N states
    logNs<-list()
    for (s in 1:nSpecies) logNs<-c(logNs, list(Un[(nlogNfromTo[s,1]:nlogNfromTo[s,2]),,drop=FALSE]))
  

   # and then the F states
   # with seasonal data, logF is redefined as the sum of seasonal F values (which is not the same as the annual F)

    logFs<-list()
    for (s in 1:nSpecies) logFs<-c(logFs, list(Uf[(nlogFfromTo[s,1]:nlogFfromTo[s,2]),,drop=FALSE]))

    timeSteps <- nYears
    stateDimF <- nlogF
    stateDimN <- nlogN
    sdLogFsta = exp(logSdLogFsta)
    varLogN = exp(logSdLogN*2)
    varLogObsCatch = exp(logSdLogObsCatch*2)
 
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

    ssb<-matrix(0,nrow=nSpecies,ncol=nYears)
    
    fq<-1L;     # first season (e.g. quarter)
    lq<-nSeasons
    nllCatch<-rep(0,nSpecies)

    for (s in (1:nSpecies)) {
      ## First take care of F
      fcor <- outer(1:stateDimF[s],
                    1:stateDimF[s],
                    function(i,j)(i!=j)*rho[s] + (i==j))
      fsd <- sdLogFsta[keyLogFsta[s,keyLogFsta[s,]>0]]
      fvar <- outer(1:stateDimF[s],
                    1:stateDimF[s],
                    function(i,j)fcor[cbind(i,j)]*fsd[i]*fsd[j])

      ans <- ans - sum(dmvnorm( diff( t(logFs[[s]])) , 0, fvar, log=TRUE))
    }
    #cat("ans after diff: ",ans,'\n')
    
    
    ## SIMULATE {
    ##   logF.col(i) = logF.col(i-1) + neg_log_densityF.simulate();
    ## }


    # Spawning Stock Biomass
    if (recSeason > spawnSeason) fa<-2L else fa<-1L    
    q<-spawnSeason
    
 
  # SSB 1. January !!   
  for (s in (1:nSpecies)) { 
     for (y in (1:nYears)) {
        ssb[s,y]<-sum(exp(logNs[[s]][,y]) *propMat[[s]][y,q,]*stockMeanWeight[[s]][y,q,])  #HUSK AT SÆTTE PropMat for ikke rekruterede aldre til 0
     }
  }  
    #for (s in (1:nSpecies)) { 
    #  for (y in (1:nYears)) {   
    # for (a in (fa:info[s,'la'])) {
    #    if (a >= info[s,"faf"]) {
    #       #HUSK proportion of F
    #       ssb[s,y]<-ssb[s,y]+ exp(logNs[[s]][a,y]) * exp(-exp(logFs[[s]][keyLogFsta[s,a] ,y]) * propF[[s]][y,q,a] -natMor[[s]][y,q,a]*propM[[s]][y,q,a])*propMat[[s]][y,q,a]*stockMeanWeight[[s]][y,q,a]
    #     } else {
    #       ssb[s,y]<-ssb[s,y]+ exp(logNs[[s]][a,y]) * exp(                                                          -natMor[[s]][y,q,a]*propM[[s]][y,q,a])*propMat[[s]][y,q,a]*stockMeanWeight[[s]][y,q,a]
    #    }
    # }
  #  }
  # }

   logssb <- log(ssb)
  
 
  for (s in (1:nSpecies)) {
   #cat('Species:',s,'\n')
    ## Now take care of N
    nvar <- outer(1:stateDimN[s], 1:stateDimN[s],
                  function(i,j) (i==j)*varLogN[ keyVarLogN[s,i]])
    predN <- numeric(stateDimN[s])
  
    for(i in 2:timeSteps) {
      if(stockRecruitmentModelCode[s]==0){ ## straight RW
          predN[recAge] = logNs[[s]][recAge, i]
      } else {
          if (stockRecruitmentModelCode[s]==1){ ##ricker
              predN[recAge] = rec_loga[s]+log(ssb[s,i-recAge])-exp(rec_logb[s])*ssb[s,i-recAge]
          }else{
              if(stockRecruitmentModelCode[s]==2){##BH
                  predN[recAge]=rec_loga[2]+log(ssb[s,i-recAge])-log(1+exp(rec_logb[s])*ssb[s,i-recAge])
              }else{
                  stop("SR model code not recognized");
              }
          }
      }

      for(j in 2:stateDimN[s]) {
        if (keyLogFsta[s,j-1]>0) {  # to take account for ages with no F
            predN[j]=logNs[[s]][j-1,i-1]-exp(logFs[[s]][(keyLogFsta[s,j-1]),i-1])-sum(natMor[[s]][i-1,,j-1]) 
        } else { 
          predN[j]=logNs[[s]][j-1,i-1]-sum(natMor[[s]][i-1,,j-1]) 
        }
      }
      if(maxAgePlusGroup[s]==1){
        predN[stateDimN[s]] = log( exp(exp(logFs[[s]][keyLogFsta[s,stateDimN[s]-1],i-1])-sum(natMor[[s]][i-1,,stateDimN[s]-1])) +
                                 exp(logNs[[s]][stateDimN[s],i-1]-exp(logFs[[s]][keyLogFsta[s,stateDimN[s]],i-1])-sum(natMor[[s]][i-1,,stateDimN[s]]))  )
      }
      
     #  cat(s,i,j,'\n')
     #  print(class(predN))
     # print(class(logNs[[s]][,i]))
      ans <- ans - dmvnorm(logNs[[s]][,i], predN, nvar, log=TRUE) ## N-Process likelihood
    }
  } #end species loop
  
  #cat("ans after N-Process likelihood: ",ans,'\n')
  
  for (s in (1:nSpecies)) {
    for(i in 1:timeSteps) {
      for(j in 1:stateDimN[s]) {

        #first season
        if (j >1 | recSeason==1 ) {  # first age might not have reitrd to the model
          logNq[[s]][i,fq,j]<-logNs[[s]][j,i]
          if((keyLogFsta[s,j])>0) {
              Zq[[s]][i,fq,j] <- exp(logFs[[s]][(keyLogFsta[s,j]),i])*seasFprop[[s]][i,fq,j]+ natMor[[s]][i,fq,j]
          } else {
            Zq[[s]][i,fq,j] <- natMor[[s]][i,fq,j]
          }
          logNbarq[[s]][i,fq,j] <- logNq[[s]][i,fq,j]-log(Zq[[s]][i,fq,j]) +log(1.0 -exp(-Zq[[s]][i,fq,j]))
          if (keyLogFsta[s,j]>0) Chat[[s]][j,i] <- exp(logNbarq[[s]][i,fq,j]+logFs[[s]][keyLogFsta[s,j],i]+log(seasFprop[[s]][i,fq,j]))
        } else Chat[[s]][j,i] = 0.1; #SNYD

        #remaining seasons
        if (nSeasons>1) for (q in 2:lq) if (j >1 | q>=recSeason ) {
          logNq[[s]][i,q,j]<- logNq[[s]][i,q-1,j]-Zq[[s]][i,q-1,j]
          if (keyLogFsta[s,j]>0) {
             Zq[[s]][i,q,j] <- exp(logFs[[s]][(keyLogFsta[s,j]),i])*seasFprop[[s]][i,q,j] + natMor[[s]][i,q,j]
          } else {
            Zq[[s]][i,q,j] <- natMor[[s]][i,q,j]
          }
          logNbarq[[s]][i,q,j] <- logNq[[s]][i,q,j]-log(Zq[[s]][i,q,j]) +log(1.0 -exp(-Zq[[s]][i,q,j]))
          if (keyLogFsta[s,j]>0)  {Chat[[s]][j,i] <- Chat[[s]][j,i]+ exp(logNbarq[[s]][i,q,j]+logFs[[s]][keyLogFsta[s,j],i]+log(seasFprop[[s]][i,q,j]))}
        }
      }
     }
  } # end species loop
  
  Chat<-lapply(Chat,log)


  ##  match to observations
 # first catches

  # first catches
  for (s in 1:nSpecies ){
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
       #cat(yy[i]-off.year,aa[i],obs[i],predObs[i],sqrt(var),dnorm(obs[i],predObs[i],sqrt(var),log=TRUE),'\n')
       ans <- ans - dnorm(obs[i],predObs[i],sqrt(var),log=TRUE)
    }
    if (Debug==1)  {  # the same and more efficient than above, but problems with convergence 
      var <- varLogObsCatch[key.v]
      nllCatch[s]<-   - sum(dnorm(obs,predObs,sqrt(var),log=TRUE))
    }
  }
  
  
# and now surveys

nFleets<-max(keySurvey[,"f"])
for (fl in  1:nFleets) {
  #cat("survey: ",fl,'\n')
  keys<-keySurvey[keySurvey[,"f"]==fl,]
  q<-keys[1,"q"]
  s<-keys[1,'s']
  for (i in 1:dim(keys)[[1]]) {
      y<-keys[i,'y']
      a<-keys[i,'a']
      keyPowerQ<-keys[i,"keyPowerQ"]
      keyCatchability<-keys[i,"keyCatchability"]
      keyVarObsSurvey<-keys[i,"keyVarObsSurvey"]
      obs.no<-keys[i,'obs.no']
      predObs<-logNq[[s]][y,q,a] - Zq[[s]][y,q,a]*sampleTimeWithinSurvey[fl]
      if(keyPowerQ>0) predObs <- predObs*exp(logQpow[keyPowerQ])
      predObs <- predObs+logCatachability[keyCatchability]

      var <- varLogObsSurvey[keyVarObsSurvey]
      #cat(i,s,fl,y,y-off.year,q,a,logSurveyObs[obs.no],predObs," var=",var,keyCatchability,keyVarObsSurvey,'\n')
      ans <- ans - dnorm(logSurveyObs[obs.no],predObs,sqrt(var),log=TRUE)
  }
} # end fleet loop
 #cat("ans: after survey ",ans,'\n')
 
# keySurvey.overview 
 
        ## SIMULATE {
        ##   obs(i,3) = exp( rnorm(predObs, sqrt(var)) ) ;
        ## }
  
   ADREPORT(ssb)
    REPORT(logNq)
    REPORT(Zq)
    REPORT(Chat)
    REPORT(nllCatch)
    
    ans
}

#Going back to our objective function f, first step is to check that you can evaluate the function as a normal R function:
# func(parameters)   # KALDET VIL PÅVIRKE KALDET TIL MakeAdFun !!!?
#An error at this point is obviously not due to RTMB.

#Next, it is useful to check that MakeADFun can be run without random effects:

# obj <- MakeADFun(func, parameters)
obj <- MakeADFun(func, parameters, random=c("Un","Uf"),silent=FALSE)

lower <- obj$par*0-Inf
upper <- obj$par*0+Inf
lower["rho"] <- 0.01
upper["rho"] <- 0.99

opt <- nlminb(obj$par, obj$fn, obj$gr, lower=lower, upper=upper)

cat("\nobjective:",opt$objective,"  convergence:",opt$convergence) # 0 indicates successful convergence.


rep<-obj$report()
str(rep,1)
lapply(rep$logNq,function(x) round(exp(x[,1,])))
lapply(rep$Zq,function(x) round(exp(x[,1,]),2))
lapply(rep$Chat,function(x) round(exp(x)))
rep$nlls

rep <- sdreport(obj)
x<-as.list(rep, "Est")
str(x)
round(exp(x$Uf),3)
round(exp(x$Un),0)

#############################
#Once obj has been constructed successfully, you should evaluate it

if (FALSE) obj$fn()

#and verify that it gives the same result as func(parameters).

if (FALSE) {  # give errors (as also seen in the SAM script from Anders)
  set.seed(1)
  chk <- checkConsistency(obj)
  chk
}

if (FALSE) {
  debug(func)
  obj <- MakeADFun(func, parameters, random=c("Un","Uf"),silent=FALSE)
}

if (FALSE) {
  set.seed(1)
  chk <- checkConsistency(obj)
  chk
}

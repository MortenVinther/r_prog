#rm(list=ls())
library(RTMB)
library(tidyverse)
sam.root<-file.path("~","cod");

rsms.root<-file.path("~","cod","RSMS");

load(file=file.path(rsms.root,"rsms_input.Rdata"),verbose=TRUE)
str(data,1)
setwd(rsms.root)

if (TRUE) {
  attach(data, warn.conflicts = FALSE)       ## discouraged but
  attach(parameters, warn.conflicts = FALSE) ## VERY handy while developing!
}
#str(data,1)

### Parameter transform
f <- function(x) {2/(1 + exp(-2 * x)) - 1}
square <- function(x){x*x}

func <- function(
                 logFpar,
                 logQpow,
                 logSdLogFsta,
                 logSdLogN,
                 logSdLogObs,
                 rec_loga,
                 rec_logb,
                 rho,
                 logScale,
                 logScaleSSB,
                 logPowSSB,
                 logSdSSB,
                 U) {

  #getAll(data, warn=FALSE)

    ## Optional (enables extra RTMB features)
    #obs<-OBS(obs) #statement tells RTMB that that obs is the response. This is needed to enable automatic simulation and residual calculations from the model object.

    ## Initialize joint negative log likelihood
    ans <- 0

    # first we have all the N states
    logNs<-list()
    for (s in 1:nSpecies) logNs<-c(logNs, list(U[(nlogNfromTo[s,1]:nlogNfromTo[s,2]),,drop=FALSE]))


   # and then the F states
   # with seasonal data, logF is redefined as the sum of seasonal F values (which is not the same as the annual F)

    logFs<-list()
    for (s in 1:nSpecies) logFs<-c(logFs, list(U[(nlogFfromTo[s,1]:nlogFfromTo[s,2])+ sum(nlogN),,drop=FALSE]))

    timeSteps <- nYears
    stateDimF <- nlogF
    stateDimN <- nlogN
    sdLogFsta = exp(logSdLogFsta)
    varLogN = exp(logSdLogN*2)
    varLogObs = exp(logSdLogObs*2)

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

    fq<-1L;     # first season (e.g. quarter)
    lq<-nSeasons



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
    ## SIMULATE {
    ##   logF.col(i) = logF.col(i-1) + neg_log_densityF.simulate();
    ## }


    # Spawning Stock Biomass
    if (recSeason > spawnSeason) fa<-2 else fa<-1
    q<-spawnSeason
    ssb<-matrix(0,ncol=nYears,nrow=nSpecies)
    for (s in (1:nSpecies)) {
      for (y in (1:nYears)) {
        for (a in (fa:info[s,'la'])) {
            if (a >= info[s,"faf"]) {
              #HUSK proportion of F
                ssb[s,y]<-ssb[s,y] + exp(logNs[[s]][a,y]) * exp(-exp(logFs[[s]][keyLogFsta[s,a] ,y]) * propF[[s]][y,q,a] -propM[[s]][y,q,a]) * propMat[[s]][y,q,a]*stockMeanWeight[[s]][y,q,a]
            } else {
                ssb[s,y]<-ssb[s,y]+ exp(logNs[[s]][a,y]) * exp(-natMor[[s]][y,q,a]*propM[[s]][y,q,a]) * propMat[[s]][y,q,a]*stockMeanWeight[[s]][y,q,a]
            }
        }
      }
    }
   logssb <- log(ssb)


  for (s in (1:nSpecies)) {
   cat('Species:',s,'\n')
    ## Now take care of N
    nvar <- outer(1:stateDimN[s], 1:stateDimN[s],
                  function(i,j) (i==j)*varLogN[ keyVarLogN[s,i]])
    predN <- numeric(stateDimN[s])

    for(i in 2:timeSteps) {
      if(stockRecruitmentModelCode[s]==0){ ## straight RW
          predN[reca] = logNs[[s]][reca, i]
      } else {
          if (stockRecruitmentModelCode[s]==1){ ##ricker
              predN[reca] = rec_loga[s]+log(ssb[s,i-reca])-exp(rec_logb[s])*ssb[s,i-reca]
          }else{
              if(stockRecruitmentModelCode[s]==2){##BH
                  predN[reca]=rec_loga[2]+log(ssb[s,i-reca])-log(1+exp(rec_logb[s])*ssb[s,i-reca])
              }else{
                  stop("SR model code not recognized");
              }
          }
      }

      for(j in 2:stateDimN[s]) {
        if (keyLogFsta[s,j-1]>0) tmpF<-logFs[[s]][(keyLogFsta[s,j-1]),i-1] else tmpF<-0 # to take account for ages with no F
        predN[j]=logNs[[s]][j-1,i-1]-exp(tmpF)-sum(natMor[[s]][i-1,,j-1])
      }
      if(maxAgePlusGroup[s]==1){
        predN[stateDimN[s]] = log( exp(exp(logFs[[s]][keyLogFsta[s,stateDimN[s]-1],i-1])-sum(natMor[[s]][i-1,,stateDimN[s]-1])) +
                                 exp(logNs[[s]][stateDimN[s],i-1]-exp(logFs[[s]][keyLogFsta[s,stateDimN[s]],i-1])-sum(natMor[[s]][i-1,,stateDimN[s]]))  )
      }
      ## SIMULATE {
          ##   logN.col(i) = predN + neg_log_densityN.simulate();
      ## }

      ans <- ans - dmvnorm(logNs[[s]][,i], predN, nvar, log=TRUE) ## N-Process likelihood
    }

     for(i in 1:timeSteps) {
      for(j in 1:stateDimN[s])

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
        for (q in 2:lq) if (j >1 | q>=recSeason ) {
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
  keyVarObsCatch # not used
debug<-FALSE
ans=0;
  # first catches
  for (s in 1:nSpecies ){
    cat('Species: ',s,'\n')
    idx<-keyCatch[,"s"]==s
    key.v<-keyCatch[idx,"keyVarObsCatch"]
    yy<-keyCatch[idx,"y"]
    aa<-keyCatch[idx,'a']
    ay<-cbind(aa,yy)
    obs<-logCatchObs[idx]
    predObs<-as.vector(Chat[[s]][ay])
    if (debug) for (i in 1:length(predObs)) {
       var <- VarObsCatch[key.v[i]]
       #cat(yy[i]-off.year,aa[i],obs[i],predObs[i],sqrt(var),dnorm(obs[i],predObs[i],sqrt(var),log=TRUE),'\n')
       ans <- ans - dnorm(obs[i],predObs[i],sqrt(var),log=TRUE)
    } else {
      var <- VarObsCatch[key.v]
      ans <- ans - sum(dnorm(obs,predObs,sqrt(var),log=TRUE))
    }
cat(ans,'\n')
  }

ans
    head(keyCatch)
    logCatchObs,
    nCatchObs

    keyFleet
    logSurveyObs


    minYear <- keys[1,1]
    predObs <- 0
    var <- 0
    for(i in 1:nobs){
        y <- keys[i,1] - minYear + 1 ## 1 based
        f <- keys[i,2] ## 1 based
        ft <- fleetTypes[f] ## 1 based
        a <- keys[i,3] - minAge + 1 ## 1 based
        if(ft==0) {   ## residual fleet
           predObs<-Chat[y,a]
        }else{
            if(ft==1){## comm fleet
                stop("Not implemented yet!!!")
            }else{
                if(ft==2){## survey
                    q<-sampleTimesSeason[f]
                    predObs<-logNq[y,q,a] - Zq[y,q,a]*sampleTimeWithin[f]
                   if(keyQpow[f,a]>(-1)){
                        predObs = predObs*exp(logQpow[keyQpow[f,a]+1L])
                   }
                    if(keyLogFpar[f,a]>(-1)){
                        predObs <- predObs+logFpar[keyLogFpar[f,a]+1L]
                    }
                }else{
                    if(ft==3){## SSB survey -- nevermind for now
                        stop("Not implemented yet!!!")
                    }else{
                        if(ft==4){## SSB survey -- nevermind for now
                            stop("Not implemented yet!!!")
                        }
                    }
                }
            }
        }
        var <- varLogObs[keyVarObs[f,a]+1L]
        ans <- ans - dnorm(obs[i],predObs,sqrt(var),log=TRUE)
        ## SIMULATE {
        ##   obs(i,3) = exp( rnorm(predObs, sqrt(var)) ) ;
        ## }
    }
  } # end species loop
    ADREPORT(ssb)
    REPORT(logNq)
    REPORT(Zq)
    REPORT(Chat)
    #REPORT(predNy)
    ans
}

## Test that we can evaluate using numeric types
environment(func) <- list2env(data)
do.call(func, parameters)

obj <- MakeADFun(function(p)do.call(func,p), parameters, random=c("U"), DLL="sam",silent=TRUE)

lower <- obj$par*0-Inf
upper <- obj$par*0+Inf
lower["rho"] <- 0.01
upper["rho"] <- 0.99

system.time(opt <- nlminb(obj$par, obj$fn, obj$gr, lower=lower, upper=upper))

cat("\nobjective:",opt$objective,"  convergence:",opt$convergence) # 0 indicates successful convergence.

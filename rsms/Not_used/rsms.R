#rm(list=ls())

load(file=file.path(rsms.root,"rsms_input.Rdata"),verbose=TRUE)
#str(parameters)
#str(data)

data$stockRecruitmentModelCode[]<-1L

data$Debug<-1L

data$doSpecies
#data$doSpecies<-c(1L,3L,7L)

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
    ssb<-matrix(0,nrow=nSpecies,ncol=nYears)
    fq<-1L;     # first season (e.g. quarter)
    lq<-nSeasons
    nlls<-matrix(0,ncol=4,nrow=nSpecies,dimnames=list(species=spNames,nll=c("catch","F","N","survey")))
                 

    for (s in (doSpecies)) {
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
    }
     
    
    ## SIMULATE {
    ##   logF.col(i) = logF.col(i-1) + neg_log_densityF.simulate();
    ## }


    # Spawning Stock Biomass
    q<-spawnSeason
    
 
  # SSB 1. January !!   
  for (s in (doSpecies)) { 
     for (y in (1:nYears)) {
        ssb[s,y]<-sum(exp(logN[[s]][,y]) *propMat[[s]][y,q,]*stockMeanWeight[[s]][y,q,])  #HUSK AT SÆTTE PropMat for ikke rekruterede aldre til 0
     }
  }  
 
  logssb <- log(ssb)
  
  for (s in (doSpecies)) {
    ## Now take care of N
    nvar <- outer(1:stateDimN[s], 1:stateDimN[s],
                  function(i,j) (i==j)*varLogN[ keyVarLogN[s,i]])
    predN <- numeric(stateDimN[s])
  
    for(i in 2:timeSteps) {
      if(stockRecruitmentModelCode[s]==0){ ## straight RW
          predN[recAge] = logN[[s]][recAge, i]
      } else {
          if (stockRecruitmentModelCode[s]==1){ ## Ricker
              predN[recAge] = rec_loga[s]+log(ssb[s,i-recAge])-exp(rec_logb[s])*ssb[s,i-recAge]
          }else{
              if(stockRecruitmentModelCode[s]==2){  ## B&H
                  predN[recAge]=rec_loga[s]+log(ssb[s,i-recAge])-log(1+exp(rec_logb[s])*ssb[s,i-recAge])
              }else{
                  stop("SR model code not recognized");
              }
          }

      }

      for(j in 2:stateDimN[s]) {
        if (keyLogFsta[s,j-1]>0) {  # to take account for ages with no F
            predN[j]=logN[[s]][j-1,i-1]-exp(logF[[s]][(keyLogFsta[s,j-1]),i-1])-sum(natMor[[s]][i-1,,j-1]) 
        } else { 
          predN[j]=logN[[s]][j-1,i-1]-sum(natMor[[s]][i-1,,j-1]) 
        }
      }
      if(maxAgePlusGroup[s]==1){
        predN[stateDimN[s]] = log( exp(logN[[s]][stateDimN[s]-1,i-1]-exp(logF[[s]][keyLogFsta[s,stateDimN[s]-1],i-1])-sum(natMor[[s]][i-1,,stateDimN[s]-1])) +
                                   exp(logN[[s]][stateDimN[s],i-1]  -exp(logF[[s]][keyLogFsta[s,stateDimN[s]],i-1])  -sum(natMor[[s]][i-1,,stateDimN[s]]))  )
      }
      ans <- ans - dmvnorm(logN[[s]][,i], predN, nvar, log=TRUE) ## N-Process likelihood
      if (Debug) nlls[s,"N"]<- nlls[s,"N"]  - dmvnorm(logN[[s]][,i], predN, nvar, log=TRUE) 
    }
  } #end species loop

  for (s in (doSpecies)) {
    for(i in 1:timeSteps) {
      for(j in 1:stateDimN[s]) {

        #first season
        if (j >1 | recSeason==1 ) {  # first age might not have enter the model yet
          logNq[[s]][i,fq,j]<-logN[[s]][j,i]
          if((keyLogFsta[s,j])>0) {
              Zq[[s]][i,fq,j] <- exp(logF[[s]][(keyLogFsta[s,j]),i])*seasFprop[[s]][i,fq,j]+ natMor[[s]][i,fq,j]
          } else {
            Zq[[s]][i,fq,j] <- natMor[[s]][i,fq,j]
          }
          logNbarq[[s]][i,fq,j] <- logNq[[s]][i,fq,j]-log(Zq[[s]][i,fq,j]) +log(1.0 -exp(-Zq[[s]][i,fq,j]))
          if (keyLogFsta[s,j]>0) Chat[[s]][j,i] <- exp(logNbarq[[s]][i,fq,j]+logF[[s]][keyLogFsta[s,j],i]+log(seasFprop[[s]][i,fq,j]))
        } else Chat[[s]][j,i] = 0.1; #SNYD for at undgå log(0), værdien bruges ikke, så OK
      
        #remaining seasons
        if (j ==1 & recSeason>1) logNq[[s]][i,recSeason,j]<-logN[[s]][j,i] 
        
        if (nSeasons>1) for (q in 2:lq) if (j >1 | q>recSeason ) {
          logNq[[s]][i,q,j]<- logNq[[s]][i,q-1,j]-Zq[[s]][i,q-1,j]
          if (keyLogFsta[s,j]>0) {
             Zq[[s]][i,q,j] <- exp(logF[[s]][(keyLogFsta[s,j]),i])*seasFprop[[s]][i,q,j] + natMor[[s]][i,q,j]
          } else {
            Zq[[s]][i,q,j] <- natMor[[s]][i,q,j]
          }
          logNbarq[[s]][i,q,j] <- logNq[[s]][i,q,j]-log(Zq[[s]][i,q,j]) +log(1.0 -exp(-Zq[[s]][i,q,j]))
          if (keyLogFsta[s,j]>0)  {Chat[[s]][j,i] <- Chat[[s]][j,i]+ exp(logNbarq[[s]][i,q,j]+logF[[s]][keyLogFsta[s,j],i]+log(seasFprop[[s]][i,q,j]))}
        }
      }
     }
  } # end species loop
  
  Chat<-lapply(Chat,log)


  ##  match to observations

  # first catches
  for (s in doSpecies ){
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
      predObs <- predObs+logCatchability[keyCatchability]

      var <- varLogObsSurvey[keyVarObsSurvey]
      ans <- ans - dnorm(logSurveyObs[obs.no],predObs,sqrt(var),log=TRUE)
      
      if (Debug==1)  {  
        nlls[s,'survey']<- nlls[s,'survey']   - dnorm(logSurveyObs[obs.no],predObs,sqrt(var),log=TRUE)
      }
  }
} # end fleet loop

# keySurvey.overview 
 
        ## SIMULATE {
        ##   obs(i,3) = exp( rnorm(predObs, sqrt(var)) ) ;
        ## }
  
    ADREPORT(ssb)
    REPORT(logNq)
    REPORT(Zq)
    REPORT(Chat)
    REPORT(nlls)
    
    ans
}

#Going back to our objective function f, first step is to check that you can evaluate the function as a normal R function:
# func(parameters)   # KALDET VIL PÅVIRKE KALDET TIL MakeAdFun !!!?
#An error at this point is obviously not due to RTMB.

# adjust if the are species/year combination with zero catches (or assumed F is very low and highly uncertain)
if (data$zeroCatchYearExists==1) {
  UfMap<-matrix(1L:(dim(parameters$Uf)[[1]]*dim(parameters$Uf)[[2]]),nrow=sum(data$nlogF),byrow=TRUE)
  for (s in 1:data$nSpecies) if (length(data$zeroCatchYear[[s]]) >0 ) {
    zy<-data$zeroCatchYear[[s]]
    fromTo<-data$nlogFfromTo[s,]
    UfMap[fromTo[1]:fromTo[2],zy]<-NA
    parameters$Uf[fromTo[1]:fromTo[2],zy]<-log(0.001)
  }
  UfMap<-factor(UfMap)
}

# fix parameters not used, due to subset of species
doFix<-(data$nSpecies > length(unique(data$doSpecies)))
if (doFix) { 
 ms<- setdiff(1:data$nSpecies,unique(data$doSpecies))
  excl_par<-function(key,pars) {
    k<-key[ms,]
    k<-k[k>0]
    map<- as.numeric(1:length(pars))
    map[k]<-NA
    factor(map)
  }
  
  excl_par2<-function(pars) {
    map<- as.numeric(1:length(pars))
    map[ms]<-NA
    factor(map)
  }
  
  excl_par3<-function(pars,type=1) {
    map<-matrix(1:(ncol(pars)*nrow(pars)),ncol=ncol(pars),nrow=nrow(pars))
    if (type==1) map[data$nlogFfromTo[ms,],]<-NA else map[data$nlogNfromTo[ms,],]<-NA
    factor(map)
  }

#  str(parameters) 
  
  stopifnot(length(parameters$logQpow)==0) # tmp fix
  
  mapFix<-list(
    Uf=excl_par3(parameters$Uf,type=1),
    Un=excl_par3(parameters$Un,type=2),
    
    logSdLogObsCatch  =excl_par(key=data$keyVarObsCatch,  pars=parameters$logSdLogObsCatch),
    logCatchability   =excl_par(key=data$keyCatchability, pars=parameters$logCatchability),
    logSdLogObsSurvey =excl_par(key=data$keyVarObsSurvey, pars=parameters$logSdLogObsSurvey),
    logSdLogFsta      =excl_par(key=data$keyLogFstaSd,    pars=parameters$logSdLogFsta),
    logSdLogN         =excl_par(key=data$keyVarLogN,      pars=parameters$logSdLogN),
    
    rho=excl_par2(parameters$rho),
    rec_loga=excl_par2(parameters$rec_loga),
    rec_logb=excl_par2(parameters$rec_logb)
  )
 str(mapFix)  
}

if (doFix) my.map<-mapFix else my.map=list() 
if (data$zeroCatchYearExists==1) my.map<-list(Uf=UfMap) else my.map=list()


  #logSdLogObsSurvey=factor(rep(NA,length(parameters$logSdLogObsSurvey))),
  # logSdLogN =factor(rep(NA,length(parameters$logSdLogN)))
  #rho =factor(rep(NA,length(parameters$rho)))



obj <- MakeADFun(func, parameters, random=c("Un","Uf"),silent=FALSE,map=my.map)

lower <- obj$par*0-Inf
upper <- obj$par*0+Inf
lower["rho"] <- 0.01
upper["rho"] <- 0.99

nl<-names(lower)

lower[nl=="logSdLogObsSurvey"]<-rep(log(0.15),length(parameters$logSdLogObsSurvey))
upper[nl=="logSdLogObsSurvey"]<-rep(log(2.0),length(parameters$logSdLogObsSurvey))

lower[nl=="logSdLogObsCatch"]<-rep(log(0.1),length(parameters$logSdLogObsCatch))
upper[nl=="logSdLogObsCatch"]<-rep(log(2.0),length(parameters$logSdLogObsCatch))

#t(rbind(lower,upper))    
opt <- nlminb(obj$par, obj$fn, obj$gr, lower=lower, upper=upper)

cat("\nobjective:",opt$objective,"  convergence:",opt$convergence) # 0 indicates successful convergence.

if (!Batch) {
  
rep<-obj$report()

#lapply(rep$logNq,function(x) (exp(x[,1,])))
#if (data$nSeasons==4)lapply(rep$logNq,function(x) (exp(x[,3,])))

#lapply(rep$Zq,function(x) round(exp(x[,1,]),2))
#lapply(rep$Chat,function(x) round(exp(x)))

cbind(rep$nlls,all=rowSums(rep$nlls))

Est<-as.list(rep, "Est", report=TRUE)

sdrep <- sdreport(obj)
obj$fn()  #  value er den samme som sum af min nnls
x<-as.list(sdrep, "Est")

round(exp(x$Uf),3)
round(exp(x$Un),0)

make_tab<-function(d,key,printIt=TRUE,roundIt=2) {
  d<-exp(d)
  tab<-key
  for (i in 1:length(d)) tab[tab[,]==i]<- d[i]
  tab
  if (printIt) {
    ptab<-tab
    ptab[ptab== -1] <-NA
    print(round(ptab,roundIt),na.print='.')
  } 
  invisible(tab)
}


make_tab(d=x$logSdLogN,key=data$keyVarLogN,roundIt=2)
make_tab(d=x$logSdLogFsta,key=data$keyLogFstaSd,roundIt=2)
make_tab(d=x$logSdLogObsCatch,key=data$keyVarObsCatch)  
make_tab(d=x$logSdLogObsSurvey,key=data$keyVarObsSurvey) 
make_tab(d=x$logCatchability,key=data$keyCatchability) 


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

cat("\nobjective:",opt$objective,"  convergence:",opt$convergence, "  # 0 indicates successful convergence\n")
}


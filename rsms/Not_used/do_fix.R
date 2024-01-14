
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


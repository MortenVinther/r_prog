lowerUpper<-function(obj,data,parameters) {
    nl<-names(obj$par)
  
  spLimit<-function(nl,key,param,l=20) {
    a<-array2DF(key) %>% mutate(Var2=NULL) %>% filter(Value>=1) %>%unique() %>% arrange(Value)
    found<-nl==param
    nl[found]<-paste(nl[found],substr(a$Var1,1,l),sep='_')
    nl
  }
  
  
  nlSp<-nl
  nlSp<-spLimit(nlSp,key=data$keyVarObsCatch,param="logSdLogObsCatch") 
  nlSp<-spLimit(nlSp,key=data$keyCatchability,param="logCatchability",l=3) 
  nlSp<-spLimit(nlSp,key=data$keyVarObsSurvey,param="logSdLogObsSurvey",l=3) 
  nlSp<-spLimit(nlSp,key=data$keyLogFstaSd,param="logSdLogFsta") 
  nlSp<-spLimit(nlSp,key=data$keyVarLogN,param="logSdLogN") 
  
  nlSp  
  
  found<-nlSp=='rho';      nlSp[found]<- paste(nlSp[found],data$spNames[data$useRho],sep='_')
  
  found<-nlSp=='rec_loga'; nlSp[found]<- paste(nlSp[found],data$spNames[data$stockRecruitmentModelCode>0],sep='_')
  
  found<-nlSp=='rec_logb'; nlSp[found]<- paste(nlSp[found],data$spNames[data$stockRecruitmentModelCode>0 & data$stockRecruitmentModelCode!=3],sep='_')
  
  
  
  lower <- obj$par*0-Inf
  upper <- obj$par*0+Inf
  lower[nl=='rho'] <- 0.01
  upper[nl=='rho'] <- 0.99
  
  lower[nl=="logSdLogObsSurvey"]<-rep(log(data$minVarObsSurvey),length(parameters$logSdLogObsSurvey))
  upper[nl=="logSdLogObsSurvey"]<-rep(log(2.0),length(parameters$logSdLogObsSurvey))
  
  lower[nl=="logSdLogObsCatch"]<-rep(log(data$minVarObsCatch),length(parameters$logSdLogObsCatch))
  upper[nl=="logSdLogObsCatch"]<-rep(log(2.0),length(parameters$logSdLogObsCatch))
  
  
  
  # N.sandeel
  #upper[nl=="logSdLogN"]<-log(c(20,0.1))
  
  
  # mac
  #upper[nl=="logSdLogN"]<-log(c(0.10))
  return(list(lower=lower,upper=upper))
}

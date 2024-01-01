rm(list=ls())
library(tidyverse)
sam.root<-file.path("~","cod");
rsms.root<-file.path("~","cod","RSMS");

if (annual) {
  
  load(file=file.path(rsms.root,"rsms_input.Rdata"),verbose=TRUE)
  
  attach(data)
  
  data$nSeasons<-1L
  data$recSeason<-1L
  
  data$propMat <-lapply(propMat,function(x) x[,1,,drop=FALSE])
  data$stockMeanWeight<- lapply(stockMeanWeight,function(x) {
    x<-apply(x,c(1,3),mean) 
    x<-array2DF(x) %>% mutate(q=1) %>% mutate_if(is.character,as.integer)
    x<-tapply(x$Value,list(x$Var1,x$q,x$Var2),sum)
  })
  data$catchMeanWeight<- lapply(catchMeanWeight,function(x) x[,1,,drop=FALSE])  # to be changed
  data$seasFprop<- lapply(seasFprop,function(x) x[,1,,drop=FALSE])  
  data$propF<- lapply(propF,function(x) x[,1,,drop=FALSE])   
  data$propM<- lapply(propM,function(x) x[,1,,drop=FALSE])
  data$natMor<- lapply(natMor,function(x) {
    x<-apply(x,c(1,3),sum) 
    x<-array2DF(x) %>% mutate(q=1) %>% mutate_if(is.character,as.integer)
    x<-tapply(x$Value,list(x$Var1,x$q,x$Var2),sum)
  })
  data$seasFprop<- lapply(seasFprop,function(x) {
    x<-apply(x,c(1,3),sum) 
    x<-array2DF(x) %>% mutate(q=1) %>% mutate_if(is.character,as.integer)
    x<-tapply(x$Value,list(x$Var1,x$q,x$Var2),sum)
  })
  
  x<-cbind(sampleTimeWithinSurvey,q=keySurvey.overview[,'q'])
  x[,"sampleTimeWithin"]<-x[,"sampleTimeWithin"]/4+(x[,"q"]-1)*0.25
  data$sampleTimeWithinSurvey<-x[,1]
  detach(data)
  
  data$keySurvey[,'q']<-1L
  
  save(data,parameters,file=file.path(rsms.root,"rsms_input.Rdata"))
}


load(file=file.path(rsms.root,"rsms_input.Rdata"),verbose=TRUE)
#str(parameters)
#str(data)

ps<-c(2L,9L)

p<-parameters
d<-data

d$info<-data$info[ps,]
d$nSpecies<-length(ps)
d$doSpecies <-ps

d$nAges<-max(d$info[,'last-age'])-data$minAge+1 
ages<-1:d$nAges
d$spNames=rownames(d$info)                
d$fbarRange<-data$fbarRange[ps,]     

d$stockRecruitmentModelCode<-data$stockRecruitmentModelCode[ps]
#$ zeroCatchYearExists      : int 1
d$zeroCatchYearExistsSp<-data$zeroCatchYearExistsSp[ps]   
d$zeroCatchYear<-data$zeroCatchYear[ps]

cut_tab<-function(tab,reNumber=FALSE) {
 tab<-tab[ps,ages]
 if (reNumber) {
   k<-sort(unique(tab[tab>0]))
   kk<-1:length(k)
   names(kk)<-k
   tab[tab>0]<-kk[as.character(tab[tab>0])]
 }
 tab
}
cutFromTo<-function(tab){
  tab<-tab[ps,,drop=FALSE]
  n<-tab[,2]-tab[,1]+1L
  tab[,2]<-cumsum(n)
  tab[,1]<-tab[,2]-n+1L
  tab
}

d$keyLogFsta<-cut_tab(data$keyLogFsta,reNumber=FALSE)          
d$keyLogFstaSd<-cut_tab(data$keyLogFstaSd,reNumber=TRUE)
d$nlogF <-max(d$keyLogFstaSd)
tab<-d$nlogFfromTo; 
d$nlogFfromTo<-cutFromTo(data$nlogFfromTo)
d$nlogNfromTo<-cutFromTo(data$nlogNfromTo)
d$nlogN<-max(d$nlogNfromTo) 
d$keyVarLogN<-cut_tab(data$keyVarLogN,reNumber=TRUE)
d$keyVarObsCatch<-cut_tab(data$keyVarObsCatch,reNumber=TRUE)
d$keyVarObsSurvey<-cut_tab(data$keyVarObsSurvey,reNumber=TRUE)
d$propMat<-data$propMat[ps]
d$stockMeanWeight<-data$stockMeanWeight[ps]
d$catchMeanWeight<-data$catchMeanWeight[ps]
d$seasFprop<-data$seasFprop[ps]
d$natMor<-data$natMor[ps]
d$propF<-data$propF[ps]
d$propM<-data$propM[ps]

d$keyCatch<-data$keyCatch[data$keyCatch[,'s'] %in% ps,]
  k<-sort(unique(d$keyCatch[,'keyVarObsCatch']))
  kk<-1:length(k)
  names(kk)<-k
  d$keyCatch[,'keyVarObsCatch']<-kk[as.character(d$keyCatch[,'keyVarObsCatch'])]
  summary(d$keyCatch)
d$logCatchObs<-d$logCatchObs[d$keyCatch[,'obs.no']]
d$keyCatch[,'obs.no']<-dim(d$keyCatch)[[1]]


$ keySurvey                : int [1:5001, 1:9] 1 2 3 4 5 6 7 8 9 10 ...
..- attr(*, "dimnames")=List of 2
$ keySurvey.overview       : int [1:34, 1:19] 1 2 3 4 5 6 7 8 9 10 ...
..- attr(*, "dimnames")=List of 2
$ logSurveyObs             : num [1:5001] 12.9 14 11 14 13.3 ...
$ sampleTimeWithinSurvey   : num [1:34, 1] 0.5 0.5 0 0.5 0.5 0.5 0 0.05 0.05 0 ...
..- attr(*, "dimnames")=List of 2

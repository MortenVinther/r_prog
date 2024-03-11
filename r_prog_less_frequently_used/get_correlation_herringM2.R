useLog<-0  # 0= no log N and no log F, 1=log N and (no log) F, 2=log N and log F

# extract covariance or correlation matrix

fit<-read.fit()
#str(fit)


sort(unique(fit$names))
getCorCov<-function(vname,fit){
    idx<-which(fit$names %in% vname)
    COR<-fit$cor[idx,idx]
    COV<-fit$cov[idx,idx]
    labels<-paste(paste("Y",fit$year[idx],"A",fit$age[idx],sep=''),sep='.')
    dimnames(COR)<-list(labels,labels)
    dimnames(COV)<-list(labels,labels)
    val<-fit$est[idx]
    return(list(cov=COV,cor=COR,val=val,species=fit$species[idx],age=fit$age[idx],vari=fit$names[idx],year=fit$year))
}

avgM2<-getCorCov(vname=c("her_avgM2"),fit=fit)
str(avgM2)
round(avgM2$cor,3)
avgM2$cov

M2<-getCorCov(vname=c("her_M2"),fit=fit)
str(M2)
round(M2$cor,3)



sort(unique(fit$names))
getCorCov<-function(vname,fit){
  idx<-which(fit$names %in% vname)
  COR<-fit$cor[idx,idx]
  COV<-fit$cov[idx,idx]
  labels<-paste(fit$names[idx],paste("Age",fit$age[idx],sep=''),sep='.')
  dimnames(COR)<-list(labels,labels)
  dimnames(COV)<-list(labels,labels)
  val<-fit$est[idx]
  return(list(cov=COV,cor=COR,val=val,species=fit$species[idx],age=fit$age[idx],vari=fit$names[idx]))
}

a<-Read.SMS.std()

save(a,M2,avgM2,file='herrringM2.Rdata',ascii=TRUE)



keepPartialM2<-FALSE
M2part<- lapply(1:dim(suitIdx)[[1]],function(yqa){
  x1<- suitIdx[yqa,]
  # print(x1)
  y<-x1$y
  q<-x1$q
  yq<-yqIdx[y,q]
  area<-x1$area
  localNq<-rbind(suppressWarnings(do.call(rbind,lapply(1:nSpecies,function(x)logNq[[x]][y,q,]))),rep(0,nAges)) # reformat prey abundance  to allow vectorization
  cat('M2 y:',y,'  q:',q,'\n')
  M2local[yq][]<-0
  do.call(rbind,lapply(1:dim(x1$data[[1]])[[1]],function(yqapa){
    x2<-x1$data[[1]][yqapa,]
    predNo<-x2$predNo
    predAge<-x2$predAge
    predW<-x2$predW
    localNq[otherFoodIdx,]<-logOtherFood[predNo]
    predAbun<-logNq[[predNo]][y,q,predAge]
    predCons<-consum[[predNo]][y,q,predAge]
    
    xx<-x2$data[[1]] %>%
      transmute(preyNo,preyAge,suit=suitability(q,pred=predNo,prey=preyNo, predSize=predW, preySize=preyW, ratio=logRatio, vulIdx=vulneraIdx),
                availFood=exp(localNq[cbind(preyNo,preyAge)]+preyW+suit),
                M2=exp(predAbun+suit)*predCons)
    availFood[[yq]][predNo,predAge]<<-sum(xx$availFood)
    xx$M2<-xx$M2/sum(xx$availFood)
    if (keepPartialM2) transmute(xx,y,q,predNo,predAge,preyNo,preyAge,M2) else transmute(xx,y,q,preyNo,preyAge,M2)
  })
  )
})


M2part<-subset(do.call(rbind,M2part),preyNo<=nSpecies) 

M2<-by(M2part,M2part$preyNo,function(x) {y<-tapply(x$M2,list(x$y,x$q,x$preyAge),sum); y[is.na(y)]<-0; y})
names(M2)<-names(preys)

#lapply(M2,function(x) round(apply(x,c(1,3),sum),3))
for (i in names(M2)) {
  dd<-dim(M2[[i]])
  d1<-1L:dd[1]; d2<-1L:dd[2]; d3<-1L:dd[3]
  MM[[i]][d1,d2,d3]<-natMor1[[i]][d1,d2,d3]+M2[[i]][d1,d2,d3]
}

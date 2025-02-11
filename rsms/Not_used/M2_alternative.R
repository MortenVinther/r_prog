
for (ayq in (1:min(dim(suitIdx)[[1]],20))) {
  y<-suitIdx[ayq,'y'][[1]]
  q<-suitIdx[ayq,'q'][[1]]
  area<-suitIdx[ayq,'area'][[1]]
  localNq<-rbind(suppressWarnings(do.call(rbind,lapply(1:nSpecies,function(x)logNq[[x]][y,q,]))),rep(0,nAges)) # reformat prey abundance
  
  cat(ayq,' y:',y,'  q:',q,'\n')
  for (ayqpa in (1:dim(suitIdx[ayq,'data'][[1]][[1]])[[1]])) {
    x<-suitIdx[ayq,'data'][[1]][[1]][ayqpa,]
    predNo<-x$predNo
    predAge<-x$predAge
    predW<-x$predW
    localNq[otherFoodIdx,]<-logOtherFood[predNo]
    predAbun<-logNq[[predNo]][y,q,predAge]
    predCons<-consum[[predNo]][y,q,predAge]
    #print(x$data[[1]])
    #if ( ayqpa<3) cat('ayqpa: ',ayqpa,' predNo:',predNo,' predAge:',predAge,'predW:',predW,'\n')
    # suitIdx[ayq,'data'][[1]][[1]][ayqpa,'data'][[1]][[1]] <- suitIdx[ayq,'data'][[1]][[1]][ayqpa,'data'][[1]][[1]] %>% rowwise() %>%
    #         mutate(suit=suitability(q,pred=predNo,prey=preyNo, predSize=predW, preySize=preyW, ratio=logRatio, vulIdx=vulneraIdx),
    #            availFood=exp(logNq[[preyNo]][y,q,preyAge]+preyW+suit),
    #            preM2=exp(predAbun+suit)*predCons
    #     )
    xx<-x$data[[1]] %>%
      transmute(preyNo,preyAge,suit=suitability(q,pred=predNo,prey=preyNo, predSize=predW, preySize=preyW, ratio=logRatio, vulIdx=vulneraIdx),
                availFood=exp(localNq[cbind(preyNo,preyAge)]+preyW+suit),
                M2=exp(predAbun+suit)*predCons)
    xx$M2<-xx$M2/sum(xx$availFood)
    
  }
}

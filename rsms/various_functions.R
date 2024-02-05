announce<-function(x) {
  cat("\nobjective:",x$objective,"  convergence:",x$convergence, "  ", x$message, "  iterations:",x$iterations, "  evaluations:",x$evaluations)
}

inputToDF<-function(data) {
  a<-rbind(
    list_rbind(lapply(data$catchMeanWeight,array2DF),names_to='s') %>% mutate(var="weca"),
    list_rbind(lapply(data$catchNumber,array2DF),names_to='s') %>% mutate(var="canum"),
    list_rbind(lapply(data$stockMeanWeight,array2DF),names_to='s') %>% mutate(var="west"),
    list_rbind(lapply(data$propMat,array2DF),names_to='s') %>% mutate(var="propMat"),
    list_rbind(lapply(data$propM,array2DF),names_to='s') %>% mutate(var="propM"),
    list_rbind(lapply(data$propF,array2DF),names_to='s') %>% mutate(var="propF"),
    list_rbind(lapply(data$natMor,array2DF),names_to='s') %>% mutate(var="M"),
    list_rbind(lapply(data$seasFprop,array2DF),names_to='s') %>% mutate(var="seasFprop")
  )
  pivot_wider(a,names_from=var,values_from=Value) %>% mutate_if(is.character,as.integer) %>%
    mutate(species=data$spNames[s], year=Var1-data$off.year, q=Var2,age=Var3-data$off.age) %>% 
    mutate(Var1=NULL,Var2=NULL,Var3=NULL)
}

#a<-inputToDF(data)


outputToDF<-function(obj) {
  rep<-obj$report()
  a<-rbind(
    list_rbind(lapply(rep$logNq,array2DF),names_to='s') %>% mutate(var="N",Value=exp(Value)),
    list_rbind(lapply(rep$Zq,array2DF),names_to='s') %>% mutate(var="Z")
  )  
  pivot_wider(a,names_from=var,values_from=Value)  %>%
    mutate(species=s, year=as.integer(y),q=as.integer(q),age=parse_number(a)-data$off.age) %>% 
    mutate(y=NULL,a=NULL,Var3=NULL,s=NULL)
} 
#  b<-outputToDF(obj)



FToDF<-function(obj,sdrep) {
  if (missing(sdrep)) sdrep <- sdreport(obj)
  x<-as.list(sdrep, what="Est")
  ff<-exp(x$Uf)
  colnames(ff)<-data$years
  FF<-NULL
  for (s in (1:data$nSpecies)) {
    i<-data$nlogFfromTo[s,]
    key<-data$keyLogFsta[s,][data$keyLogFsta[s,]>0]
    faf<-data$info[s,'faf']; la<-data$info[s,'la']
    fff<-ff[i[1]:i[2],,drop=FALSE][key,]
    rownames(fff)<-(faf:la)-data$off.age
    fff<-array2DF(fff) %>% mutate(s=s,year=as.integer(Var2),age=as.integer(Var1),Var1=NULL,Var2=NULL) %>%rename(FF=Value)
    FF<-rbind(FF,fff)
  }
  FF
}
#  FF<-FToDF(obj,sdrep=sdrep)


#b<-left_join(inputToDF(data),outputToDF(obj),by = join_by(species, year, q, age))

  

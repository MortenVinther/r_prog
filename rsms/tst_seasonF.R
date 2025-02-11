
for (s in (1:data$nSpecies)){
  x<-data$seasFprop[[s]]
  y<-ftable(apply(x,c(1,3),sum))
  cat(data$spNames[s],'\n')
  print(y)
}


for (s in (1:data$nSpecies)){
  x<-data$seasFprop[[s]]
  y<-ftable(apply(x,c(1,3),function(x){sum(x)>0}))
  cat(data$spNames[s],'\n')
  print(y)
}

ageQCatchExist<-lapply(data$seasFprop,function(x) apply(x,c(1,3),function(x){sum(x)>0}))


N<-filter(sms$rep$res,quarter==1) %>% mutate(Catch=NULL, wCatch=NULL)
nPred<-filter(sms$rep$resAnnual) %>% mutate(Chat=NULL, Yield=NULL, M2=NULL)  
a<-left_join(N,nPred,by = join_by(year, age, species)) %>% mutate(processZ=log(N/predN))

filter(a,species=='SSA')

a<-filter(sms$rep$res,species=='SSA' &year==1999)
a

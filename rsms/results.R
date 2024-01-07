
rep<-obj$report()

#lapply(rep$logNq,function(x) (exp(x[,1,])))
#if (data$nSeasons==4)lapply(rep$logNq,function(x) (exp(x[,3,])))

#lapply(rep$Zq,function(x) round(exp(x[,1,]),2))
#lapply(rep$Chat,function(x) round(exp(x)))

cbind(rep$nlls,all=rowSums(rep$nlls))

Est<-as.list(rep, "Est", report=TRUE)

sdrep <- sdreport(obj)
# obj$fn()  #  value er den samme som sum af min nnls
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

# collect data to summary plot

convert_var<-function(x) {
  xx<-lapply(x,function(x) {
    as.data.frame(x)  %>% mutate(year=as.numeric(rownames(x))) %>% 
      pivot_longer(!year, names_to = "Age_idx", values_to = "N")
  })
  for (s in names(x)) xx[[s]]$species<-s
  xx<-do.call(rbind,xx)
  xx$Age<-as.integer(matrix(unlist(strsplit(xx$Age_idx,' ')),ncol=2,byrow=TRUE)[,2]) -data$off.age
  xx  
}

if (data$nSeasons==1) N<-lapply(rep$logNq,function(x) (exp(x[,1,])))
if (data$nSeasons>1) N<-lapply(rep$logNq,function(x) (exp(x[,3,])))
Recruit<-convert_var(N) %>% filter(Age==0) %>% rename(Species=species,Year=year)

sdrep <- sdreport(obj)
x<-as.list(sdrep, "Est", report=TRUE)

ssb<-x$ssb
colnames(ssb)<-data$years
rownames(ssb)<-data$spNames
ssb<-array2DF(ssb); colnames(ssb)<-c('Species','Year','SSB')
ssb<- ssb %>% mutate(Year=as.numeric(Year))

x<-as.list(sdrep, "Est", report=FALSE)
ff<-exp(x$Uf)

avg_F<-NULL
for (s in (1:data$nSpecies)) {
  i<-data$nlogFfromTo[s,]
  key<-data$keyLogFsta[s,][data$keyLogFsta[s,]>0]
  fff<-ff[i[1]:i[2],,drop=FALSE][key,]
  data$fbarRange[s,]
  ageF<-data$fbarRange[s,]+data$off.age-data$info[s,'faf']+1
  
  avg_F<-rbind(avg_F,data.frame(Year=data$year,Species=data$spNames[s],Species.n=s,mean.F=apply(fff[ageF[1]:ageF[2],],c(2),mean)))
}

rsms<-left_join(left_join(ssb,Recruit),avg_F) %>% select(Species, Year, SSB, N, mean.F) %>% rename(Rec=N) %>% mutate(source='rsms')


sms<-Read.summary.table(dir=file.path(root,SMSenv),read.init.function=TRUE) %>% select(Species,Year,Rec,SSB,mean.F)
sms$source='sms'


rsp<-unique(rsms$Species)

b<-rbind(rsms,filter(sms,Species %in% rsp )) %>% filter(Year %in% data$years) %>% mutate(Rec=Rec/1000)
b<-pivot_longer(b,cols=c(Rec,SSB,mean.F),names_to='variable') %>% mutate_if(is_character, as.factor)


for (s in (rsp)) {
  bb=filter(b,Species==s)
  
  p<-ggplot(data=bb, aes(x=Year, y=value, group=source)) +
    geom_line(aes(linetype=source,col=source))+
    geom_point(aes(shape=source,col=source))+
    facet_grid(variable ~ ., scales="free_y")+
    ggtitle(unlist(bb[1,'Species']))
  print(p)
  cat('press return to see the next plot:')
  readLines(n=1)
  cat('\n')
}



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



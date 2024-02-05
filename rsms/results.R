t(t(opt$par))
rep<-obj$report()

overw<-as.data.frame(data$keySurvey.overview) %>%  rownames_to_column(var = "fleet") %>% select(f,fleet,s,type)

survResid<-as.data.frame(data$keySurvey) %>% as_tibble()  %>% 
  mutate(predict=rep$predSurveyObs,obs=data$logSurveyObs,Residual=predict-obs,year=y-data$off.year,Age=factor(a-data$off.age),species=data$spNames[s]) %>% 
  left_join(.,overw,by = join_by(f, s)) %>% select(s,species,f,fleet,type,year,q,Age,s,predict,obs,Residual) %>%
  mutate(cols=if_else(Residual<0,'negativ','positiv')) 

plt<-by(survResid,survResid$s,function(x){
  plt<-x %>% ggplot(aes(year, Age,color= cols,fill=cols,size = sqrt(abs(Residual)))) +
    geom_point(shape = 21,alpha=0.75) +
    scale_color_manual(values = c("blue", "red")) +
    labs(x = "", y = "Age",title=x[1,'species']) +
   # facet_grid(rows =vars(fleet), scales="free_y")
  facet_wrap(vars(fleet),ncol=1)+
    theme_bw() + 
    theme(strip.text = element_text(size = 10,margin = margin(0.01,0,0.01,0, "cm")))
  print(plt)
})
plt


N<-lapply(rep$logNq,function(x) (cbind(a1=x[,data$recSeason,1],x[,1,-1]))) # N in recruiting season for the first age and age 1 jan for the rest
# Process N
resid<-lapply(data$spNames,function(x) {
  y<-t(N[[x]])-rep$predN[[x]]; 
  rownames(y)<-(1:dim(y)[[1]])-data$off.age
  y<-y[,-1]
  attr(y,'species')<-x; 
  y
})

lapply(resid,function(x) {
  as.data.frame(x) %>%
    rownames_to_column(var = "Age") %>%
    gather(key, Residual, -Age) %>% as_tibble() %>%
    mutate(Age=factor(as.integer(Age)),year=as.integer(key),cols=if_else(Residual<0,'negativ','positiv')) %>% #mutate(cols=factor(cols)) %>%
    ggplot(aes(year, Age,color= cols,fill=cols,size = sqrt(abs(Residual)))) +
    geom_point(shape = 21,alpha=0.75) +
    scale_color_manual(values = c("blue", "red")) +
    labs(x = "", y = "Age",title=attr(x,'species')) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5))
})
  
#lapply(rep$Zq,function(x) round(exp(x[,1,]),2))
#lapply(rep$Chat,function(x) round(exp(x)))

cbind(rep$nlls,all=rowSums(rep$nlls))

Est<-as.list(rep, "Est", report=TRUE)

sdrep <- sdreport(obj)

sdrep$pdHess 
# obj$fn()  #  value er den samme som sum af min nnls


x<-as.list(sdrep, what="Est")

#summary(sdrep, "random")                      ## Only random effects
#summary(sdrep, "fixed", p.value = TRUE)       ## Only non-random effects
#summary(sdrep, "fixed", p.value = FALSE)    ## Only non-random effects
#summary(sdrep, "report")                      ## Only report
# 
# ssb<-summary(sdrep, "report") 
# plot(ssb[,"Estimate"],type='l')
# lines(ssb[,"Estimate"]+2*ssb[,"Std. Error"],type='l',col='red')
# lines(ssb[,"Estimate"]-2*ssb[,"Std. Error"],type='l',col='red')


parVal<-as.list(sdrep, what="Est")
parSd<-as.list(sdrep, what="Std")

make_tab<-function(d,key,printIt=TRUE,roundIt=2,addSd=TRUE,expIt=TRUE) {
  cat(d,'\n')
  dd<-parVal[[d]]
  if (expIt) dd<-exp(dd)
  tab<-key
  for (i in 1:length(dd)) tab[tab[,]==i]<- dd[i]
  if (addSd) {
    dd<-parSd[[d]]
    if (expIt) dd<-exp(dd)
    tabSd<-key
    for (i in 1:length(dd)) tabSd[tabSd[,]==i]<- dd[i]
    tabSd
  }
  
  if (printIt) {
    ptab<-tab
    ptab[ptab== -1] <-NA
    print(round(ptab,roundIt),na.print='.')
    if (addSd) {
      ptab<-tabSd
      ptab[ptab== -1] <-NA
      print(round(ptab,roundIt),na.print='.')
    }
    
  } 
  invisible(tab)
}

addSd<-TRUE;expIt<-T
make_tab(d="logSdLogN",key=data$keyVarLogN,roundIt=2,addSd=addSd,expIt=expIt)


#round(exp(x$Uf),3)
#round(exp(x$Un),0)

x$rec_loga; exp(x$rec_loga)
x$rec_logb; exp(x$rec_logb)

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

addSd<-TRUE
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

#sdrep <- sdreport(obj)
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

SMSenv<-"rsms_input"
sms<-Read.summary.table(dir=file.path(root,SMSenv),read.init.function=TRUE) %>% select(Species,Year,Rec,SSB,mean.F)
sms$source='sms'

ices<-Read.summary.table(dir=file.path(root,SMSenv), infile="summary_table_raw_ICES.out",read.init.function=TRUE) %>% select(Species,Year,Rec,SSB,mean.F)
ices$source='ICES'

sms<-rbind(sms,ices)

rsp<-unique(rsms$Species)

b<-rbind(rsms,filter(sms,Species %in% rsp )) %>% filter(Year %in% data$years) %>% mutate(Rec=Rec/1000)
b<-pivot_longer(b,cols=c(Rec,SSB,mean.F),names_to='variable') %>% mutate_if(is_character, as.factor)


for (s in (rsp)) {
  bb=filter(b,Species==s)
  
  p<-ggplot(data=bb, aes(x=Year, y=value, group=source)) +
    geom_line(aes(linetype=source,col=source))+
    geom_point(aes(shape=source,col=source))+
    facet_grid(variable ~ ., scales="free_y")+
    ylim(0,NA)+
    ggtitle(unlist(bb[1,'Species']))
  print(p)
  cat('press return to see the next plot:')
  readLines(n=1)
  cat('\n')
}


if (FALSE) for (s in (rsp)) {
  bb=filter(b,Species==s & source=='rsms')
  
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

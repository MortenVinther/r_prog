

#sdrep <- sdreport(obj)
x<-as.list(sdrep, "Est", report=TRUE)
#xstd<-as.list(sdrep, "Std", report=TRUE)

ssb<-x$ssb
colnames(ssb)<-data$years
rownames(ssb)<-data$spNames
ssb<-array2DF(ssb); colnames(ssb)<-c('Species','Year','SSB')
ssb<- ssb %>% mutate(Year=as.numeric(Year))

st<-obj$report()
N<-lapply(st$logNq,function(x) (x[,data$recSeason,1])) # N in recruiting season for the first age 
N<-exp(do.call(rbind,N))
N<-array2DF(N); colnames(N)<-c('Species','Year','Recruit')
N<- N %>% mutate(Year=as.numeric(Year))

ssbRec<-left_join(ssb,N)
stockRecruitmentModelCode<-data$info[,"SSB/R"]


opt$par

rec_loga<-opt$par[grep("rec_loga",names(opt$par))]
if ("rec_loga" %in% names(my.map)) rec_loga<-  rec_loga[my.map$rec_loga]
rec_loga

rec_logb<-opt$par[grep("rec_logb",names(opt$par))]
if ("rec_logb" %in% names(my.map)) rec_logb<-  rec_b[my.map$rec_logb]
rec_logb

recPars<-data.frame(species=names(stockRecruitmentModelCode),species.n=1:length(stockRecruitmentModelCode),model=stockRecruitmentModelCode,rec_loga,rec_logb)




SSB_R<-function(model,a,b,s,species,ssb) {
  
  if(model==0){    ## straight RW
    rec = 0; #logN[[s]][a, y-1]
  } else {
    if (model==1){ ## Ricker
      rec<-a+log(ssb)-exp(b)*ssb
    } else {
      if(model==2){  ## B&H
        rec<-a+log(ssb)-log(1+exp(b)*ssb)
      } else {
        if(model==3){  ## GM
          rec<-a
        } else {
          stop(paste0("SR model code ",model," not recognized"))
        }
      }
    }
  }
  return(rec)
}

SSBmodRec<-lapply(1:length(stockRecruitmentModelCode),function(s) {
  SSB<-filter(ssbRec,Species==recPars[s,'species'])$SSB
  SSB<-seq(10,max(SSB),max(SSB)/100)
  data.frame(Species=recPars[s,'species'],model=recPars[s,'model'],SSB=SSB,Recruit=exp(SSB_R(model=recPars[s,'model'],a=recPars[s,'rec_loga'],b=recPars[s,'rec_logb'],s=recPars[s,'species.n'],
        species=recPars[s,'species'],ssb=SSB)))
})

SSBmodRec<-do.call(rbind,SSBmodRec) %>% filter(model>0) %>% mutate(Recruit=Recruit/1000,SSB=SSB/1000)

ggplot(data = ssbRec, mapping = aes(x =SSB/1000, y = Recruit/1000)) +
  geom_point() +
  geom_line(mapping = aes(x =SSB,y =Recruit), data=SSBmodRec,  colour = 'red',size=1)+
  ylim(0,NA)+xlim(0,NA)+ylab('Recruits (million)')+xlab('SSB (1000 t)')+

facet_wrap(vars(Species), scales="free")



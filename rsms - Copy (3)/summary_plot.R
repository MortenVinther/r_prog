summary_plot<-function(obj,data,sdrep,out,incl_ICES_plot=TRUE,printPlot=FALSE,pLabel='') {
  # collect data to summary plot
  
  rep<-obj$report()
  
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
  if (data$nSeasons>1) N<-lapply(rep$logNq,function(x) (exp(x[,data$recSeason,])))
  Recruit<-convert_var(N) %>% filter(Age==0) %>% rename(Species=species,Year=year)
  
  ssb<-out[['ssb']]
  colnames(ssb)<-data$years
  rownames(ssb)<-data$spNames
  ssb<-array2DF(ssb); colnames(ssb)<-c('Species','Year','SSB')
  ssb<- ssb %>% mutate(Year=as.numeric(Year))
  
  #x<-as.list(sdrep, "Est", report=FALSE)
  #ff<-exp(x$Uf)
  ff<- FToDF(obj,sdrep,data) %>% mutate(species=data$spNames[s])
  
  avg_F<-NULL
  FAges<-data.frame(data$fbarRange,species=data$spNames)
  
  avg_F<-left_join(ff,FAges,by = join_by(species)) %>% filter(age>=first.age & age<=last.age) %>%
    group_by(s,species,year) %>% mutate(mean.F=mean(FF),Species.n=s) %>% ungroup() %>%
    transmute(Species=species,Year=year,mean.F)
  
  rsms<-left_join(left_join(ssb,Recruit),avg_F) %>% select(Species, Year, SSB, N, mean.F) %>% rename(Rec=N) %>% mutate(source='rsms')
  
  SMSenv<-my.stock.dir
  sms<-Read.summary.table(dir=file.path(root,SMSenv),read.init.function=TRUE) %>% select(Species,Year,Rec,SSB,mean.F)
  sms$source='sms'
  
  if (incl_ICES_plot) {
    ices<-Read.summary.table(dir=file.path(root,SMSenv), infile="summary_table_raw_ICES.out",read.init.function=TRUE) %>% select(Species,Year,Rec,SSB,mean.F)
    ices$source='ICES'
    sms<-rbind(sms,ices)
  }
  #xtabs(~Species+source,data=sms)
  rsp<-unique(rsms$Species)
  
  b<-rbind(rsms,filter(sms,Species %in% rsp )) %>% filter(Year %in% data$years) %>% mutate(Rec=Rec/1000)
  b<-pivot_longer(b,cols=c(Rec,SSB,mean.F),names_to='variable') %>% mutate_if(is_character, as.factor)
  
  
  plts<-lapply(rsp, function(s)  {
    bb=filter(b,Species==s)
    
    p<-ggplot(data=bb, aes(x=Year, y=value, group=source)) +
      geom_line(aes(linetype=source,col=source))+
      geom_point(aes(shape=source,col=source))+
      facet_grid(variable ~ ., scales="free_y")+
      ylim(0,NA)+
      ggtitle(paste(unlist(bb[1,'Species']),pLabel))
    if (printPlot) print(p)
    #cat('press return to see the next plot:')
    #readLines(n=1)
    #cat('\n')
    return(p)
   })
  
  
  return(plts)
}

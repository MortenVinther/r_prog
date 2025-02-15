
my_theme <- function() {
  theme_bw()+
    theme( #panel.grid.major = element_line(linetype = "blank"),
      panel.grid.minor = element_line(linetype = "blank"),
      #axis.text.x = element_text(angle = tAngle, vjust = 0.5),
      legend.title = element_blank(),
      legend.position = "none",
    )  
}

plotSSB<-function(x,sp,mult=1E-3,tit='SSB') {
  x<-filter(x,s==sp) %>% mutate(SSB=SSB*mult)
  a<-ggplot(data=x, aes(x=year, y=SSB,col=run,shape=run,group=run)) +
    geom_line()+
    geom_point()+
    ylim(0,NA)+
    labs(x='Year',y='(1000 tonnes)',title=tit)+
    my_theme()
  a
}


plotSSBRibbon<-function(x,sp,mult=1E-3,tit='SSB') {
  x<-filter(x,s==sp ) %>% mutate(mid=mid*mult,low=low*mult,high=high*mult)
  a<-ggplot(data=x, aes(x=year,y=mid,ymin=low, ymax=high,fill=run,col=run,shape=run,group=run)) +
    geom_line()+
    geom_point()+
    geom_ribbon(alpha=0.25) +
    ylim(0,NA)+
    labs(x='Year',y='(1000 tonnes)',title=tit)+
    my_theme()
  a
}


plotLabels<-function(tit='Hej',labels){
  nlabels<-length(labels)
  ll<-max(nchar(labels))
  ggplot(data.frame(x=1,y=1,run=factor(labels,labels)),
         aes(x,y,col=run,shape=run))+
    geom_line(linewidth = 1.5)+
    geom_point(aes(shape=run),size=4)+
    xlim(2,3)+ylim(2,3)+
    #labs(color= tit)+
    #guides( size = FALSE,shape=FALSE)+
    labs(title=tit)+
    theme_void()+
    theme(plot.title = element_text(color="black", size=25, face="bold",margin=margin(t=20, r=0,b =0, l=15, unit="pt")),
          legend.title = element_blank(),
          legend.text=element_text(size=17-ll/5),
          legend.position = "inside",
          legend.position.inside=c(0.15,0.4+0.04*nlabels),
          legend.justification = c("left"),
          legend.box.margin = margin(15, 15, 15, 15),
          legend.key.size = unit(3, "lines"),
          legend.key.height = unit(1.3, "lines"),                
          legend.key.width = unit(2, "lines"), 
          legend.box.background = element_rect(color="black", linewidth=1))
}  

plotRecruits<-function(x,sp,mult=1E-6,tit='Recruits') {
  x<-filter(x,s==sp) %>% mutate(recruit=recruit*mult)
  
  a<-ggplot(data=x, aes(x=year, y=recruit,col=run,shape=run,group=run)) +
    geom_line()+
    geom_point()+
    ylim(0,NA)+
    labs(x='Year',y='(billions)',title=tit)+
    my_theme() 
  a
}

plotM2<-function(x,sp) {
  xx<-filter(x,s==sp) 
  showAges<-sort(unique(x$age))
  p<-lapply(showAges,function(a){
    tit<-paste('Age:',a)
    xxx<-filter(xx,age==a)%>% mutate(age=paste('age:',age))
    #print(xtabs(M2~run+year,data=xxx))
    pp<-ggplot(data=xxx, aes(x=year, y=M2,col=run,shape=run,group=run)) +
    geom_line()+
    geom_point()+
    ylim(0,NA)+
    labs(x='Year',y=' ',title=tit)+
    #  facet_wrap(~age)+
    my_theme()
   })
}


plotFage<-function(x,sp) {
  xx<-filter(x,s==sp) 
  showAges<-sort(unique(x$age))
  p<-lapply(showAges,function(a){
    tit<-paste('Age:',a)
    xxx<-filter(xx,age==a)%>% mutate(age=paste('age:',age))
    #print(xtabs(M2~run+year,data=xxx))
    pp<-ggplot(data=xxx, aes(x=year, y=FisQ,col=run,shape=run,group=run)) +
      geom_line()+
      geom_point()+
      ylim(0,NA)+
      labs(x='Year',y=' ',title=tit)+
      #  facet_wrap(~age)+
      my_theme()
  })
}

plotNage<-function(x,sp) {
  xx<-filter(x,s==sp) 
  showAges<-sort(unique(x$age))
  p<-lapply(showAges,function(a){
    tit<-paste('Age:',a)
    xxx<-filter(xx,age==a)%>% mutate(age=paste('age:',age))
    #print(xtabs(M2~run+year,data=xxx))
    pp<-ggplot(data=xxx, aes(x=year, y=N,col=run,shape=run,group=run)) +
      geom_line()+
      geom_point()+
      ylim(0,NA)+
      labs(x='Year',y=' ',title=tit)+
      #  facet_wrap(~age)+
      my_theme()
  })
}


plotRecruitsRibbon<-function(x,sp,mult=1E-6,tit='Recruits') {
  x<-filter(x,s==sp ) %>% mutate(mid=mid*mult,low=low*mult,high=high*mult)
  a<-ggplot(data=x, aes(x=year,y=mid,ymin=low, ymax=high,fill=run,col=run,shape=run,group=run)) +
    geom_line()+
    geom_point()+
    geom_ribbon(alpha=0.25) +
    ylim(0,NA)+
    labs(x='Year',y='(billions)',title=tit)+
    my_theme()
  a
}

plotF<-function(x,sp,tit='F') {
  x<-filter(x,s==sp) 
  a<-ggplot(data=x, aes(x=year, y=Fbar,col=run,shape=run,group=run)) +
    geom_line()+
    geom_point()+
    ylim(0,NA)+
    labs(x='Year',y=' ',title=tit)+
    my_theme() 
  a
}


plotFRibbon<-function(x,sp,tit='F') {
  x<-filter(x,s==sp ) 
  a<-ggplot(data=x, aes(x=year,y=mid,ymin=low, ymax=high,fill=run,col=run,shape=run,group=run)) +
    geom_line()+
    geom_point()+
    geom_ribbon(alpha=0.25) +
    ylim(0,NA)+
    labs(x='Year',y=' ',title=tit)+
    my_theme()
  a
}


# 1. Summary plots (Recruit, SSB and F), option "compSummary", 
# 2. Summary pots with uncertainties, option "compSummaryConf" 
# 3. M2 by age, option "compM2" input rep$res

if (FALSE) {
  #parameters to function call
  Type<-"compSummaryConf"
  #Type<-"compSummary"
  showSpecies<-c(1,2,3)
  showSpecies<-1:10
  inpRdata<-list('Single','Multi')
  labels=c('Single sp','Multi sp')
  
  #Type<-"compSummary"
  #inpRdata<-list('Single','ICES_single_sp')
  #labels=c('Single sp','ICES')
  
  outFormat<-c('screen','pdf','png') [1]
  longSpNames<-TRUE
  showAges=0:5
}

####
plotCompareRunSummary<-function(Type=c("compSummaryConf","compSummary","compM2","compF","compN"),showSpecies=1:100,
                        inpRdata=list('Single','Multi'),
                        labels=c('Single sp','Multi sp'),
                        outFormat=c('screen','pdf','png')[1],
                        showAges=0:4,
                        multN=0.001,
                        ncol=2,
                        fileLabel='pl',
                        longSpNames=TRUE){
   
  if (Type=="compSummary") {
    x<-do.call(rbind,lapply(1:length(labels),function(i){
      load(file=file.path(data.path,paste0(inpRdata[i],".Rdata")),verbose=F)
      sms$rep$resSummary %>% mutate(year=as.numeric(year),run=labels[i]) %>% filter(s %in% showSpecies)
    })) 
  } else if (Type=="compSummaryConf"){
    x<-do.call(rbind,lapply(1:length(labels),function(i){
      load(file=file.path(data.path,paste0(inpRdata[i],".Rdata")))
      out<-data.frame(value=sms$sdrep$value,sd=sms$sdrep$sd,variable=names(sms$sdrep$value)) 
      nvar<-3
      cbind(expand.grid(s=1:sms$data$nSpecies,year=sms$data$years,dummy=1:nvar),out) %>% filter(s %in% showSpecies) %>%
         mutate(run=labels[i],species=sms$data$spNames[s],mid=exp(value),low=exp(value-2*sd),high=exp(value+2*sd))
  
    }))
  } else if (Type=='compM2'){
    x<-do.call(rbind,lapply(1:length(labels),function(i){
      load(file=file.path(data.path,paste0(inpRdata[i],".Rdata")))
      out<-sms$rep$res %>% transmute(year,age,s,species,run=labels[i],M2) %>% 
        filter(s %in% showSpecies & age %in% showAges &M2>0) %>% group_by(s,species,year,age,run) %>%
        summarize(M2=sum(M2),.groups='drop') %>% group_by(s,age) %>% mutate(maxM2=max(M2)) %>% filter(maxM2>1E-5) %>%
        mutate(maxM2=NULL)
    }))
    
  } else if (Type=='compF'){
  x<-do.call(rbind,lapply(1:length(labels),function(i){
    load(file=file.path(data.path,paste0(inpRdata[i],".Rdata")))
    if (c('F') %in% colnames(sms$rep$res)) sms$rep$res$FisQ<-sms$rep$res$F
    if (c('FF') %in% colnames(sms$rep$res)) sms$rep$res$FisQ<-sms$rep$res$FF
    out<-sms$rep$res %>% transmute(year,age,s,species,run=labels[i],FisQ) %>% 
      filter(s %in% showSpecies & age %in% showAges & FisQ>0) %>% group_by(s,species,year,age,run) %>%
      summarize(FisQ=sum(FisQ),.groups='drop') %>% group_by(s,age) %>% mutate(maxF=max(FisQ)) %>% filter(maxF>1E-5) %>%
      mutate(maxF=NULL)
  }))
  } else if (Type=='compN'){
    x<-do.call(rbind,lapply(1:length(labels),function(i){
      load(file=file.path(data.path,paste0(inpRdata[i],".Rdata")))
      out<-sms$rep$res %>% transmute(year,quarter,age,s,species,run=labels[i],N=N*multN) %>% 
        filter(s %in% showSpecies & age %in% showAges & N>0 & (quarter==data$fqa[age+data$off.age])) %>% 
        mutate(quarter=NULL)
    }))
  }
  #3xtabs(M2~paste(s,species)+run,data=out)
  #xtabs(N~paste(s,species)+age,data=sms$rep$res)
  x<- x %>% mutate(run=factor(run,labels))
  
  load(file=file.path(data.path,paste0(inpRdata[1],".Rdata")))
  allSpNamesLong<-sms$data$allSpNamesLong
  names(allSpNamesLong)<-sms$data$allSpNames
  showSpecies<-showSpecies[showSpecies %in% sort(unique(x$s))]
  
  cleanup()
  for (mysp in showSpecies) {
     tit<-allSpNamesLong[mysp]
     titFile<-tit
     if (Type %in% c("compM2")) tit<-paste('M2\n',tit)
     if (Type %in% c("compF"))  tit<-paste('F\n',tit)
     if (Type %in% c("compN"))  tit<-paste('N\n',tit)
     myLab<-plotLabels(tit,labels)
   
     if (Type=="compSummary") {
       pSSB<-plotSSB(x,sp=mysp)
       pF<-plotF(x,sp=mysp)
       pRecruit<-plotRecruits(x,sp=mysp)
       
     } else if (Type=="compSummaryConf"){
        pSSB<-plotSSBRibbon( x=filter(x,variable=='SSB'),sp=mysp,mult=1E-3)
        pRecruit<-plotRecruitsRibbon(x=filter(x,variable=='recruit'),sp=mysp,mult=1E-6)
        pF<-plotFRibbon(x=filter(x,variable=='FBAR'),sp=mysp)
    
    }else if (Type=="compM2"){
       pAge<-plotM2(x,sp=mysp) 
  
    } else if (Type=="compF"){
       pAge<-plotFage(x,sp=mysp) 
    
    } else if (Type=="compN"){
      pAge<-plotNage(x,sp=mysp) 
  } 
   if (Type %in% c("compSummary","compSummaryConf")) pl<-suppressWarnings(print(ggarrange(myLab,pSSB,pRecruit,pF, ncol = 2, nrow = 2)))
   if (Type %in% c("compM2","compF","compN")) {
       nc<-ncol
       nr<-(length(showAges)+1) %/% nc + (length(showAges)+1) %% nc
       pl<-suppressWarnings(print(ggarrange(myLab,plotlist=lapply(pAge,print), ncol = nc, nrow = nr)))
   }
    if (outFormat=='screen'){
     print(pl)
    } else {
       ggexport(pl, filename = paste0(fileLabel,'_',Type,'_',titFile,'.',outFormat),width = 900,height = 600, pointsize = 16)
       cleanup()
    }
  }
}


plotSeasoData<-function(x,tit='F',ncols=2) {
  a<-ggplot(data=x, aes(x=year, y=value,col=quarter,shape=quarter,group=quarter)) +
    geom_line()+
    geom_point()+
    ylim(0,NA)+
    labs(x='Year',y=' ',title=tit)+
    my_theme() +  theme(legend.title = element_blank(),    legend.position = "right")+
    facet_wrap(vars(age),ncol=ncols, scales="free_y")
    
  a
}



plotSeasonalData<-function(inp,Type=c("N","F","C","M","M1",'"M2','Z',"WEST","WECA","propMat","seasFprop","FiProp")[2],
                                CombineSeason=FALSE,
                                showSpecies=1:100,
                                outFormat=c('screen','pdf','png')[1],
                                showAges=0:4,
                                multN=0.001,
                                ncols=2,
                                fileLabel='pl',
                                cummulate=FALSE,  
                                longSpNames=TRUE){
 
  a<-c("N","F",   "C",    "M","M1",'"M2','Z',"WEST","WECA","propMat","seasFprop","FiProp")
  b<-c("N","FisQ","canum","M","M1",'"M2','Z','west',"weca","propMat","seasFprop","FiProp")
  multiplier<-rep(1,length(a)); multiplier[1]<-multN
  labels<-c("Stock numbers","Fishing mortality","Catch numbers at age","M","M1",'"M2','Z, total Mortality','Weight at age in the stock',"Weight at age in the catch","propMat","seasFprop","proportion of F")
  names(b)<-a
  names(labels)<-a
  names(multiplier)<-a
  
  load(file=file.path(data.path,paste0(inp,".Rdata")),verbose=T)
  recSeason<-sms[['data']]$recSeason
  spNames<-sms[['data']]$allSpNames
  spNamesLong<-sms[['data']]$allSpNamesLong
  cat(spNamesLong,'\n')
  x<-sms[['rep']]$res %>% filter(s %in% showSpecies & age %in% showAges) 
  x$value<-unlist(x[,b[Type]])
  x<-x  %>% transmute(species,s,year,quarter=factor(quarter),age=paste('Age',age),value=value*multiplier[Type])
  
  
  if (CombineSeason) {
    if (Type=='N') filter(x,(a>1 & Quarter==1) | (a==1 & quarter==recSeason)) 
    ## more to come ....  
  }
  if (cummulate) {
    x<-x %>% arrange(s,year,age,quarter) %>% group_by(s,year,age) %>% mutate(value=cumsum(value)) 
  }
  cleanup()
  p<-by(x,x$s,function(x){
    pl<-plotSeasoData(x,tit=paste0(spNamesLong[unlist(x[1,'s'])],', ',labels[Type]),ncol=ncols) 
    if (outFormat=='screen'){
      print(pl)
    } else {
      ggexport(pl, filename = paste0(fileLabel,'_',Type,'_',titFile,'.',outFormat),width = 900,height = 600, pointsize = 16)
      cleanup()
    }
  })
 
}  
 



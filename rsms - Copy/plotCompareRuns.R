
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
  M2ages<-sort(unique(x$age))
  p<-lapply(M2ages,function(a){
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
  M2ages=0:5
}

####
plotCompareRunSummary<-function(Type=c("compSummaryConf","compSummary","compM2"),showSpecies=1:100,
                        inpRdata=list('Single','Multi'),
                        labels=c('Single sp','Multi sp'),
                        outFormat=c('screen','pdf','png')[1],
                        M2ages=C(0:4),ncol=2,
                        fileLabel='pl',
                        longSpNames=TRUE){
   
  if (Type=="compSummary") {
    x<-do.call(rbind,lapply(1:length(labels),function(i){
      load(file=file.path(data.path,paste0(inpRdata[i],".Rdata")))
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
        filter(s %in% showSpecies & age %in% M2ages &M2>0) %>% group_by(s,species,year,age,run) %>%
        summarize(M2=sum(M2),.groups='drop') %>% group_by(s,age) %>% mutate(maxM2=max(M2)) %>% filter(maxM2>1E-5) %>%
        mutate(maxM2=NULL)
    }))
  }
  #3xtabs(M2~paste(s,species)+run,data=out)
  #xtabs(N~paste(s,species)+age,data=sms$rep$res)
  
  x<- x %>% mutate(run=factor(run,labels))
  
  #load(file=file.path(data.path,paste0(inpRdata[1],".Rdata")))
  allSpNamesLong<-sms$data$allSpNamesLong
  names(allSpNamesLong)<-sms$data$allSpNames
  showSpecies<-showSpecies[showSpecies %in% sort(unique(x$s))]
  
  cleanup()
  for (mysp in showSpecies) {
     tit<-allSpNamesLong[mysp]
     myLab<-plotLabels(tit=if_else(Type %in% c("compM2"),paste('M2\n',tit),tit),labels)
   
     if (Type=="compSummary") {
       pSSB<-plotSSB(x,sp=mysp)
       pF<-plotF(x,sp=mysp)
       pRecruit<-plotRecruits(x,sp=mysp)
       
     } else if (Type=="compSummaryConf"){
        pSSB<-plotSSBRibbon( x=filter(x,variable=='SSB'),sp=mysp,mult=1E-3)
        pRecruit<-plotRecruitsRibbon(x=filter(x,variable=='recruit'),sp=mysp,mult=1E-6)
        pF<-plotFRibbon(x=filter(x,variable=='FBAR'),sp=mysp)
    
    }else if (Type=="compM2"){
       pM2<-plotM2(x,sp=mysp) 
  
    }
   if (Type %in% c("compSummary","compSummaryConf")) pl<-suppressWarnings(print(ggarrange(myLab,pSSB,pRecruit,pF, ncol = 2, nrow = 2)))
   if (Type %in% c("compM2")) {
       nc<-ncol
       nr<-(length(M2ages)+1) %/% nc + (length(M2ages)+1) %% nc
       pl<-suppressWarnings(print(ggarrange(myLab,plotlist=lapply(pM2,print), ncol = nc, nrow = nr)))
   }
    if (outFormat=='screen'){
     print(pl)
    } else {
       ggexport(pl, filename = paste0(fileLabel,'_',Type,'_',tit,'.',outFormat),width = 900,height = 600, pointsize = 16)
       cleanup()
    }
  }
}



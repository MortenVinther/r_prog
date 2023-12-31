do_parameters<-function(tables,figures) {
  # read SMS.std and present results
  
  #remove output variables
  #remo<-c('M2_sd','hist_SSB','avg_F')
  remo<-c('hist_SSB','avg_F')
  
  a<-Read.SMS.std(excludeNotUsed=TRUE,remove=remo)
  
  
  # sort parameter uncertainties
  a$absCV.round<-abs(a$CV.round)
  b<-a[order(abs(a$CV.round)),]
  b<-subset(b,select=c(-area,-predator,-prey, -fleet,-index, -used))
  b<-b[order(abs(b$CV.round),b$species),]
  b$CV.Pct<-round(b$CV.round)
  
  
  rownames(b)<-NULL
  b<-subset(b,select=c( name,par,value,std, CV.Pct, species, year, quarter, age))
  b$species.name<-' 0 predation'
  b[b$species>0,"species.name"]<-paste(formatC(b[b$species>0,"species"],width=2),sp.names[b[b$species>0,"species"]])

  if (file.exists( 'par_names.csv')) p<-read.csv('par_names.csv') else p<-read.csv('../HTML_source/par_names.csv')
  
  b<-merge(b,p,by='par')
  b<-droplevels(b)
  a<-xtabs(~par_text+group,data=b)
  a<-rbind(a,all=colSums(a))
  a<-cbind(a,all=rowSums(a))

  
  
  xtab(a, caption=paste("Table ", "Parameter overview. Number of estimated parameters by group of data"), cornername='Parameter',
       file=file.path(tables,'_Parameter_overview.html'), dec=rep(0,dim(a)[[2]]), width='"100%"')
  
  
  
  a<-xtabs(~par_text+species.name,data=b)
  a<-rbind(a,all=colSums(a))
  a<-cbind(a,all=rowSums(a))
  colnames(a)<-trimws(gsub('[[:digit:]]+', '', colnames(a)))

    xtab(a, caption=paste("Table ", "Parameter overview. Number of estimated parameters by species"), cornername='Parameter',
       file=file.path(tables,'_Parameter_species_overview.html'), dec=rep(0,dim(a)[[2]]), width='"100%"')
  
  #plot of year effect in F
  a<-subset(b,par=="F_y_ini")
  a<-subset(a,std<10000)
  library(ggplot2)
  
  p<- ggplot(a, aes(x=year, y=value)) + 
    geom_bar(stat="identity", color="red", position=position_dodge()) +
    #geom_line() +  geom_point()+
    labs(y = "F year effect +-2*sd")+
    facet_wrap(~ species.name, ncol=2,scales='free_y')+
    geom_errorbar(aes(ymin=value-std, ymax=value+std), width=.2,
                  position=position_dodge(.9))
  ggsave(file.path(figures,"uncertanty_F_year_effect.png"),width = 17, height = 22, units = "cm")
  
  
  #plot of recruits
  a<-subset(b,par=="log_rec")

  p<- ggplot(a, aes(x=year, y=value)) + 
    #geom_bar(stat="identity", color="red", position=position_dodge()) +
    geom_line() +  geom_point()+
    labs(y = "log(recruitment) +- 2*sd")+
    facet_wrap(~ species.name, ncol=2,scales='free_y')+
    geom_errorbar(aes(ymin=value-2*std, ymax=value+2*std), width=.2,
                  position=position_dodge(.9))
  ggsave(file.path(figures,"uncertanty_recruit.png"),width = 17, height = 22, units = "cm")
  
   
  #plot of N first year
  a<-subset(b,par=="log_rec_older" | (par=="log_rec" & year==SMS.control@first.year))
  a[a$age==-1,'age']<-SMS.control@first.age
  #a$agec<-formatC(a$age,flag="0",digits=0,format="d",width=2)
  p<- ggplot(a, aes(x=age, y=value)) + 
    #geom_bar(stat="identity", color="red", position=position_dodge()) +
    geom_line() +  geom_point()+
    labs(y = "log(Stock numbers) +- 2*sd",x='Age')+
    facet_wrap(~ species.name, ncol=2,scales='free_y')+
  
  scale_x_continuous(limits = c(SMS.control@first.age, SMS.control@max.age.all),breaks = scales::breaks_width(2))+
    geom_errorbar(aes(ymin=value-2*std, ymax=value+2*std), width=.2,
                  position=position_dodge(.9))
  ggsave(file.path(figures,"uncertanty_N_first_year.png"),width = 17, height = 22, units = "cm")
  
  
  
  ############## predation parameters
  a<-Read.SMS.std(excludeNotUsed=TRUE,remove=remo)
  sort(unique(a$name))
  names<-a$name
  
  withIndex<-c(grep("Stom_var",names),grep("var_size_ratio_ini",names),grep("init_pref_size_ratio",names))
  b<-rbind(subset(a,name %in% c("vulnera","init_stl_other_suit_slope","init_season_overlap")),a[withIndex,])
  b<-b[order(b$index),]
  b$name<-as.character(b$name)
  b$CV.round<-abs(b$CV.round)
  names<-b$name
   #sort(unique(names))
   
   
  found<-grep("Stom_var",names)
  b[found,]
  if (dim(b[found,])[[1]]>0) {
    b[found,'predator']<-readr::parse_number(b[found,'name'])
    b[found,'prey']<- -1
    b[found,'name']<-'Stom_var'
  }

  b[b$predator>0,'Predator']<-paste(formatC(b[b$predator>0,'predator'],width=2),sp.names[b[b$predator>0,'predator']],sep='_')
  #b[b$predator>0,'Predator']<-factor(sp.names[b[b$predator>0,'predator']],levels=sp.names)
  
  b[b$prey>0,'Prey']<-paste(formatC(b[b$prey>0,'prey'],width=2),sp.names[b[b$prey>0,'prey']],sep='_')
  
  bb<-subset(b,name=='vulnera')
  a<-round(tapply(bb$CV.round,list(bb$Predator,bb$Prey),sum))
 # print(a, na.print = "")
  
  colnames(a)<-trimws(gsub('_','',gsub('[[:digit:]]+', '', colnames(a))))
  rownames(a)<-trimws(gsub('_','',gsub('[[:digit:]]+', '', rownames(a))))
  
 if (nrow(a)>1) xtab(a, caption=paste("Table ", "Parameter overview. CV of predator - prey vulnerability parameter"), cornername='Predator',
       file=file.path(tables,'_Parameter_overview_vulnerab_CV.html'), dec=rep(0,dim(a)[[2]]), width='"100%"')
  
  
  other<-unique(subset(bb,select=c(Predator,predator)))
  other$value<-1
  other$std<-0
  other$Prey<-'0_Other food'
  other$prey<-0
  bb<-subset(bb,select=c(Predator,predator,Prey,prey,value,std)) 
  bb<-rbind(bb,other)
  
  a<-round(tapply(bb$value,list(bb$Predator,bb$Prey),sum),2)
  #print(a, na.print = "")
  
  colnames(a)<-trimws(gsub('_','',gsub('[[:digit:]]+', '', colnames(a))))
  rownames(a)<-trimws(gsub('_','',gsub('[[:digit:]]+', '', rownames(a))))
  
  if (nrow(a) >1) xtab(a, caption=paste("Table ", "Parameter overview. Predator - prey vulnerability parameter"), cornername='Predator',
       file=file.path(tables,'_Parameter_overview_vulnerab.html'), dec=rep(2,dim(a)[[2]]), width='"100%"')
  
  summary(bb)
  head(bb)
  #sort(unique(bb$predator));  sort(unique(bb$Predator))
  #sort(unique(bb$prey)); sort(unique(bb$Prey))
  
  bb$Predator<-factor(sp.names[bb$predator],levels=sp.names) %>% droplevels()
  sp.names2<-c("Other food",sp.names)
  #sp.names2
  
  bb$Prey<-factor(sp.names2[bb$prey+1],levels=sp.names2) %>% droplevels()
  
    p<- ggplot(bb, aes(x=Prey, y=value)) + 
    geom_bar(stat="identity", color="red", position=position_dodge()) +
   # geom_line() +  geom_point()+
    labs(y = "Vulnerability +- 1*sd")+
    facet_wrap(~ Predator, ncol=2,scales='free_y')+
    theme(axis.text.x = element_text(angle=90))+
    geom_errorbar(aes(ymin=value-1*std, ymax=value+1*std), width=.2,colour = "red",
                  position=position_dodge(.9))
  ggsave(file.path(figures, "vulnerability.png"),width = 17, height = 22, units = "cm")

  
  
  for (pp in sort(unique(bb$Predator)) ) {
    p<- ggplot(subset(bb,Predator==pp), aes(x=Prey, y=value)) + 
      geom_bar(stat="identity", color="red", position=position_dodge()) +
      # geom_line() +  geom_point()+
      labs(y = "Vulnerability +- 1*sd")+
     # facet_wrap(~ Predator, ncol=2,scales='free_y')+
      ggtitle(pp)+
      theme(axis.text.x = element_text(angle=90))+
      geom_errorbar(aes(ymin=value-1*std, ymax=value+1*std), width=.2,position=position_dodge(.9),colour = "red")
    ggsave(file.path(figures,'Vulnerability',paste0("vulnerability_",pp,".png")),width = 10, height = 7, units = "cm")
  }

  
  if (FALSE) {
    bb<-subset(b,name=='init_stl_other_suit_slope')
    a<-round(t(t(tapply(bb$CV.round,list(bb$Predator),sum))))
    print(a, na.print = "")
    
    bb<-subset(b,name=='init_season_overlap')
    print(round(t(t(tapply(bb$CV.round,list(bb$Predator),sum)))), na.print = "")
  } 
  
  ### catchability
  
  a<-Read.SMS.std(excludeNotUsed=TRUE,remove=remo)
  a<-subset(a,name=='qq_ini')
  
  b<-readFleetInfo() 
  a<-merge(x=a,y=b, all.y=TRUE)
  a<-subset(a,select=c(species,age,fleet,fleetName,value,std))
  a<-a[order(a$species,a$fleet,a$age),]
  #b<-aggregate(cbind(mv2=value)~species+fleet,data=a,mean,na.rm=TRUE)
  #a<-merge(a,b)

  for (i in (2:dim(a)[[1]])) {if (is.na(a[i,"value"])) a[i,"value"]<-a[i-1,"value"]}
  a$Species<-sp.names[a$species]
  
  for (pp in sort(unique(a$species)) ) {
    p<- ggplot(subset(a,species==pp), aes(x=age, y=value)) + 
    #geom_bar(stat="identity", color="red", position=position_dodge()) +
    geom_line() +  geom_point()+
    labs(y = "Catchability+- 2*sd",x='Age')+
    facet_wrap(~ fleetName, ncol=2,scales='free_y')+
    
    scale_x_continuous(breaks = scales::breaks_width(1),minor_breaks = NULL)+
    #  scale_x_discrete()+
    geom_errorbar(aes(ymin=value-2*std, ymax=value+2*std), width=.2,
                  position=position_dodge(.9))
    ggsave(file.path(figures,"Catchability",paste0("uncertanty_catchability_",sp.names[pp],".png")),width = 13, height = 10, units = "cm")
  }
  
  ## survey variance
  a<-Read.SMS.std(excludeNotUsed=TRUE,remove=remo)
  a<-subset(a,name=='qq_s2_ini')
  
  b<-readFleetInfo() 
  a<-merge(x=a,y=b, all.y=TRUE)
  a<-subset(a,select=c(species,age,fleet,fleetName,value,std))
  a<-a[order(a$species,a$fleet,a$age),]
  

  for (i in (2:dim(a)[[1]])) {if (is.na(a[i,"value"])) a[i,"value"]<-a[i-1,"value"]}
  a$Species<-sp.names[a$species]
  
  for (pp in sort(unique(a$species)) ) {
    p<- ggplot(subset(a,species==pp), aes(x=age, y=value)) + 
      #geom_bar(stat="identity", color="red", position=position_dodge()) +
      geom_line() +  geom_point()+
      labs(y = "Observation variance+- 2*sd",x='Age')+
      facet_wrap(~ fleetName, ncol=2,scales='free_y')+
      
      scale_x_continuous(breaks = scales::breaks_width(1),minor_breaks = NULL)+
      #  scale_x_discrete()+
      geom_errorbar(aes(ymin=value-2*std, ymax=value+2*std), width=.2)
    ggsave(file.path(figures,"Survey_observations",paste0("surv_variance_",sp.names[pp],".png")),width = 13, height = 10, units = "cm")
  }
  
}



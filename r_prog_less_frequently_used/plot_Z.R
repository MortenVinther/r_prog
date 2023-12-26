##########################################################################
do_Z<-function(data_in=data.path,label=' ',  first_cohort=1973, last_cohort=2014,inclSPn=16:22,age_ranges=NULL) {
  cleanup()


 age_z<-matrix(data=0,nrow=2,ncol=length(inclSPn))  
 colnames(age_z)<-sp.names[inclSPn]
 rownames(age_z)<-c("first_age","last_age")
 if (is.null(age_ranges)) {
   age_z[1,]<-SMS.control@avg.F.ages[sp.names[inclSPn],1]
   age_z[2,]<-SMS.control@species.info[inclSPn,1]- SMS.control@species.info[inclSPn,'+group']
  } else age_z<-age_ranges

 print(age_z)
  Init.function()
  
  dat<-Read.summary.data(read.init.function=F,dir=data_in) %>% filter(Species.n >=first.VPA & Species.n %in% inclSPn )
  head(dat)
  
  dat$DM1<-dat$M1*dat$N.bar
  dat$DM2<-dat$M2*dat$N.bar
  dat$DM<-dat$DM1+dat$DM2
  dat$DF<-dat$F*dat$N.bar
  dat$DZ<-dat$Z*dat$N.bar
  dat$sumM1<-dat$M1
  dat$sumM2<-dat$M2
  dat$sumZ<-dat$Z
  dat$sumM<-dat$sumM1+dat$sumM2
  
  dat<-aggregate(cbind(DM,DM1,DM2,DF,DZ,Z,sumM,sumM1,sumM2,sumZ,C.obs)~Year+Age+Species.n+Species,sum,na.rm=T,data=dat)
  
  dat$M<-dat$DM/dat$DZ*dat$Z
  dat$M1<-dat$DM1/dat$DZ*dat$Z
  dat$M2<-dat$DM2/dat$DZ*dat$Z
  dat$F<-dat$DF/dat$DZ*dat$Z
  
  dat$cohort<-dat$Year-dat$Age 
  dat<-as_tibble(dat) %>% dplyr::select(-DM1,-DM2,-DF,-DZ,-DM) %>%arrange(Species.n,cohort,Age)
  
  WK<-dat%>%  dplyr::select(Year,Age,Species,M,M1,M2,F,Z,sumM,sumM1,sumM2,sumZ) %>% mutate(run=label) %>% arrange(run,Species,Year,Age)
  write_csv(WK,file=file.path(data_in,'WKBBALTPEL.csv'))
  
  #dat<- dat %>% arrange(Species.n,Species,cohort,Age) %>% 
  #    mutate(catch_Z=log(lead(C.obs)/C.obs),use=cohort==lead(cohort)) %>%
   #   filter(use & cohort>=first_cohort & cohort<=last_cohort) %>%dplyr::select(-use)
 
  
  dat<- dat %>% arrange(Species.n,Species,cohort,Age) %>% 
    mutate(catch_Z=log(lead(C.obs)/C.obs),use=cohort==lead(cohort)) %>%
   filter(cohort>=first_cohort & cohort<=last_cohort) 
  
  
  do_z_plot <-function(spNo=2) {
    spName<-sp.names[spNo]
    cat('doing ',spName,'\n')
    d2<-filter(dat,Species.n==spNo & Age>0 & C.obs>0)
    
    first_age_z<-age_z[1,spName]
    last_age_z<-age_z[2,spName]
    n_plot<-last_age_z-first_age_z+1
    
    ra<-range(log(d2$C.obs))
    if (ra[2]-ra[1] >5) ra_step<-2 else ra_step<-1
    ra<-seq(floor(ra[1]),ceiling(ra[2]+6)+1,ra_step)
    png(file.path(data_in,paste0('Z_histo_',spName,'.png')),width = 480, height = 480, units = "px",pointsize=12)
    z<-filter(d2,Age>=first_age_z & Age <last_age_z & use) ;median_Z<-median(z$Z,na.rm=T); hist(z$Z,main=spName,xlab='Z'); abline(v=median_Z, col='red',lwd=2) 
    
    p<- ggplot(d2, aes(x=Age, y=log(C.obs)),lwd=5) + geom_point()+geom_line()+
        facet_wrap(vars(cohort))+theme_bw()+ 
        theme(panel.grid.major.y = element_blank(),  panel.grid.minor.y = element_blank())+
        ylab('log(catch numbers)')+labs(title=paste(spName, "median Z (ages",first_age_z,"-",last_age_z-1,") =",round(median_Z,2),"   Slope of red lines = ",-round(median_Z,1)))+
        lapply(ra,function(x)geom_abline(intercept=x,slope=-median_Z,col='red',lty=3,lwd=0.75))
    png(file.path(data_in,paste0('catch_curve_',spName,'.png')),width = 900, height = 1000, units = "px",pointsize=16)
    print(p)
    
     p<- ggplot(filter(d2,Age>=first_age_z & Age <last_age_z & use), aes(x=Year, y=-catch_Z)) + geom_point()+
       geom_line()+geom_smooth(method="loess",formula=y ~ x)+
        facet_wrap(vars(paste('Age',Age)))+theme_bw()+ 
        ylab('Z from catch curve analysis')
     if (n_plot <4) {
       png(file.path(data_in,paste0('catch_curve_Z_',spName,'.png')),width = 900, height = 500, units = "px",pointsize=16)
     } else png(paste0('catch_curve_Z_',spName,'.png'),width = 900, height = 900, units = "px",pointsize=16)
     
      print(p)
   
    d3<- d2 %>%  pivot_longer(cols=c("M","M1","M2","F","Z"),names_to='Mortality')
    
     
    p<- ggplot(filter(d3,Mortality %in% c('M1','M2','F') & Age>=first_age_z & Age<last_age_z), aes(x=Year, y=value, fill=Mortality)) +
       geom_bar(stat="identity")  +
       geom_line(aes(x=Year, y=-catch_Z),stat="identity",lwd=0.75)+
       facet_wrap(vars(paste('Age',Age)))+theme_bw()+ 
       theme(panel.grid.major.y = element_blank(),  panel.grid.minor.y = element_blank())+
       ylab('Mortality')
    
    if (n_plot <4) {
      png(file.path(data_in,paste0('Z_',spName,'.png')),width = 900, height = 500, units = "px",pointsize=16)
    } else png(file.path(data_in,paste0('Z_',spName,'.png')),width = 900, height = 900, units = "px",pointsize=16)
    
    print(p)
    cleanup() 
       
  }
  for (i in (inclSPn)) do_z_plot(i)
}


if (FALSE) { #test
  #inclSPn<-16:17
  inclSPn<-2:3
  age_z<-matrix(data=0,nrow=2,ncol=length(inclSPn))  
  colnames(age_z)<-as.character(inclSPn)
  rownames(age_z)<-c("first_age","last_age")
  age_z[1,]<-SMS.control@avg.F.ages[sp.names[inclSPn],1]
  age_z[2,]<-SMS.control@species.info[inclSPn,1]- SMS.control@species.info[inclSPn,'+group']
  age_z
  do_Z(data_in=data.path,label=' ',  first_cohort=1973, last_cohort=2014,inclSPn=16:27,age_ranges=NULL) 
}

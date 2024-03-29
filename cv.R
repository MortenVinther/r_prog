

a<-Read.SMS.std()
gw<-1200
gh<-800


aa<-subset(a,name=='her_M2')
head(aa)
round(tapply(aa$CV.round,list(aa$year,paste(aa$species,sp.names[aa$species],sep='.')),sum,na.rm=T),1)
trellis.device(device='png',file=file.path(data.path,paste0('CV_M2','.png')),width = gw, height = gh)
print(xyplot(CV.round~year|paste0('Age:', aa$age),data=aa,type=c('a','g'),
             panel = function(x,y,...){
               #   panel.refline(h = seq(10,50,10))
               panel.xyplot(x,y,...)
             },
             scales = list(x=list(relation="same"),
                           y=list(relation='same')),
             xlab='year', ylab='CV of M2 (%)',lwd=2,
             auto.key =   list(space = "right", points = FALSE, lines = TRUE)))

cleanup() 

aa<-subset(a,name=='her_M2')

aaa<-rbind(aa %>% transmute(species,year,age,value,source='mean'),
      aa %>% transmute(species,year,age,value=value+2*std,source='+2*std'),
      aa %>% transmute(species,year,age,value=value-2*std,source='-2*std'))
png(filename = 'M2_std.png',width = 900, height = 900, units = "px", pointsize = 18,)
aaa<-aaa %>%  mutate(source=factor(source,levels=c('+2*std','mean','-2*std')))
ap<-  ggplot(data=aaa, aes(x=year, y=value, group=source)) +
  geom_line(aes(linetype=source))+
    scale_linetype_manual(values=c("dotted","solid", "dotted"))+
 
  facet_wrap(vars(paste0('Age:',age)),ncol=2,scales="free_y")+
    ylim(0,NA)  +ylab('M2')+
    theme_bw() 
print(ap)
cleanup()

aa<-subset(a,name=='her_M2')
aa<-filter(aa,year==2022)
aaa<-rbind(aa %>% transmute(species,year,age,value,source='mean'),
           aa %>% transmute(species,year,age,value=value+2*std,source='+2*std'),
           aa %>% transmute(species,year,age,value=value-2*std,source='-2*std'))

png(filename = 'M2_std_2022.png',width = 500, height = 500, units = "px", pointsize = 18,)
aaa<-aaa %>%  mutate(source=factor(source,levels=c('+2*std','mean','-2*std')))
ap<-  ggplot(data=aaa, aes(x=age, y=value, group=source)) +
  geom_line(aes(linetype=source))+
  scale_linetype_manual(values=c("dotted","solid", "dotted"))+
  
  ylim(0,NA)  +ylab('M2')+
  theme_bw() 
print(ap)
cleanup()


aa<-subset(a,name=='her_avgM2')

aaa<-rbind(aa %>% transmute(species,age,value,source='mean'),
           aa %>% transmute(species,age,value=value+2*std,source='+2*std'),
           aa %>% transmute(species,age,value=value-2*std,source='-2*std'))
png(filename = 'M2_std_avg.png',width = 500, height = 500, units = "px", pointsize = 18,)
aaa<-aaa %>%  mutate(source=factor(source,levels=c('+2*std','mean','-2*std')))
ap<-  ggplot(data=aaa, aes(x=age, y=value, group=source)) +
  geom_line(aes(linetype=source))+
  scale_linetype_manual(values=c("dotted","solid", "dotted"))+

  ylim(0,NA)  +ylab('M2')+
  theme_bw() 
print(ap)
cleanup()




aa<-subset(a,name=='hist_SSB')
round(tapply(aa$CV.round,list(aa$year,paste(aa$species,sp.names[aa$species],sep='.')),sum,na.rm=T),1)
trellis.device(device='png',file=file.path(data.path,paste0('CV_SSB','.png')),width = gw, height = gh)
print(xyplot(CV.round~year|paste(aa$species,sp.names[aa$species],sep='.'),data=aa,type=c('a','g'),
             panel = function(x,y,...){
            #   panel.refline(h = seq(10,50,10))
               panel.xyplot(x,y,...)
             },
             scales = list(x=list(relation="same"),
                           y=list(relation='same')),
             xlab='year', ylab='CV of SSB (%)',lwd=2,
             auto.key =   list(space = "right", points = FALSE, lines = TRUE)))
   
cleanup() 


aa<-subset(a,name=='avg_F' & species>0)
round(tapply(aa$CV.round,list(aa$year,paste(aa$species,sp.names[aa$species],sep='.')),sum,na.rm=T),1)

trellis.device(device='png',file=file.path(data.path,paste0('CV_avg_F','.png')),width = gw, height = gh)
print(xyplot(CV.round~year|paste(aa$species,sp.names[aa$species],sep='.'),data=aa,type=c('a','g'),
             panel = function(x,y,...){
               #  panel.refline(h = seq(10,100,10))
               panel.xyplot(x,y,...)
             },
             scales = list(x=list(relation="same"),
                           y=list(relation='same')),
             xlab='year', ylab='CV of mean F (%)',lwd=2,
             auto.key =   list(space = "right", points = FALSE, lines = TRUE)))

cleanup() 



aa<-subset(a,name=='log_rec') # used std as CV as it is the log recruitment
aa[aa$std>=0.5,'std']<-0.5
round(tapply(aa$std*100,list(aa$year,paste(aa$species,sp.names[aa$species],sep='.')),sum,na.rm=T),1)
trellis.device(device='png',file=file.path(data.path,paste0('CV_Recruit','.png')),width = gw, height = gh)
print(xyplot(std*100~year|paste(aa$species,sp.names[aa$species],sep='.'),data=aa,type=c('a','g'),
             panel = function(x,y,...){
               #   panel.refline(h = seq(10,50,10))
               panel.xyplot(x,y,...)
             },
             scales = list(x=list(relation="same"),
                           y=list(relation='same')),
             xlab='year', ylab='CV of recruitmet (%) truncated at 50%',lwd=2,
             auto.key =   list(space = "right", points = FALSE, lines = TRUE)))

cleanup() 


aa<-subset(a,name=='M2_sd0')
if (dim(aa)[[1]]>0) {
  round(tapply(aa$CV.round,list(aa$year,paste(aa$species,sp.names[aa$species],sep='.')),sum,na.rm=T),1)
  trellis.device(device='png',file=file.path(data.path,paste0('CV_M2_age_0','.png')),width = gw, height = gh)
  print(xyplot(CV.round~year|paste(aa$species,sp.names[aa$species],sep='.'),data=aa,type=c('a','g'),
               panel = function(x,y,...){
             #    panel.refline(h = seq(5,30,5))
                 panel.xyplot(x,y,...)
               },
               scales = list(x=list(relation="same"),
                             y=list(relation='same')),
               xlab='year', ylab='CV of M2 for age 0 (%)',lwd=2,
               auto.key =   list(space = "right", points = FALSE, lines = TRUE)))
  
  cleanup() 
}

aa<-subset(a,name=='M2_sd1')
if (dim(aa)[[1]]>0) {
  round(tapply(aa$CV.round,list(aa$year,paste(aa$species,sp.names[aa$species],sep='.')),sum,na.rm=T),1)
  trellis.device(device='png',file=file.path(data.path,paste0('CV_M2_age_1','.png')),width = gw, height = gh)
  print(xyplot(CV.round~year|paste(aa$species,sp.names[aa$species],sep='.'),data=aa,type=c('a','g'),
               panel = function(x,y,...){
               # panel.refline(h = seq(2,30,2))
                 panel.xyplot(x,y,...)
               },
               scales = list(x=list(relation="same"),
                             y=list(relation='same')),
               xlab='year', ylab='CV of M2 for age 1 (%)',lwd=2,
               auto.key =   list(space = "right", points = FALSE, lines = TRUE)))
  
  cleanup() 
}


aa<-subset(a,name=='M2_sd2')
if (dim(aa)[[1]]>0) {
  round(tapply(aa$CV.round,list(aa$year,paste(aa$species,sp.names[aa$species],sep='.')),sum,na.rm=T),1)
  trellis.device(device='png',file=file.path(data.path,paste0('CV_M2_age_2','.png')),width = gw, height = gh)
  print(xyplot(CV.round~year|paste(aa$species,sp.names[aa$species],sep='.'),data=aa,type=c('a','g'),
               panel = function(x,y,...){
               #  panel.refline(h = seq(2,30,2))
                 panel.xyplot(x,y,...)
               },
               scales = list(x=list(relation="same"),
                             y=list(relation='same')),
               xlab='year', ylab='CV of M2 for age 1 (%)',lwd=2,
               auto.key =   list(space = "right", points = FALSE, lines = TRUE)))
  
  cleanup() 
}


#############################################################################
# CV of SSB, F or recruitment, calculated from mean and standard deviations

plot.CV<-function(paper=FALSE,nox=1,noy=1,var.name='hist_SSB',legends=T,legend.x=NULL,legend.y=NULL,ylim=NULL) {
tit<-F
tmp<-Read.SMS.std()

if (var.name=='hist_SSB') {in.var.name<- "hist_SSB"; titl<-"CV of SSB";}
if (var.name=='F') {in.var.name<- "avg_F"; titl<-"CV of mean F";}
if (var.name=='Rec') {in.var.name<- "rec_sd"; titl<-"CV of recruitment";}
if (var.name=='M2_sd0') {in.var.name<- 'M2_sd0'; titl<-"CV of M2 age 0";}
if (var.name=='M2_sd1') {in.var.name<- 'M2_sd1'; titl<-"CV of M2 age 1";}
if (var.name=='M2_sd2') {in.var.name<- 'M2_sd2'; titl<-"CV of M2 age 2";}

a<-droplevels(subset(tmp,name==in.var.name,drop=TRUE))
if (is.null(a$species)) a$species<-1
a<-data.frame(Species=sp.names[a$species],Species.no=a$species,Year=a$year,
               value=a$value,std=a$std,CV=a$std/a$value*100)

CV.obs<-tapply(a$CV,list(a$Year,a$Species),sum,na.rm=T)
print(t(round(t(CV.obs),1)))

if (paper) dev<-"png" else dev<-"screen"
noxy<-nox*noy

year.obs<-as.numeric(dimnames(CV.obs)[[1]])
sp.nam<-dimnames(CV.obs)[[2]]
max.CV.obs<-max(CV.obs,na.rm=T)
spno<-sort(unique(a$Species.no))
sp.indx<-0
for (sp in spno) {
  sp.indx<-sp.indx+1
  #if (paper) { 
    if (sp==spno[1]) {                         
    newplot(dev,nox,noy,filename=paste("CV_",in.var.name,sep=''),Portrait=T,w8=8,w11=9);
      plot(year.obs,CV.obs[,sp.indx],xlab=" ",ylab="CV (%)",type='b',ylim=c(0,ifelse(is.null(ylim),max.CV.obs,ylim)),lty=sp.indx,col=sp.indx,pch=sp.indx,lwd=2 )
      if (legends) legend(x=ifelse(is.null(legend.x),min(year.obs),legend.x),y=ifelse(is.null(legend.y),1,legend.y),
                           sp.nam,col=1:length(spno),ncol=4,pch=1:nsp,text.col=1,lwd=2,lty=1:length(spno))
    }
    else lines(year.obs,CV.obs[,sp.indx],col=sp.indx,type='b',pch=sp.indx,lwd=2,lty=sp.indx)
  #}
}
if (paper) cleanup()
}

if (FALSE) {
  plot.CV(paper=T,nox=1,noy=1,var.name='hist_SSB',legends=T,legend.y=39)
  plot.CV(paper=T,nox=1,noy=1,var.name='F',legends=T,legend.y=55)
  plot.CV(paper=T,nox=1,noy=1,var.name='Rec',legends=T,legend.y=150)
  plot.CV(paper=T,nox=1,noy=1,var.name='Rec',legends=T,legend.y=5,ylim=35)
  
  plot.CV(paper=T,nox=1,noy=1,var.name='M2_sd0',legend.y=30)
  plot.CV(paper=T,nox=1,noy=1,var.name='M2_sd1')
  plot.CV(paper=T,nox=1,noy=1,var.name='M2_sd2')
}


#############################################################################
# CV of SSB and F bar on the same graph. Only really works for single species runs
#
# needs an update!!!!
if (FALSE) {
l.model.y<-SMS.control@last.year.model
a<-read.table(file.path(data.path,"sms.std"),comment.char = "#",header=FALSE,skip=1) 
tmp<-data.frame(index=a$V1,name=a$V2, value=a$V3,  std=a$V4,CV=a$V4/a$V3*100 )

obs<-subset(tmp,name=="hist_SSB",drop=TRUE)
yobs<-nrow(obs)/(nsp-first.VPA+1)

a<-data.frame(Species=rep(sp.names[first.VPA:length(sp.names)],each=yobs),
              Year=rep(seq(l.model.y-yobs+1,l.model.y),(nsp-first.VPA+1)),
               value=obs$value,std=obs$std,CV=obs$CV)

CV.ssb<-tapply(a$CV,list(a$Year,a$Species),sum,na.rm=T)
year.ssb<-as.numeric(dimnames(CV.ssb)[[1]])

obs<-subset(tmp,name=="avg_F",drop=TRUE)
yobs<-nrow(obs)/(nsp-first.VPA+1)

a<-data.frame(Species=rep(sp.names[first.VPA:length(sp.names)],each=yobs),
              Year=rep(seq(l.model.y-yobs+1,l.model.y),(nsp-first.VPA+1)),
               value=obs$value,std=obs$std,CV=obs$CV)

CV.F<-tapply(a$CV,list(a$Year,a$Species),sum,na.rm=T)
year.F<-as.numeric(dimnames(CV.F)[[1]])
nox<-1; noy<-2;
paper<-FALSE

if (paper) dev<-"wmf" else dev<-"screen"
noxy<-nox*noy

#Maximum observed CV for setting graph limits
max.CV <-max(CV.F,CV.ssb)
if(nsp>1) {
        warning("Can't plot SSB and F on same graph for more than one species.")
   } else {
        newplot(dev,nox,noy,filename="CV_F,SSB");
        par(mar=c(4,4,2,3))
        plot(range(year.ssb,year.F),c(0,max.CV),xlab="",ylab="CV (%)",type='n',ylim=c(0,max.CV))
        lines(year.F,CV.F[,1],xlab=" ",ylab="CV (%)",type='b',ylim=c(0,max.CV),lty=1,col=1,pch=1,lwd=2 )
        lines(year.ssb,CV.ssb[,1],type="b",lty=1,col=1,pch=2,lwd=2 )
        legend("topleft",legend=c("F.bar","SSB"),col=1,ncol=2,pch=1:2,text.col=1,lwd=2,lty=1)
        grid()
}
if (paper) cleanup()

N<-subset(tmp,name=="term_N",drop=TRUE)
N.next<-subset(tmp,name=="term2_N",drop=TRUE)
N.next[1,'CV']<-NA

if (dim(N)[[1]]>0 | dim(N.next)[[1]]>0) {  
  control<-read.FLSMS.control()
  fa<-control@first.age
  la<-control@species.info[1,'last-age']
  
  a<-data.frame(age=seq(fa,la),N=N$value,          CV.N=N$CV)
  b<-data.frame(age=seq(fa+1,la),N.next=N.next$value,CV.N.next=N.next$CV)
  ab<-merge(a,b)
  max.CV<-max(ab$CV.N,ab$CV.N.next,na.rm=T)
  plot(ab$age,ab$CV.N,xlab="age",ylab="CV (%)",type='n',ylim=c(0,max.CV))
  lines(ab$age,ab$CV.N,xlab=" ",ylab="CV (%)",type='b',ylim=c(0,max.CV),lty=1,col=1,pch=1,lwd=2 )
  lines(ab$age,ab$CV.N.next,type="b",lty=1,col=1,pch=2,lwd=2 )
  legend("top",legend=c("N 1. Jan 2008","N 1. Jan 2009"),col=1,ncol=2,pch=1:2,text.col=1,lwd=2,lty=1)
  grid()
  if (paper) cleanup()
}

}
########################################################################



plot.CV<-function(paper=FALSE,nox=1,noy=1,var.name='hist_SSB',newPlot=T,tit=T,outfile='Uncertanties',w8=8,w11=9,pointsize=12) {

  if (paper) dev<-"png" else dev<-"screen"
  noxy<-nox*noy
  
  if (plot.ref.lines) ref<-Read.reference.points()

  tmp<-Read.SMS.std()
  tmp$name[tmp$name=="next_SSB"]<-'hist_SSB'
  a<-subset(tmp,name==var.name & species>0 ,drop=TRUE)
  a$Species<-sp.names[a$species]
  a<-subset(a,Species!="Saithe",drop=TRUE)
  
  if (var.name=='hist_SSB')  ytitl<-"SSB (1000 t)"  else if (var.name=='avg_F') ytitl<-"mean F"  else if (var.name=='rec_sd') ytitl<-"Recruitment (10^6)" else if (var.name=='M2_sd0' | var.name=='M2_sd1' | var.name=='M2_sd2') ytitl<-"M2"  else stop("name must be hist_SSB or avg_F")
  
  if (var.name=='rec_sd') {
      a$value<-a$value/1000
      a$std<-a$std/1000
  }
  if (var.name=='hist_SSB') {
      a$value<-a$value/1000
      a$std<-a$std/1000
  }
# tapply(a$value,list(a$species,a$year),sum)    
  if (newPlot) newplot(dev,nox,noy,filename=outfile,w8=w8,w11=w11,Portrait=TRUE,pointsize=pointsize);
       par(mar=c(3,5,1,2))  #c(bottom, left, top, right)
 
  by(a,list(a$species),function(x) {
      minval<-min(x$value-2*x$std*2,0)
      maxval<-max(x$value+2*x$std*2)
      plot( x$year,x$value,xlab='',ylab=ytitl,ylim=c(minval,maxval),type='b')
      if (tit) title(main=x[1,]$Species)

       # print(aaa[2,'Cod',,1])
      lines(x$year,x$value-2*x$std,lty=2)
      lines(x$year,x$value+2*x$std,lty=2)
     # lines(x$year,aaa[2,x[1,"Species"],as.character(1975:2000),2],col=2,lty=5) # from plot_M2_uncertanty_MCMC
     # lines(x$year,aaa[3,x[1,"Species"],as.character(1975:2000),2],col=2,lty=5) # from plot_M2_uncertanty_MCMC
     if (plot.ref.lines) {
       sp<-x$species[1]
       if (var.name=='avg_F') {
        if (ref[sp,"Flim"]>0) abline(h=ref[sp,"Flim"],lty=2,lwd=2)
        if (ref[sp,"Fpa"]>0) abline(h=ref[sp,"Fpa"],lty=3,lwd=2)
       }
       if (var.name=='hist_SSB') {
          Blim<-ref[sp,"Blim"]/1000; Bpa<-ref[sp,"Bpa"]/1000
          if (Blim>0) abline(h=Blim,lty=2,lwd=2)
          if (Bpa>0) abline(h=Bpa,lty=3,lwd=2)
       }
     }  

  })  
  if (dev!="screen") cleanup();
}


#
#plot.CV(paper=FALSE,nox=2,noy=2,var.name='M2_sd0') 

plot.ref.lines<-F
plot.CV(paper=TRUE,nox=2,noy=4,var.name='M2_sd0',outfile='Uncertanties_M2_age0')
plot.CV(paper=TRUE,nox=2,noy=4,var.name='M2_sd1',outfile='Uncertanties_M2_age1')
plot.CV(paper=TRUE,nox=2,noy=4,var.name='M2_sd2',outfile='Uncertanties_M2_age2')

plot.CV(paper=T,nox=1,noy=3,var.name='avg_F',newPlot=T,tit=T)
plot.CV(paper=T,nox=1,noy=3,var.name='hist_SSB',newPlot=F,tit=F)
plot.CV(paper=T,nox=1,noy=3,var.name='rec_sd',newPlot=F,tit=F)


####################################
if (FALSE) {
  tmp<-Read.SMS.std()
  
  a<-subset(tmp,name %in% c("hist_SSB","avg_F","rec_sd") ,drop=TRUE)
  a<-data.frame(Species=sp.names[a$species],Year=a$year, variable=a$name,
                 value=a$value,std=a$std,CV=a$std/a$value*100)
                 
                 
  a$titl<- ifelse (a$variable=='hist_SSB',"SSB", ifelse(a$variable=='avg_F',"Average F",ifelse(a$variable=='rec_sd', "Recruitment",'error')))
  
  pp<-function(x){
    cleanup()
    trellis.device(device = "windows",
                   color = T, width=6, height=8,pointsize = 10,
                   new = TRUE, retain = FALSE)
  
    xyplot( CV~Year|titl,groups=Species, data=x,
      type='b',lwd=1 , layout=c(1,3), ylab='Coefficient of Variation (%)',
       between = list(y = c(1, 1),x = c(1, 1)),
       strip = strip.custom( bg='white'),par.strip.text=list(cex=0.9, lines=1.7),
       auto.key = list(space = "bottom", points = T, lines = F,cex=0.9, columns = 3,title='Species') ,
  
        scales = list(x = list( cex=0.8,relation='same'), y= list(cex=0.8),alternating = 1,relation='free')
    )
  }
  
  pp(a)
  
  
  pp(droplevels(subset(a,Year>=1985 & (Species %in% c('Cod','Haddock','Whiting','Saithe')))))
  pp(droplevels(subset(a,Year>=1985 & !(Species %in% c('Cod','Haddock','Whiting','Saithe','Plaice','Sole')))))
  
  pp(droplevels(subset(a,Year>=1985 & !(Species %in% c('Cod','Haddock','Whiting','Saithe','Plaice','Sole')))))
}


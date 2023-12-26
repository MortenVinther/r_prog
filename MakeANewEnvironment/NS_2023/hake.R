# Hake
species<-'HKE'
sp<-'Hake'
fdir<-file.path( root,"Data_NorthSea","input_NS_2023",'ByStock',species)
setwd(fdir)

load("ss3_R_output.RData",verbose=TRUE) # from WGBIE Sharepoint
summary(output)
myreplist<-output; rm(output)

#output$natage

str(myreplist,1,list.len=200)
nl<-myreplist$lbinspop    # length bins used for population
write.csv(nl,file='lengthBin.csv',row.names =FALSE)


convL<-function(x,v.names='n') {
  # test:  x<-N
  x<-subset(x,x[,"Beg/Mid"]=='M')
  x<-subset(x,select=c("Yr","Seas",as.character(nl)))

  x<-reshape(x,direction='long',varying=list(3:(2+length(nl))),v.names=v.names)
  x$Length<-nl[x$time]
  x<-subset(x,select=c("Yr","Seas",'Length',v.names))
  return(x)
}



# Data for SMS (all by length classess)
str(myreplist$inputs)
str(myreplist,1,list.len=300)
myreplist$Growth_Parameters
nl<-myreplist$lbinspop    # length bins used for population
write.csv(nl,file='lengthBin.csv',row.names =FALSE)

#myreplist$mean_body_wt

convL<-function(x,v.names='n') {
  # test:  x<-N
  # there are three BirthSeas ?

  x<-subset(x,x[,"Beg/Mid"]=='M')

  x<-subset(x,select=c("Yr","Seas",as.character(nl)))

  x<-reshape(x,direction='long',varying=list(3:(2+length(nl))),v.names=v.names)
  x$Length<-nl[x$time]
  x<-subset(x,select=c("Yr","Seas",'Length',v.names))
  return(x)
}

subset(myreplist$natlen,Yr==1998 & Seas==1)

# head(myreplist$natlen)
# summary(myreplist$natlen)
n<-convL(x=myreplist$natlen,v.names='n')
head(n)
n$w<-myreplist$Growth_Parameters[1,'WtLen1']*n$Length**myreplist$Growth_Parameters[1,'WtLen2']

head(subset(n,Length==50 & Yr==1998) ,20)


# tst<-subset(n,Length==50 & Seas==1)  # there are 2 BirthSeas and 2 sex'es ?
# subset(tst,Yr==1998)

n$age<-cut(n$Length,breaks=c(0,25,30,35,40,50,60,70,80,100,150),include.lowest =TRUE)
#n$age<-ifelse(n$Length<25,1,ifelse(n$Length<60,2,3))

n$nl<-n$Length*n$n
n$nw<-n$w*n$n
head(subset(n,Length==50))

nn<-aggregate(cbind(n,nl,nw) ~ Yr+Seas+age,data=n,sum)
nn$ml<-nn$nl/nn$n
nn$mw<-nn$nw/nn$n
nn$bio<-nn$n*nn$mw
head(nn)

if (FALSE) {
  tapply(nn$ml,list(nn$Yr,nn$Seas,nn$age),sum)
  tapply(nn$mw,list(nn$Yr,nn$Seas,nn$age),sum)
  round(tapply(nn$bio,list(nn$Yr,nn$Seas),sum),0)
}

#write.csv(nb,file='SS3_results.csv',row.names =FALSE)


#proportion in the North Sea
a<-read.csv(file=file.path(root,exchangeDir,'ByStock',species,"Hake_catch_WGBIE_2023.csv"),skip=3)
a<- filter(a,year>=1974)
a[is.na(a)]<-0
a
tst<-cbind(a %>% mutate(NorthSea=(L_4)/LTotal) %>% select(year,NorthSea),a %>% mutate(NorthSea=(L_4+L_3)/LTotal) %>%
             select(NSK=NorthSea)) %>% filter(year>=2013)
tst

#proportion i area 4 of landings in 3456 combined, data from 2013
p4<-filter(a,year>=2013) %>% transmute(year=year,p4=L_4/(L_3+L_4+L_5+L_6)) %>% summarize(p4=mean(p4))
p4
a<-mutate(a,L_4=if_else(year>=2013,L_4,L_3456*as.numeric(p4))) %>%
  transmute(year=year,NorthSea=L_4/LTotal)

a
agem<- a %>% select(year,NorthSea)
cleanup()

filename=file.path(paste0(sp,'_prop_NS'))
newplot(dev='png',nox=1,noy=1,Portrait=TRUE,filename=filename,dir=fdir,w11=6,w8=8);
#par(mar=c(3,5,3,2))
plot(x=a$year,y=a$NorthSea,ylim=c(0,0.3),ylab='Proportion within the North Sea',xlab='Year',type='b')
cleanup()

load("Table_2.RData",verbose=TRUE)  # from WGBIE Sharepoint 2023_hke.27.3a46-8abd_assessment_data.zip

names(tab2)
sort(unique(tab2$area))
sort(unique(tab2$season))
sort(unique(tab2$seasontype))
sort(unique(tab2$category))

tab2[tab2$seasontype=='Year','season']<-9
ns<-c( "27.4","27.4.a","27.4.b","27.4.c" )
tab2$area2<-'outside'
tab2[tab2$area %in% ns,'area2']<-"North Sea"

round(xtabs(canum*weca~ year+area2,data=tab2)/1000000000,1)


tab3<-filter(tab2, category=='Landings' & season !=9  ) %>%
    mutate(caton=canum*weca/1000000) %>% group_by(year,season,area2) %>% summarize(caton=sum(caton,na.rm=TRUE)) %>%ungroup()
tab3

tab4<-tab3 %>% group_by(year,season) %>% mutate(perc=caton/sum(caton))
b<-xtabs(perc~year+season,data=tab4)
round(b*100,1)
tst2<-tab4 %>% group_by(year,area2) %>% summarize(caton=sum(caton)) %>% group_by(year) %>% mutate(perc=caton/sum(caton)) %>% ungroup() %>% filter(area2=='North Sea')
tst2
b<-xtabs(perc~year+season,data=filter(tab4,area2=='North Sea'))
round(b*100,1)

per1<-as.data.frame(b)

c22<-read_delim(file="intercatch_2022.txt",delim='\t')
sort(unique(c22$Area))
c22$area2<-'outside'
c22[c22$Area %in% ns,'area2']<-"North Sea"
c22

tab33<-filter(c22, CatchCategory=='Landings' & Season %in% (1:4) & area2=='North Sea') %>% group_by(Year,Season) %>% summarize(caton=sum(CATON,na.rm=TRUE)) %>%ungroup()
tab33
sum(tab33$caton)/1000
tab44<-tab33 %>% group_by(Year) %>% mutate(perc=caton/sum(caton))
tab44
b<-xtabs(perc~Year+Season,data=tab44)
round(b*100,1)
per2<-as.data.frame(b) %>% rename(year=Year,season=Season)

pp<-rbind(per1,per2) %>% rename(percent=Freq)
pp$year<-as.character(pp$year)

ppp<-xtabs(percent~year+season,data=pp )
ppp
pgem<-data.frame(ppp)
ppp<-cbind(ppp,percent_NSK=subset(tst,year>=2013)$NSK)

ppp<-rbind(ppp,Average=colMeans(ppp))
round(ppp*100,1)

sc<-c(ppp['Average',1:4]/ppp['Average',5],1)
ppp<-rbind(ppp,scaled=sc)
xtab(ppp*100,caption='Proportion af annual catch',file='Q_dist.html',dec=rep(1,ncol(ppp)))

sc[1:4]/100


n<- select(nn,year=Yr,quarter=Seas,age,n,mw,ml)
str(n)
str(pgem)

pp<-pgem %>% mutate(year=as.integer(as.character(year)),quarter=as.integer(as.character(season))) %>%rename(perc=Freq)
head(pp)
scgem<-data.frame(quarter=1L:4L,prop=sc[1:4]/100)
b<-expand.grid(year=1974L:2012L,quarter=1L:4L)
b<-left_join(left_join(b,scgem),agem) %>%mutate(factor=prop*NorthSea*100)

n1<-left_join(x=filter(n,year<2013),y=b,by = join_by(year, quarter)) %>% mutate(n=n*factor) %>% select(year,quarter,age,n,mw,ml)
n2<-left_join(x=filter(n,year>=2013),y=pp, by = join_by(year, quarter)) %>% mutate(n=n*perc)  %>% select(year,quarter,age,n,mw,ml)

a<-filter(n1,year==1978)
n0<-NULL
for (y in (1974:1977)) {
  a$year<-y
  n0<-rbind(a,n0)
}

hake_n<-rbind(n0,n1,n2) %>% mutate(bio=n*mw)

# change into "ages", but delete the 0-25 cm length group
hake_n$age<-unclass(hake_n$age)-1
hake_n<-subset(hake_n,age>=1)
summary(hake_n)

ohake <- hake_n %>% transmute(year=year,quarter=quarter,species=species,age=age,sub_area=1,WSEA=mw,N=n)
write.csv(ohake,file='HAKE_N_WSEA.csv',row.names=FALSE)

ohake <- hake_n %>% transmute(year=year,quarter=quarter,species=species,age=age,sub_area=1,mean_l=ml,SMS_area=1)
write.csv(ohake,file='HAKE_length.csv',row.names=FALSE)


pl<- hake_n %>% group_by(year,quarter) %>% summarize(BIO=sum(bio))

cleanup()
#dev<-"print"
dev<-'png'

filename=file.path(paste0(sp,'_BIO.png'))
newplot(dev,nox=2,noy=2,Portrait=TRUE,filename=filename,dir=fdir,w11=8);
#par(mar=c(3,5,3,2))

by(pl,list(pl$quarter),function(x){
  b<- tapply(x$BIO/1000,list(x$year),sum)
  barplot(b,space=0.0,ylab='Biomass (1000 tonnes)')
  title(main=paste(species,' Q:', x[1,'quarter'],sep=''))

})
cleanup()


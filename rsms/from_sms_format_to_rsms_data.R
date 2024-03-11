From_SMS_format_to_rsms<-function(otherPredExist=TRUE,catchMultiplier=1,dir=data.path,sms.dat="rsms.dat") {
la<-SMS.control@max.age.all
fa<-SMS.control@first.age
rec.season <-SMS.control@rec.season
years<-c(1,1)
years[1]<-SMS.control@first.year
years[2]<-SMS.control@last.year
ny<-years[2]-years[1]+1
npr<-sum(SMS.control@species.info[,'predator']>=1)
nsp<-SMS.control@no.species
nq<-SMS.control@last.season
noAreas<-SMS.control@no.areas

#############  catch data
scanData<-function(name='canum.in'){
  cat('Reading ',name,'\n')
  x<-scan(file.path(dir,name),comment.char='#',quiet=TRUE)
  stopifnot(tail(x,1)== -999)
  x<-head(x,-1)
}
CATCHN<-scanData('canum.in')
WCATCH<-scanData('weca.in')
Prop.landed<-scanData('proportion_landed.in'); Prop.landed<-Prop.landed[1:length(WCATCH)]

b<-expand.grid(sub_area=1:noAreas,species.n=first.VPA:nsp,year=years[1]:years[2],quarter=1:nq,age=fa:la)
b<-b[order(b$sub_area,b$species.n,b$year,b$quarter,b$age),]
b<-data.frame(b,CATCHN=CATCHN,WCATCH=WCATCH,PROP_CAT=Prop.landed)
b<-subset(b,select=c(year,species.n,quarter,sub_area,age,WCATCH,CATCHN,PROP_CAT))
b$CATCHN<-catchMultiplier*b$CATCHN

pf<-Read.summary.data(dir=dir) %>%  filter((Age>0 |Quarter>2) & Species.n>=first.VPA) %>% select(Species.n,Year,Quarter,Age,"F") %>% rename(FF="F") %>%
  group_by(Species.n, Year,Age)%>%  mutate(propF=FF/sum(FF),FF=NULL) %>% as_tibble() %>% mutate(propF=if_else(is.na(propF),0,propF)) %>%
  rename(year=Year,species.n=Species.n,age=Age,quarter=Quarter)

b<-left_join(b,pf,by = join_by(year, species.n, quarter, age))
b[is.na(b$propF),'propF']<-0

out<-list(catch=b)
############## bio data
WSEA<-scanData('west.in');WSEA<-WSEA[((first.VPA-1)*noAreas*ny*(la-fa+1)*nq+1):length(WSEA)]
PROPMAT<-scanData('propmat.in')
#M<-scanData('natmor.in')
M<-scanData('natmorm1m2.out')
M1<-scanData('natmor1.in')
PROP_M2<-scanData('n_proportion_m2.in')

b<-expand.grid(sub_area=1:noAreas,species.n=first.VPA:nsp,year=years[1]:(years[2]),quarter=1:nq,age=fa:la)
b<-b[order(b$sub_area,b$species.n,b$year,b$quarter,b$age),]
length(WSEA);length(PROPMAT);length(M);length(PROP_M2)
b<-data.frame(b,WSEA=WSEA, PROPMAT=PROPMAT,M=M,M1=M1,PROP_M2=PROP_M2)
b<-subset(b,quarter>=rec.season | age>fa,select=c(year,species.n,quarter,age,sub_area,WSEA,PROPMAT,M,M1,PROP_M2))

out<-c(out,list(bio=b))
################  other_sp

if (otherPredExist) {
  WSEA<-scanData('west.in')
  WSEA<-WSEA[1:((first.VPA-1)*noAreas*ny*(la-fa+1)*nq)]
  length(WSEA)
  N<-scanData('other_pred_n.in')
  b<-expand.grid(sub_area=1:noAreas,species.n=1:(first.VPA-1),year=years[1]:years[2],quarter=1:nq,age=fa:la)
  b<-b[order(b$sub_area,b$species.n,b$year,b$quarter,b$age),]
  length(N)
  if ( length(WSEA) !=length(N) ) stop('differen number of observations in west.in for "other predators") and other_pred_n.in')

  b<-data.frame(b,WSEA=WSEA, N=N)
  b<-subset(b,quarter>=rec.season | age>fa,select=c(year,species.n,species.n,age,quarter,sub_area,WSEA,N))
  b$N<-b$N*catchMultiplier

  out<-c(out,list(other=b))

  ###################### mean l

  l<-scanData('lsea.in')
  b<-expand.grid(SMS_area=1:noAreas,species.n=1:nsp,year=years[1]:years[2],quarter=1:nq,age=fa:la)
   b<-b[order(b$SMS_area,b$species.n,b$year,b$quarter,b$age),]

  b<-data.frame(b,mean_l=l)
  b<-subset(b,quarter>=rec.season | age>fa,select=c(year,species.n,age,quarter,SMS_area,mean_l))
  out<-c(out,list(mean_l=b))

  ###########################  consum

  CONSUM<-scanData('consum.in')
  b<-expand.grid(SMS_area=1:noAreas,species.n=1:npr,year=years[1]:years[2],quarter=1:nq,age=fa:la)
  b<-b[order(b$SMS_area,b$species.n,b$year,b$quarter,b$age),]

  b<-data.frame(b,CONSUM=CONSUM)
  b<-subset(b,quarter>=rec.season | age>fa,select=c(year,species.n,quarter,age,SMS_area,CONSUM))
  out<-c(out,list(consum=b))
}
 return(out)
}

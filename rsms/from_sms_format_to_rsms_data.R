From_SMS_format_to_rsms<-function(otherPredExist=TRUE,catchMultiplier=1,dir=data.path) {
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

CATCHN<-head(scan(file.path(dir,'canum.in'),comment.char='#'),-1)
WCATCH<-head(scan(file.path(dir,'weca.in'),comment.char='#'),-1)
Prop.landed<-head(scan(file.path(dir,'proportion_landed.in'),comment.char='#'),-1)
Prop.landed<-Prop.landed[1:length(WCATCH)]
b<-expand.grid(sub_area=1:noAreas,species.n=first.VPA:nsp,year=years[1]:years[2],quarter=1:nq,age=fa:la)

b<-b[order(b$sub_area,b$species.n,b$year,b$quarter,b$age),]
b<-data.frame(b,CATCHN=CATCHN,WCATCH=WCATCH,PROP_CAT=Prop.landed)
b<-subset(b,select=c(year,species.n,quarter,sub_area,age,WCATCH,CATCHN,PROP_CAT))
b$CATCHN<-catchMultiplier*b$CATCHN
out<-list(catch=b)
############## bio data

WSEA<-head(scan(file.path(dir,'west.in'),comment.char='#'),-1)
WSEA<-WSEA[((first.VPA-1)*noAreas*ny*(la-fa+1)*nq+1):length(WSEA)]
PROPMAT<-head(scan(file.path(dir,'propmat.in'),comment.char='#'),-1)
M<-head(scan(file.path(dir,'natmor.in'),comment.char='#'),-1)
M1<-head(scan(file.path(dir,'natmor1.in'),comment.char='#'),-1)
PROP_M2<-head(scan(file.path(dir,'n_proportion_m2.in'),comment.char='#'),-1)

b<-expand.grid(sub_area=1:noAreas,species.n=first.VPA:nsp,year=years[1]:(years[2]),quarter=1:nq,age=fa:la)
b<-b[order(b$sub_area,b$species.n,b$year,b$quarter,b$age),]

b<-data.frame(b,WSEA=WSEA, PROPMAT=PROPMAT,M=M,M1=M1,PROP_M2=PROP_M2)
b<-subset(b,quarter>=rec.season | age>fa,select=c(year,species.n,quarter,age,sub_area,WSEA,PROPMAT,M,M1,PROP_M2))

out<-c(out,list(bio=b))
################  other_sp

if (otherPredExist) {
  WSEA<-head(scan(file.path(dir,'west.in'),comment.char='#'),-1)
  WSEA<-WSEA[1:((first.VPA-1)*noAreas*ny*(la-fa+1)*nq)]
  length(WSEA)
  N<-head(scan(file.path(dir,'other_pred_n.in'),comment.char='#'),-1)
  b<-expand.grid(sub_area=1:noAreas,species.n=1:(first.VPA-1),year=years[1]:years[2],quarter=1:nq,age=fa:la)
  b<-b[order(b$sub_area,b$species.n,b$year,b$quarter,b$age),]
  length(N)
  if ( length(WSEA) !=length(N) ) stop('differen number of observations in west.in for "other predators") and other_pred_n.in')

  b<-data.frame(b,WSEA=WSEA, N=N)
  b<-subset(b,quarter>=rec.season | age>fa,select=c(year,species.n,species.n,age,quarter,sub_area,WSEA,N))
  b$N<-b$N*catchMultiplier

  out<-c(out,list(other=b))

  ###################### mean l

  l<-head(scan(file.path(dir,'lsea.in'),comment.char='#'),-1)
  b<-expand.grid(SMS_area=1:noAreas,species.n=1:nsp,year=years[1]:years[2],quarter=1:nq,age=fa:la)
   b<-b[order(b$SMS_area,b$species.n,b$year,b$quarter,b$age),]

  b<-data.frame(b,mean_l=l)
  b<-subset(b,quarter>=rec.season | age>fa,select=c(year,species.n,age,quarter,SMS_area,mean_l))
  out<-c(out,list(mean_l=b))

  ###########################  consum

  CONSUM<-head(scan(file.path(dir,'consum.in'),comment.char='#'),-1)
  b<-expand.grid(SMS_area=1:noAreas,species.n=1:npr,year=years[1]:years[2],quarter=1:nq,age=fa:la)
  b<-b[order(b$SMS_area,b$species.n,b$year,b$quarter,b$age),]

  b<-data.frame(b,CONSUM=CONSUM)
  b<-subset(b,quarter>=rec.season | age>fa,select=c(year,species.n,quarter,age,SMS_area,CONSUM))
  out<-c(out,list(consum=b))
}
 return(out)
}

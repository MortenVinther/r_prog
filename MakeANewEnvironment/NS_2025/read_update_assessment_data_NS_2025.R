#library(remotes)
#install_github("ices-tools-prod/icesTAF")
#library(icesTAF)
# ?icesTAF

old.dir<-getwd()

SMS<-SMS.control  # option with a shorter name

fa<-SMS@first.age


# old catches
OldCat<-read.csv(file=file.path(root,exchangeDir,paste0('VPA_Ca01_',old_key_label,'.csv')))
summary(OldCat)

SMS@species.info


NewCat<-OldCat

NewCat$sub_area<-NewCat$cat<-NULL

OldBio<-read.csv(file=file.path(root,exchangeDir,paste0('VPA_Bi01_',old_key_label,'.csv')))
#any(duplicated(select(OldBio, year,species,quarter, age,sub_area )))

# copy old bio to newbio as a start
NewBio<-subset(OldBio,year>=firstY)
b<-subset(NewBio,year==lastYold)
for (y in ((lastYold+1):lastY)) {
  b$year<-y
  NewBio<-rbind(NewBio,b)
}
NewBio$ICES_WSEA<- -1

names(NewBio)
#  any(duplicated(select(NewBio, year,species,quarter, age,sub_area )))


check_bio<-function(sp='COD',y=1992,a=3){
  subset(NewBio,species==sp & year==y  & age==a)
}

check_bio(sp="COD",y=2022,a=5)


if (add_cod_Mmagic) {
 # The 2021 cod benchmark introduced an ad hoc adjustment on natural mortality of ages 3+ from 2011 to mimic a 15% for migration to the West of Scotland area, which is not included in the assessment area.  M1 used for SMS is adjusted using the same approach:
   # 〖M1〗_(adj,a,y)=〖M1〗_(a,y)-ln⁡〖(1-α)〗
  criteria<-NewBio$species=='COD' & NewBio$year>=2011 & NewBio$age>=3
  NewBio[criteria,'M1']<-NewBio[criteria,'M1'] -log(1-0.15/4)
}


check_bio(sp="COD",y=2022,a=5)



source(file.path(root.prog,'r_prog','MakeANewEnvironment','read_csv_fish.R'))

arr2dfName <- function(arr,idx=NULL,name='y') {
  if(is.null(dimnames(arr))){dimnames(arr)<-lapply(dim(arr), FUN=function(x)1:x)}
  dn <- dimnames(arr)
  #if (any(unlist(lapply(dn,is.null)))) stop('Length of dimnames must equal length of dimension.')
  for (i in 1:length(dim(arr))) if (is.null(dn[[i]])) dn[[i]]<-as.character(1:(dim(arr)[i]))
  if(is.null(names(dn))){names(dn)<-paste('index', 1:length(dim(arr)), sep=".")}
  ans <- cbind(expand.grid(dn,stringsAsFactors = FALSE),as.vector(arr))
  if (is.null(idx)) idx<-c(colnames(ans)[-ncol(ans)])
  colnames(ans)<-c(idx,name)
  return(as.data.frame(ans))
}


doCatch<-function(fa=fa,la=la,cn=cn,cw=cw,ln=ln,species='AAA') {
  a<-merge( arr2dfName(cn,idx=c('year','age'),name='CATCHN'), arr2dfName(cw,idx=c('year','age'),name='WCATCH'))

  b<-arr2dfName(ln,idx=c('year','age'),name='lanN')
  a<-merge(a,b)

  a[as.numeric(a$age)>la,'age']<-la
  a$sop<-a$CATCHN*a$WCATCH
  a<-aggregate(cbind(CATCHN,sop,lanN)~year+age,sum,na.rm=T,data=a)
  a$WCATCH<-a$sop/a$CATCHN
  a[is.na(a$WCATCH),'WCATCH']<-0

  a$PROP_CAT<-a$lanN/a$CATCHN
  a[is.na(a$PROP_CAT),'PROP_CAT']<-1
  a[a$PROP_CAT>1,'PROP_CAT']<-1
  a$species<-species
  a$quarter<-1
  a$sop<-a$lanN<-NULL
  a$age<-as.numeric(a$age)
  a$year<-as.numeric(a$year)

  save(a,file=paste0(species,'_catch.Rdata'))
}


doBio01<-function(fa=fa,la=la,mo=mo,sw=sw,nm=nm,species='AAA') {
  cat('Species:',species,'\n')
  a<-inner_join(arr2dfName(mo,idx=c('year','age'),name="PROPMAT"), arr2dfName(nm,idx=c('year','age'),name='M'),by = join_by(year, age))
  #sort(unique(a$year))
  a<-subset(a,year<=as.character(lastY))
  a<-inner_join(a,arr2dfName(sw,idx=c('year','age'),name='ICES_WSEA'),by = join_by(year, age))
  a<-subset(a,as.numeric(a$age)<=la & year<=as.character(lastY))
  a$species<-species
  a$quarter<-9
  a$age<-as.numeric(a$age)
  a$year<-as.numeric(a$year)
  #print(subset(a,year==2022 & age==5))
  save(a,file=paste0(species,'_bio01.Rdata'))
}


updBIO<-function(sp_name=sp_name,sp=species,fdir=fdir,firstY=firstY,q=FALSE) {
  load(file=file.path(fdir,paste0(sp,'_bio01.Rdata')),verbose=F) #load a

  a<-subset(a,year>=firstY & year <=lastY,select=c(species,year,quarter,age,PROPMAT,M,ICES_WSEA))
  tst2<-select(a,year,age); tst2 %>% mutate(dup=duplicated(tst2)) %>%filter(dup)
  filter(a,year==2020)
  newYear<-sort(unique(a$year))
  if (!q) a$quarter<-NULL
  b<-subset(NewBio,(species==sp & year %in% newYear))

  if (q) a$newM<-a$M else a$newM<-a$M/4

  a$newPROPMAT<-a$PROPMAT;
  a$M<-a$PROPMAT<-NULL
  a$newICES_WSEA<-a$ICES_WSEA; a$ICES_WSEA<-NULL
  if (q) ab<-merge(x=a,y=b,by=c("species","year","quarter","age"), all.y=T)

  if (!q) ab<-merge(x=a,y=b,by=c("species","year","age"), all.y=T)
  ab$PROPMAT<-ifelse(is.na(ab$newPROPMAT),ab$PROPMAT,ab$newPROPMAT)
  ab$M<-ifelse(is.na(ab$newM),ab$M,ab$newM)
  ab$ICES_WSEA<-ifelse(is.na(ab$newICES_WSEA),ab$ICES_WSEA,ab$newICES_WSEA)
  ab$newPROPMAT<-ab$newM<-ab$newICES_WSEA<-NULL
  NewBio<-subset(NewBio,!(species==sp & year %in% newYear))

  NewBio<-rbind(ab,NewBio)
  dim(NewBio)
  return(NewBio)
}


updCATCH<-function(sp_name=sp_name,sp=species,fdir=fdir,firstY=firstY) {
  load(file=file.path(fdir,paste0(sp,'_catch.Rdata')),verbose=F) #load a
  a<-subset(a,year>=firstY & year<=lastY)
  newYear<-sort(unique(a$year))
  NewCat<-subset(NewCat,!(species==sp & year %in% newYear))
  NewCat<-rbind(a,NewCat)
  return(NewCat)
}

plot_bio<-function(oldD,oldLabel=paste(lastYold+1,'run',sep='-'),
                   newD,newLabel=paste(lastY+1,'run',sep='-'),
                   addD=NA,sp,years=firstY:lastY) {
  oldD$label<-oldLabel
  newD$label<-newLabel
  oldD<-subset(oldD,species==sp & year %in% years, select=c(year,species,quarter,age,WSEA, PROPMAT,M, M1,PROP_M2,label))
  newD<-subset(newD,species==sp & year %in% years, select=c(year,species,quarter,age,WSEA, PROPMAT,M, M1,PROP_M2,label))
  a<-rbind(oldD,newD)
  if (!is.na(addD)) a<-rbind(a,add)
  a<-as_tibble(a) %>% mutate(aq=paste0("Q",quarter," Age:", formatC(age,width=2)),year=as.double(year))

  p<-ggplot(a, aes(x=year, y=WSEA, group=label, color=label)) +
    geom_line()+
    geom_point()+
    scale_x_continuous(name="Year") +
    scale_y_continuous(name="Weight (Kg)")+

    facet_wrap(vars(aq))
  p
}




# cod, from TAF
# I have not changed the stock area and I use the (still updated) North Sea catches, extract new survey indices (with the old area) and invent stok mean weight etc.
species<-'COD'  # age 1 to 15, 1963-2022
sp_name<-'Cod'
fdir <-file.path(root,exchangeDir,'ByStock',species)
fdir
# I have just downloaded the taf manually!
#outf <-file.path(root,exchangeDir,'ByStock',species,'codTAF.zip')
#download.file("https://github.com/ices-taf/2023_cod.27.46a7d20_assessment",destfile=outf,mode="wb")

# install.packages("C:/_C_drev/SMS-git/Data_NorthSea/input_NS_2023/ByStock/COD/stockassessment_0.12.2.zip", repos = NULL, type = "win.binary") #from ICES WGNSSK Sharepoint
library(stockassessment)


# install.packages("C:/_C_drev/SMS-git/Data_NorthSea/input_NS_2023/ByStock/COD/multiStockassessment_0.3.3.zip", repos = NULL, type = "win.binary") #from ICES WGNSSK Sharepoint
#library(multiStockassessment)


codTAF<-file.path(fdir,"codTAF")
codTAF
setwd(codTAF)


#extract new indices from the same GAM models as used by WGNSSK (Takes a night to run)
if (FALSE) {
  source(file.path(fdir,"data_prepareSurveyQ1.R"))
  source(file.path(fdir,"data_prepareSurveyQ1.R"))
  source(file.path(fdir,"data_indicesQ1.R"))
  source(file.path(fdir,"data_indicesQ34.R"))
}
# produces
# q34index0_7_Old_area_MV.txt, old area indices
# NorthwestQ1MV.txt # quarter 1, but only the North Sea part of NorthWest
# SouthQ1MV.txt
# VikingQ1MV.txt

setwd(old.dir)


source(file.path(fdir,"SMS-data","datascript.R"))
la<-SMS@species.info[sp_name,'last-age']
cat(species,' last age:',la,'\n')
fdir

dat<-datascript(cutage=la,low=1,fdir=file.path(root,exchangeDir,'ByStock',species,'SMS-data'))

setwd(fdir)
doCatch(fa=fa,la=la,cn=dat$cn,cw=dat$cw,ln=dat$ln,species=species)


if (FALSE) {
  #unallocated catches
  load(file=paste0(species,'_catch.Rdata'),verbose=T)

  Cscale<-data.frame(species='COD',year=1993:2005, scale=c(0.94, 1.04, 1.18, 1.03, 0.87, 0.74, 0.89, 0.93, 1.27, 0.87, 1.50, 1.12, 0.96))
  a<-merge(x=a,y=Cscale,all.x=TRUE)
  a$CATCHN<-ifelse(is.na(a$scale),a$CATCHN,a$CATCHN*a$scale)
  a$scale<-NULL
  save(a,file=paste0(species,'_catch.Rdata'))
  #apply(dat$cn*dat$cw,1,sum)
}

doBio01(fa=fa,la=la,mo=dat$mo,sw=head(dat$sw,-1),nm=dat$nm,species=species)

all.surveys<-list(COD=dat$surveys)
#update catch  and bio data
#NewCat<-subset(NewCat,!(species==species & year %in% newYear))
#NewCat<-rbind(a,NewCat)

NewCat<-updCATCH(sp_name=sp_name,sp=species,fdir=fdir,firstY=firstY)

#tst<-subset(NewCat,species==species); tapply(tst$CATCHN*tst$WCATCH,list(tst$year),sum)
check_bio(sp=species,y=2022,a=5)

NewBio<-updBIO(q=F,sp_name=sp_name,sp=species,fdir=fdir,firstY=firstY)
check_bio(sp=species,y=2022,a=5)

subset(NewBio,species=='COD' & year==2022 & quarter==2)
#


#plot_bio(oldD=OldBio,oldLabel=paste(lastYold+1,'run',sep='-'),
#         newD=NewBio,newLabel=paste(lastY+1,'run',sep='-'),
#         addD=NA,sp=sp,years=2015:2019)


############################


# Whiting, Data from stockassessment.org
species<-'WHG' # 1978	2022	ages 0	6+
sp_name<-'Whiting'
fdir <-file.path(root,exchangeDir,'ByStock',species,'SMS-data')
# single species assessment uses 6+, SMS used 8+.
source(file.path(fdir,"datascript.R"))

la<-SMS@species.info[sp_name,'last-age']
la<-6
cat(species,' last age:',la,'\n')
dat<-datascript(fdir=fdir)
str(dat,1)
if (TRUE) {
  tmp<-dat$sw
  # raise the very low mean weights
  lowW<-apply(dat$sw,2,function(x) quantile(x,probs=0.10))
  for (i in (1:dim(dat$sw)[[2]])) dat$sw[dat$sw[,i]<lowW[i],i ] <-lowW[i]
  dat$sw
  dat$sw-tmp
}

setwd(fdir)
doCatch(fa=fa,la=la,cn=dat$cn,cw=dat$cw,ln=dat$ln,species=species)
#load(file=paste0(species,'_catch.Rdata'),verbose=T)

doBio01(fa=fa,la=la,mo=dat$mo,sw=dat$sw,nm=dat$nm,species=species)
all.surveys[[species]]<-dat$surveys
names(all.surveys)

#update catch  and bio data
NewCat<-updCATCH(sp_name=sp_name,sp=species,fdir=fdir,firstY=firstY)
#subset(NewCat,species=='WHG' & year==1974 & quarter==1)
check_bio(sp=species)
NewBio<-updBIO(q=F,sp_name=sp_name,sp=species,fdir=fdir,firstY=firstY)
check_bio(sp=species)



# Haddock, Data from stockassessment.org
species<-'HAD'  #1965	2022	ages 0	15
sp_name<-'Haddock'
fdir <-file.path(root,exchangeDir,'ByStock',species,'SMS-data')
# single species assessment uses 8+, SMS uses 10+.
source(file.path(fdir,"datascript.R"))
la<-SMS@species.info[sp_name,'last-age']
cat(species,' last age:',la,'\n')

dat<-datascript(cutage=la,low=0,fdir=fdir)

setwd(fdir)
doCatch(fa=fa,la=la,cn=dat$cn,cw=dat$cw,ln=dat$ln,species=species)
doBio01(fa=fa,la=la,mo=dat$mo,sw=dat$sw,nm=dat$nm,species=species)
all.surveys[[species]]<-dat$surveys
names(all.surveys)

#update catch  and bio data
NewCat<-updCATCH(sp_name=sp_name,sp=species,fdir=fdir,firstY=firstY)
check_bio(sp=species)
NewBio<-updBIO(q=F,sp_name=sp_name,sp=species,fdir=fdir,firstY=firstY)
check_bio(sp=species)


#Saithe, from TAF (updated in 2025)
species<-'POK'  # 1967	2022, ages 3	10
sp_name<-'Saithe'

fdir <-file.path(root,exchangeDir,'ByStock',species,'SMS-data')

source(file.path(fdir,"datascript.R"))
la<-SMS@species.info[sp_name,'last-age']
cat(species,' last age:',la,'\n')

dat<-datascript(cutage=la,fdir=fdir)

setwd(fdir)
doCatch(fa=fa,la=la,cn=dat$cn,cw=dat$cw,ln=dat$ln,species=species)
doBio01(fa=fa,la=la,mo=dat$mo,sw=dat$sw,nm=dat$nm,species=species)
all.surveys[[species]]<-dat$surveys
names(all.surveys)

#update catch  and bio data
NewCat<-updCATCH(sp_name=sp_name,sp=species,fdir=fdir,firstY=firstY)
check_bio(sp=species)
NewBio<-updBIO(q=F,sp_name=sp_name,sp=species,fdir=fdir,firstY=firstY)
check_bio(sp=species)

#tst<-subset(NewCat,NewCat$species=='POK'); unique(tst$species); tapply(tst$CATCHN*tst$WCATCH,list(tst$year),sum)
#tst<-subset(NewCat,NewCat$species=='POK'); tapply(tst$CATCHN,list(tst$year,tst$age),sum)




# Mackerel from SharePoint- 
species<-'MAC'  # 1980	2022, ages 0	12
sp_name<-'Mackerel'

fdir <-file.path(root,exchangeDir,'ByStock',species,'SMS-data')
# single species assessment uses 12+, SMS uses 10+.
source(file.path(fdir,"datascript.R"))
la<-SMS@species.info['Mackerel','last-age']
cat(species,' last age:',la,'\n')
dat<-datascript(cutage=la,fdir=fdir)

setwd(fdir)
doCatch(fa=fa,la=la,cn=dat$cn,cw=dat$cw,ln=dat$ln,species=species)
doBio01(fa=fa,la=la,mo=dat$mo,sw=dat$sw,nm=dat$nm,species=species)
all.surveys[[species]]<-dat$surveys
names(all.surveys)

#update catch  and bio data
NewCat<-updCATCH(sp_name=sp_name,sp=species,fdir=fdir,firstY=firstY)
check_bio(sp=species)
NewBio<-updBIO(q=F,sp_name=sp_name,sp=species,fdir=fdir,firstY=firstY)
check_bio(sp=species)


load(file.path(fdir,"..","model fit.RData"),verbose=TRUE)


ntab <- ntable(fit.new)
ntab[,'10']<- rowSums(ntab[,as.character(10:12)])
ntab<-ntab[,as.character(0:10)]
#ftab <- faytable(fit.new)
###


# Herring from ices Sharepoint
species<-'HER' # 1947	2022	ages 0	8
sp_name<-'Herring'

fdir <-file.path(root,exchangeDir,'ByStock',species)
# install.packages("FLCore", repos="http://flr-project.org/R")
#install.packages("FLSAM", repos="http://flr-project.org/R")
#library(FLSAM)  # we dont need it

if (FALSE) { # have to be run manually (error message), as FLSAM cannot be loaded
  load(file.path(fdir,"NSAS_HAWG2023_sf.Rdata"),verbose=TRUE)
  # Loading objects:  NSH NSH.tun NSH.ctrl NSH.sam
  outf<-file.path(root,exchangeDir,'ByStock',species,'SMS-data','her')
writeFLStock(NSH,output.file=outf)
 save(NSH.tun,file=file.path(fdir,"NSAS_HAWG2023_sf_tun.Rdata"))
} else {
  load(file=file.path(fdir,"NSAS_HAWG2023_sf_tun.Rdata"),verbose=TRUE)
}

outf<-file.path(root,exchangeDir,'ByStock',species,'SMS-data','index.dat')
cat('NS_herring\n',file=outf)

cat(length(NSH.tun)+100-4,'\n',append=TRUE,file=outf)

inames<-NSH.tun@names
i<-1
for (idc in NSH.tun) if (i<=4) {  # exclude LAI survey (larvae index, but used as SSB index, which SMS cannot use)
  # i<-8; idc<-NSH.tun[[i]]
  cat(inames[i],'\n',append=TRUE,file=outf)
  cat(inames[i],'\n')
  ran<-idc@range
  cat(ran['minyear'],ran['maxyear'],'\n',append=TRUE,file=outf)
  cat(1,1,ran['startf'],ran['endf'],'\n',append=TRUE,file=outf)
  cat(ran['min'],ran['max'],'\n',append=TRUE,file=outf)
  if (dim(idc)[1]>1) idx<-t(idc@index@.Data[,,,,,]) else idx<-idc@index@.Data[,,,,,]
  idx<-cbind(idc@effort@.Data[,,,,,],idx)
  write.table(idx,row.names=FALSE,col.names=FALSE,append=TRUE,file=outf)
  i<-i+1
}


fdir <-file.path(root,exchangeDir,'ByStock',species,'SMS-data')

source(file.path(fdir,"datascript.R"))
la<-SMS@species.info['Herring','last-age']
cat(species,' last age:',la,'\n')
dat<-datascript(cutage=la,fdir=fdir)
dat$cn[is.na(dat$cn)]<-0

setwd(fdir)
doCatch(fa=fa,la=la,cn=dat$cn,cw=dat$cw,ln=dat$ln,species=species)
doBio01(fa=fa,la=la,mo=dat$mo,sw=dat$sw,nm=dat$nm,species=species)
all.surveys[[species]]<-dat$surveys
names(all.surveys)

# update with quarterly data

# new quarterly catch data
a<-read.csv(file.path(root,exchangeDir,'ByStock',species,'canum_weca_2005_2022.csv'))
head(a);dim(a)
b<-reshape(a,direction='long',varying = list(4:21),v.names='value')
head(b)
b$year<-2022-b$time+1
b1<-subset(b,type=='CANUM',select=c(year,quarter,age,value))
b1$CATCHN<-b1$value; b1$value<-NULL
b2<-subset(b,type=='WECA',select=c(year,quarter,age,value))
b2$WCATCH<-b2$value;  b2$value<-NULL

b<-merge(b1,b2)
b<-subset(b,!(age==0 & quarter %in% c(1,2)))
b$CATCHN<-b$CATCHN*1000

b[b$age>la,'age']<-la
b<-aggregate(cbind(WCATCH,CATCHN,cw=CATCHN*WCATCH)~year+quarter+age,data=b,sum)
b$WCATCH<-b$cw/b$CATCHN
b$cw<-NULL



#round(ftable(tapply(b$CATCHN,list(b$year,b$quarter,b$age),sum))  )
b$species<-species
b$PROP_CAT<-1
b

qy<-unique(b$year) # years with quarterly data

load(file=paste0(species,'_catch.Rdata'),verbose=T)
a<-subset(a,!(year %in% qy) & quarter==1 & year>=firstY,select=c(year,age,CATCHN))
a$sumC<-a$CATCHN; a$CATCHN<-NULL;

a2<-subset(NewCat,!(year %in% qy) & species=='HER')  #Old quarterly data
a3<-aggregate(cbind(sumCq=CATCHN)~year+age,data=a2,sum)

a<-merge(a,a3)
a$fac<-a$sumC/a$sumCq
a$sumC<-a$sumCq<-NULL

a<-merge(x=a2,y=a,x.all=TRUE)
a$CATCHN<-a$CATCHN*a$fac; a$fac<-NULL
a<-rbind(b,a)


round(ftable(tapply(a$CATCHN,list(a$year,a$quarter,a$age),sum))  )

save(a,file=paste0(species,'_catch.Rdata'))

#update catch  and bio data
tail(subset(NewCat,species=='HER'& year==2017 & quarter==1),10)
NewCat<-updCATCH(sp_name=sp_name,sp=species,fdir=fdir,firstY=firstY)
tail(subset(NewCat,species=='HER'& year==2017 & quarter==1),10)
summary(filter(NewCat,species=='HER'))


check_bio(sp=species,y=2016)
NewBio<-updBIO(q=FALSE,sp_name=sp_name,sp=species,fdir=fdir,firstY=firstY)
check_bio(sp=species,y=2016,a=9)

delher<-NewBio$species=='HER' & NewBio$age>8
NewBio[delher,'WSEA']<- -1
NewBio[delher,'ICES_WSEA']<- -1
NewBio[delher,'PROP_M2']<- -1
NewBio[delher,'M']<- -1
NewBio[delher,'M1']<--1
NewBio[delher,'PROPMAT']<--1

if (TRUE) {
  ###### update Sandeel stocks
  species<-'sandeel'
  a<-read.csv(file.path(root,exchangeDir,'ByStock',species,'Total_catch_in_numbers_and_mean_weight_length_2022.csv')) %>% filter(aar!= 1982)

  a0<-a; a0$age<-0; a0$WCATCH<-a$mw0; a0$CATCHN<-a$n0; a0$MEANL<-a$ml0;
  a1<-a; a1$age<-1; a1$WCATCH<-a$mw1; a1$CATCHN<-a$n1; a1$MEANL<-a$ml1;
  a2<-a; a2$age<-2; a2$WCATCH<-a$mw2; a2$CATCHN<-a$n2; a2$MEANL<-a$ml2;
  a3<-a; a3$age<-3; a3$WCATCH<-a$mw3; a3$CATCHN<-a$n3; a3$MEANL<-a$ml3;
  a4<-a; a4$age<-4; a4$WCATCH<-a$mw4; a4$CATCHN<-a$n4; a4$MEANL<-a$ml4;
  a<-rbind(a0,a1,a2,a3,a4)
  a<-bind_rows(a0,a1,a2,a3,a4) %>% as_tibble()
  a
  sort(unique(a$area))

  a<- a %>%mutate(species=if_else (area  %in% c('N','SH'),'NSA','SSA')) %>% select(species,aar,hy,age,n_samples,WCATCH,CATCHN,MEANL)
  summary(a)
  a <- a %>%
    group_by(species,aar,hy,age) %>%
    summarize(WCATCH=weighted.mean(x=WCATCH, w=CATCHN+0.01,na.rm = FALSE ) ,
              MEANL=weighted.mean(x=MEANL, w=CATCHN+0.01,na.rm = FALSE ),
              CATCHN=sum(CATCHN),n_samples=sum(n_samples)) %>% ungroup()

  a %>% filter(species=='SSA') %>% arrange(desc(WCATCH))
  a %>% filter(species=='NSA') %>% arrange(desc(WCATCH))

  tst<-a%>%filter(species=='NSA') %>%transmute(species, aar,hy,age,WCATCH,CATCHN)
  round(ftable(tapply(tst$CATCHN,list(tst$aar,tst$hy,tst$age),sum)),0)


  any(a %>% select(species,aar,hy,age) %>% duplicated)


  san_length<-data.frame(year=a$aar,quarter=a$hy+1,species=a$species,age=a$age,n_samples=a$n_samples,MEANL=a$MEANL)
  san_length[ is.na(san_length$n_samples),'n_samples']<-2
  head(san_length)
  summary(san_length)

  a<-data.frame(year=a$aar,species=a$species,age=a$age,PROP_CAT=1,quarter=a$hy+1,CATCHN=a$CATCHN,WCATCH=a$WCATCH,PROP_CAT=1)
  tst<- a %>% select(species,year,quarter,age)
  any(duplicated(tst))

  aa<-subset(a,!(age==0 & quarter %in% c(1,2)),select=c( species, quarter, year, age, CATCHN, WCATCH, PROP_CAT))

  sp<-'NSA'
  cat(sp,'\n')
  a<-subset(aa,species==sp)
  newYear<-sort(unique(a$year))
  NewCat<-subset(NewCat,!(species==sp & year %in% newYear))
  NewCat<-rbind(a,NewCat)


  id<-read.csv(file.path(root,exchangeDir,'ByStock',species,'North_2023.csv'))
  y<-id$Year
  id2<-as.matrix(id[2:3])
  rownames(id2)<-y
  colnames(id2)<-0:1
  attr(id2,"time")<-c(0.7,0.9)
  id2
  all.surveys[[sp]]<-list(NSA_drgde=id2)
  names(all.surveys[[sp]])

  sp<-'SSA'
  cat(sp,'\n')
  a<-subset(aa,species==sp)
  newYear<-sort(unique(a$year))
  NewCat<-subset(NewCat,!(species==sp & year %in% newYear))
  NewCat<-rbind(a,NewCat)

  NewCat[NewCat$species %in% c('NSA','SSA'),'PROP_CAT']<-1

  id<-read.csv(file.path(root,exchangeDir,'ByStock',species,'South_2023.csv'))
  y<-id$Year
  id2<-as.matrix(id[2:3])
  rownames(id2)<-y
  colnames(id2)<-0:1
  attr(id2,"time")<-c(0.7,0.9)
  all.surveys[[sp]]<- list(SSA_drgde=id2)


  # bio data both sandeels

  a<-subset(aa, WCATCH >0) %>%  select(species, quarter,year, age,WCATCH)
  summary(a$year)
  if (FALSE) {
   #old data to get the period before 1983
    old<-read.csv("C:/_C_drev/SMS-git/Data_NorthSea/SMS-input-2014-key-run-2014/VPA_Ca01.IN") %>% filter(species %in% c('NSA','SSA') & year<=1982) %>%
         select(species, quarter,year, age,WCATCH)

    a<-rbind(old,a)
  }
  old<-NewCat %>% filter(species %in% c('NSA','SSA') & year<=1982 & age<=4) %>%
       select(species, quarter,year, age,WCATCH)
  a<-rbind(old,a)

  tst<-a%>% filter(quarter %in% c(2,3) & age>0) %>% mutate(q=paste0(species,"_Q",quarter),Age=paste('Age',age))
  summary(tst)
  sort(unique(tst$q))

 png(filename=file.path(root,exchangeDir,'ByStock',species,"sandeel_mw_catch.png"),width=700, height = 700);
  #par(mar=c(3,5,3,2))
  ggplot(tst,aes(x=year,y=WCATCH))+geom_point()+geom_line()+geom_smooth(color='red')+facet_grid(cols=vars(Age),row=vars(q),scales="free_y")+
   ylab("Mean weight in the catch (kg)")
  cleanup()

 # pp<-by(tst,list(tst$q,tst$age),function(x) data.frame(x,newW=predict(stats::loess(WCATCH~year,data=subset(x,WCATCH>0)))))
  pp<-by(tst,list(tst$q,tst$age),function(x) data.frame(subset(x,WCATCH>0),newW=predict(stats::loess(WCATCH~year,data=subset(x,WCATCH>0)))))
  

  a1p<-do.call(rbind,pp) %>% select(species, quarter,year, age,newW)
  q1<-filter(a1p,quarter ==2 ) %>% mutate(quarter=1,CATCHN=NULL,newW=newW*0.79)
  q4<-filter(a1p,quarter ==3 ) %>% mutate(quarter=4,CATCHN=NULL,newW=newW*0.89)
  q1234<-rbind(q1,q4,a1p)
  summary(q1234)
  ffy<-filter(q1234,year==min(q1234$year)) %>% mutate(year=NULL)
  #misy<-data.frame(year=c(firstY:(min(q1234$year-1))))
  #q1234<-rbind(q1234,cross_join(misy,ffy))

  a0<-aa %>% filter(quarter==3 & age==0) %>% group_by(species) %>% summarize(newW=mean(WCATCH)) %>% ungroup() %>% mutate(quarter=4)
  a03<-a0 %>% mutate(quarter=3,newW=newW*0.5)
  a0<-rbind(a0,a03) %>% mutate(age=0)
  #a0<-cross_join(a0,expand.grid(year=firstY:lastY))
  a0<-cross_join(a0,expand.grid(year=(min(q1234$year)+1):lastY))
  san<-rbind(a0,q1234)

  png(filename=file.path(root,exchangeDir,'ByStock',species,"sandeel_NSA_mw_stock.png"),width=700, height = 700)
  a<-filter(san,species=='NSA') %>% mutate(q=paste0("Q ",quarter),Age=paste('Age',age))
  ggplot(a,aes(x=year,y=newW))+geom_point()+geom_line()+facet_grid(rows=vars(Age),cols=vars(q),scales="fixed")+
    ylab("Mean weight in the sea (kg)")
  cleanup()

  png(filename=file.path(root,exchangeDir,'ByStock',species,"sandeel_SSA_mw_stock.png"),width=700, height = 700)
  a<-filter(san,species=='SSA') %>% mutate(q=paste0("Q ",quarter),Age=paste('Age',age))
  ggplot(a,aes(x=year,y=newW))+geom_point()+geom_line()+facet_grid(rows=vars(Age),cols=vars(q),scales="fixed")+
    ylab("Mean weight in the sea (kg)")
  cleanup()

  b<-subset(san,select=c(species,age,year,quarter,newW))
  check_bio(sp=sp,y=2002,a=3)

  NewBio<-left_join(x=NewBio,y=b,by = join_by(species, year, age, quarter))

  check_bio(sp=sp,y=2002,a=3)
  NewBio[!is.na(NewBio$newW),'WSEA']<-NewBio[!is.na(NewBio$newW) ,'newW']
  NewBio$newW<-NULL
  check_bio(sp=sp,y=2002,a=3)


  png(filename=file.path(root,exchangeDir,'ByStock',species,"sandeel_NSA_mw_stock_full.png"),width=700, height = 700)
  a<-filter(NewBio,species=='NSA' & age<=4 & WSEA>0) %>% mutate(q=paste0("Q ",quarter),Age=paste('Age',age))
  ggplot(a,aes(x=year,y=WSEA))+geom_point()+geom_line()+facet_grid(rows=vars(Age),cols=vars(q),scales="fixed")+
    ylab("Mean weight in the sea (kg)")
  cleanup()


}


#Sprat from SMS assessment
species<-'SPR'
fdir <-file.path(root,exchangeDir,'ByStock',species,'SMS-data')
sp_name<-'Sprat'

la<-SMS@species.info[sp_name,'last-age']
cat(species,' last age:',la,'\n')

setwd(fdir)
a<-read.table(file.path(fdir,"summary.out"),header=TRUE)

a<-data.frame(year=a$Year,species=species,age=a$Age,PROP_CAT=1,quarter=a$Quarter,CATCHN=a$C.obs,WCATCH=a$weca,M=a$M,ICES_WSEA=a$west,PROPMAT=a$propmat)
s1<-round(tapply(a$CATCHN*a$WCATCH,list(a$year,a$quarter),sum))
qShift<-function(a) {
  q12<-a$quarter %in% c(1,2)
  q34<-a$quarter %in% c(3,4)
  a[q12,'quarter']<-a[q12,'quarter']+2
  a[q34,'quarter']<-a[q34,'quarter']-2
  a[q34,'year']<-a[q34,'year']+1
  a[q34,'age']<-a[q34,'age']+1
  return(a)
}
a<-qShift(a)

s2<-round(tapply(a$CATCHN*a$WCATCH,list(a$year,a$quarter),sum))
cbind(s1,s2)
aa<-a
a<-subset(aa,year<=lastY ,select=c( year, species, age, PROP_CAT, quarter,WCATCH,CATCHN))

# plus group
a[a$age==4,'age']<-3
#a<-subset(a,age<=3)

filter(a,WCATCH<0)
a2<-aggregate(cbind(WCATCH,CATCHN,WCATCH*CATCHN)~year+ species+ age+ PROP_CAT+ quarter,data=a,sum)
a2$WCATCH2<-a2$V3/a2$CATCHN
filter(a2,WCATCH2<0)

subset(a2,age==3 &year==1975)
a[a$CATCHN<100,'CATCHN']<-0
a[a$CATCHN<100,'WCATCH']<-0

# SPR year: 1974 Q1 & Q2 from old SMS
ad<-data.frame(year=1974,species='SPR',age=rep(0:3,2),PROP_CAT=1,quarter=rep(c(1,2),each=4),
  CATCHN=c(0.00, 6272325.12,790695.82,85688.62 ,0.00 ,  2222931.24,    155564.40 ,    13872.68),
  WCATCH=c(   0.0000, 0.0060,  0.0129,0.0428 ,0.0000 , 0.0055, 0.0113, 0.0390))

a<-rbind(a,ad)
filter(a,WCATCH<0)
summary(a)
dim(a)
a2<-aggregate(cbind(WCATCH,CATCHN,WCATCH*CATCHN)~year+ species+ age+ PROP_CAT+ quarter,data=a,sum)
a2$WCATCH2<-a2$V3/a2$CATCHN
dim(a2)
filter(a2,WCATCH2<0)
filter(a2,is.na(WCATCH2))
a2[is.na(a2$WCATCH2),'WCATCH2']<-0
a2$WCATCH<-a2$WCATCH2; a2$WCATCH2<-NULL;a2$V3<-NULL

a<-a2

save(a,file=paste0(species,'_catch.Rdata'))

a<-subset(aa,year<=lastY & age <=la,select=c( year, quarter, species, age, M, ICES_WSEA,PROPMAT))


#update catch  and bio data
NewCat<-updCATCH(sp_name=sp_name,sp=species,fdir=fdir,firstY=firstY)

if (FALSE) {
 tst<-filter(NewCat,species=='SPR')
 filter(tst,WCATCH<=0 & CATCHN>0)
 filter(tst,WCATCH<=0)
}


# special proportion mature
ap<-subset(a,quarter==3,select=c(year, species, age,PROPMAT))
ap$newPROPMAT<-ap$PROPMAT
ap$PROPMAT<-NULL
a<-merge(x=a,y=ap,all.x=TRUE)
head(subset(a,quarter==1))
a$PROPMAT<-a$newPROPMAT
a$newPROPMAT<-NULL
save(a,file=paste0(species,'_bio01.Rdata'))
#

check_bio(sp=species)
NewBio<-updBIO(q=TRUE,sp_name=sp_name,sp=species,fdir=fdir,firstY=firstY)
check_bio(sp=species)

#Sprat survey
a<-read.table(file.path(fdir,'catch_survey_residuals.out'),header=TRUE)
a<-subset(a,fleet>0)
head(a)
a<-data.frame(year=a$Year,quarter=a$Quarter,species=species,age=a$Age,fleet=a$fleet, obs=a$observed/1000)
unique(paste(a$fleet,'Q:',a$quarter))
a<-qShift(a)

head(a)
#data.frame(fleet=1:3,fname=c("IBTS_Q1","IBTS_Q3","Acoustic"))
fl<-data.frame(fleet=1:4,fname=c("IBTS_Q1_Rec","IBTS_Q1","IBTS_Q3","Acoustic"))
fl
#fl<-subset(fl,fleet !=1)
#fl
a<-right_join(a,fl) %>% mutate(fleet=fleet-1)
summary(a)
unique(paste(a$fleet,'Q:',a$quarter))
ftable(round(tapply(a$obs,list(a$fleet,a$year,a$age),sum)))

fl2<-by(a,a$fname,function(x) {
  y<-tapply(x$obs,list(x$year,x$age),sum)
  attr(y,"time")<-c(0,1)
  return(y)
})

all.surveys[[species]]<-fl2




#Norway pout  from seasonal SAM assessment
species<-'NOP'
sp_name<-'Norway pout'
fdir <-file.path(root,exchangeDir,'ByStock',species,'data')
b1<-read.table(file.path(fdir,"aux_data.dat"),header=TRUE)
b1<-subset(b1,age<=3)
b2<-read.table(file.path(fdir,"obs.dat"),header=TRUE)
b2$age<-b2$ageFrom; b2$ageFrom<-b2$ageTo<-NULL
fleet<-subset(b2,fleet>1)
catch<-subset(b2,fleet==1)

catch$cn<-catch$obs; catch$obs<-catch$fleet<-NULL
b<-merge(b1,catch)
b$quarter<- (b$t1 %% 1)*4 +1
b$year<-(b$t1 %/% 1)

head(b)

fdir <-file.path(root,exchangeDir,'ByStock',species,'SMS-data')
setwd(fdir)
a<-data.frame(year=b$year,age=b$age,CATCHN=b$cn*1000,WCATCH=b$CW/1000,PROP_CAT=1,species=species,quarter=b$quarter)
#head(a)
#filter(a,CATCHN>0 & WCATCH==0)
a[a$CATCHN>0 & a$WCATCH==0,'WCATCH']<-0.050

save(a,file=paste0(species,'_catch.Rdata'))

a<-data.frame(year=b$year,age=b$age,M=b$M,ICES_WSEA=b$SW/1000,PROPMAT=b$PM,species=species,quarter=b$quarter)
save(a,file=paste0(species,'_bio01.Rdata'))

fleet$quarter<- (fleet$t1 %% 1)*4 +1
fleet$year<-(fleet$t1 %/% 1)
tapply(fleet$obs,list(fleet$fleet,fleet$quarter),sum,na.rm=TRUE)
fleet2<-subset(fleet, fleet==3 & quarter==1 & age %in% c(1,2,3) & year >=1984) ; fleet2$surv<-'IBTSQ1'
fleet3<-subset(fleet, fleet==4 & quarter==3 & age %in% c(0,1) & year>=1992)   ; fleet3$surv<-'EGFSQ3'
fleet4<-subset(fleet, fleet==5 & quarter==3 & age %in% c(0,1) & year >=1998)   ; fleet4$surv<-'SGFSQ3'
fleet5<-subset(fleet, fleet==6 & quarter==3 & age %in% c(2,3) & year >=1991)   ; fleet5$surv<-'IBTSQ3'

fl<-rbind(rbind(fleet2,fleet3),rbind(fleet4,fleet5))

fl2<-by(fl,fl$surv,function(x) {
  y<-tapply(x$obs,list(x$year,x$age),sum)
  attr(y,"time")<-c(0,1)
  return(y)
})

all.surveys[[species]]<-fl2




#update catch  and bio data
NewCat<-updCATCH(sp_name=sp_name,sp=species,fdir=fdir,firstY=firstY)
subset(NewBio,species=='NOP' & year==1992 & quarter==1 & age==3)
NewBio<-updBIO(q=T,sp_name=sp_name,sp=species,fdir=fdir,firstY=firstY)
subset(NewBio,species=='NOP' & year==1992 & quarter==1 & age==3)


# Plaice, from TAF
species<-'PLE'
outf<-file.path(root,exchangeDir,'ByStock',species,'SMS-data','ple')

fdir<-file.path(root,exchangeDir,'ByStock','PLE','SMS-data')

source(file.path(fdir,"datascript.R"))
la<-SMS@species.info['Plaice','last-age']
dat<-datascript(cutage=la,fdir=fdir)

setwd(fdir)
doCatch(fa=fa,la=la,cn=dat$cn,cw=dat$cw,ln=dat$ln,species=species)
doBio01(fa=fa,la=la,mo=dat$mo,sw=dat$sw,nm=dat$nm,species=species)
all.surveys[[species]]<-dat$surveys
names(all.surveys)

#update catch  and bio data
NewCat<-updCATCH(sp_name=sp_name,sp=species,fdir=fdir,firstY=firstY)
NewBio<-updBIO(q=F,sp_name=sp_name,sp=species,fdir=fdir,firstY=firstY)
#subset(NewBio,species=='PLE' & year==2010 & quarter==2)

#change M1 to M for plaice
NewBio<-NewBio %>% mutate(M1=if_else(species=='PLE',M,M1))
pleM1age1<-mean(subset(NewBio,species=='PLE' & age==1 )$M)
NewBio[NewBio$species=='PLE' & NewBio$age==0 & NewBio$quarter %in% c(3,4),'M1']<-pleM1age1



# Sole,from WGNSSK Sharepoint
species<-'SOL'
outf<-file.path(root,exchangeDir,'ByStock',species,'SMS-data','SOL')
if (FALSE) {
  load(file.path(root,exchangeDir,'ByStock',species,'data.Rda'),verbose=TRUE)

  writeFLStock(stock,output.file=outf)


  outf<-file.path(root,exchangeDir,'ByStock',species,'SMS-data','index.dat')
  cat('Sole IV\n',file=outf)

  cat(length(indices)+100,'\n',append=TRUE,file=outf)

  inames<-indices@names
  i<-1
  for (idc in indices) {
    #idc<-ass.indices[[1]]
    cat(inames[i],'\n',append=TRUE,file=outf)

    ran<-idc@range
    cat(ran['minyear'],ran['maxyear'],'\n',append=TRUE,file=outf)
    cat(1,1,ran['startf'],ran['endf'],'\n',append=TRUE,file=outf)
    cat(ran['min'],ran['max'],'\n',append=TRUE,file=outf)
    idc@effort@.Data[,,,,,]<-1  # seems missing
    idx<-cbind(idc@effort@.Data[,,,,,],t(idc@index@.Data[,,,,,]))
    write.table(idx,row.names=FALSE,col.names=FALSE,append=TRUE,file=outf)
    i<-i+1
  }
}
fdir<-file.path(root,exchangeDir,'ByStock',species,'SMS-data')

source(file.path(fdir,"datascript.R"))
la<-SMS@species.info['Sole','last-age']
dat<-datascript(cutage=la,fdir=fdir)

setwd(fdir)
doCatch(fa=fa,la=la,cn=dat$cn,cw=dat$cw,ln=dat$ln,species=species)
doBio01(fa=fa,la=la,mo=dat$mo,sw=dat$sw,nm=dat$nm,species=species)
#dat$surveys[[3]]<-dat$surveys[[3]]*1000
all.surveys[[species]]<-dat$surveys
names(all.surveys)

#update catch  and bio data
NewCat<-updCATCH(sp_name=sp_name,sp=species,fdir=fdir,firstY=firstY)
NewBio<-updBIO(q=F,sp_name=sp_name,sp=species,fdir=fdir,firstY=firstY)

####################


#sink(file.path(finalExchangeDir,'surveys.dat'))
# lapply(all.surveys,print)
#sink()


#check

# ftable(round(tapply(NewCat$PROP_CAT,list(NewCat$species,NewCat$year,formatC(NewCat$age,digits=3,width=2,flag="0")),mean),2))
#ftable(round(tapply(NewCat$PROP_CAT,list(NewCat$species,NewCat$quarter,NewCat$year,formatC(NewCat$age,digits=0,width=2,flag="0")),sum)))

round(tapply(NewCat$CATCHN*NewCat$WCATCH,list(NewCat$year,NewCat$species),sum))

filter(NewCat,species=='HER' & year==1978)

NewCat$oldNew<-'new';
OldCat$oldNew<-'Old'

if (makeGraphs) {
  a<-rbind(NewCat,subset(OldCat,select=names(NewCat)))
  a$yield<-a$CATCHN*a$WCATCH
  b<-aggregate(yield~species+year+oldNew,sum,na.rm=T,data=a)
  ftable(round(tapply(b$yield,list(b$species,b$year,b$oldNew),sum)))
  head(a)

  #b<-subset(b,!is.na(yield))
  trellis.device(device='png',file=file.path(OutputDir,paste0('CatchWeight_all_sp','.png')),width = 1000, height = 800)

  print(xyplot(yield/1000~as.numeric(year)|species,groups=oldNew,data=b,type='a',scales = list(x=list(relation="same"),y=list(relation='free')),xlab='year', ylab='Catch weight (1000 tonnes)',lwd=2,
               auto.key =   list(space = "right", points = FALSE, lines = TRUE)))

  cleanup()

  b<-aggregate(CATCHN~species+year+age+oldNew,sum,na.rm=T,data=a)
  b2<-aggregate(CATCHN~species+age+oldNew,sum,na.rm=T,data=b)
  b2$tot<-b2$CATCHN; b2$CATCHN<-NULL
  b<-subset(merge(b,b2),tot>0)


  by(b,b$species,function(x) {
    trellis.device(device='png',file=file.path(OutputDir,paste0('CatchN_',x[1,'species'],'.png')),width = 1000, height = 800)
    print(xyplot(CATCHN/1000~as.numeric(year)|
                   paste(species,'Age:',formatC(age, digits = 0, width = 2, format = "f", flag = "0")),groups=oldNew,data=x,type='a',scales = list(x=list(relation="same"),y=list(relation='free')),xlab='year',ylab='Catch number (millions)',lwd=2,
                 auto.key =   list(space = "right", points = FALSE, lines = TRUE)))
    cleanup()

  })


  aa<-subset(a,CATCHN>0)
  aa$SOP<-aa$WCATCH*aa$CATCHN
  b<-aggregate(SOP    ~species+year+age+oldNew,sum,na.rm=T,data=aa)
  b2<-aggregate(CATCHN~species+year+age+oldNew,sum,na.rm=T,data=aa)
  b<-merge(b,b2)
  head(b)
  b$WCATCH<-b$SOP/b$CATCHN


  by(b,b$species,function(x) {
    trellis.device(device='png',file=file.path(OutputDir,paste0('CatchMW_',x[1,'species'],'.png')),width = 1000, height = 800)
    print(xyplot(WCATCH~as.numeric(year)|
                   paste(species,'Age:',formatC(age, digits = 0, width = 2, format = "f", flag = "0")),groups=oldNew,data=x,type='a',scales = "free",xlab='year',ylab='Annual mean weight (kg)',lwd=2,
                 auto.key =   list(space = "right", points = FALSE, lines = TRUE)))
    cleanup()

  })
}


if (makeGraphs) {
  NewBio$oldNew<-'new';
  # head(subset(NewBio,quarter==1 & species=='COD' & age==10),200)
  OldBio$oldNew<-'Old'
  OldBio$PROP_CAT<-1
  OldBio$ICES_WSEA<- -1
  a<-rbind(NewBio,subset(OldBio,select=names(NewBio)))
  b<-subset(a,quarter==1)
  b<-aggregate(PROPMAT~species+year+age+oldNew,sum,na.rm=T,data=b)

  subset(b,species=='NOP' & age==10 )
  #cod10
  #summary(cod10)
  b<-subset(b,PROPMAT>=0)

  by(b,b$species,function(x) {
    trellis.device(device='png',file=file.path(OutputDir,paste0('PropMat_',x[1,'species'],'.png')),width = 1000, height = 800)
    print(xyplot(PROPMAT~as.numeric(year)|
                   paste(species,formatC(age, digits = 0, width = 2, format = "f", flag = "0")),groups=oldNew,data=x,type='a',lwd=2,
                 scales = list(x=list(relation="same"),y=list(relation='free')),xlab='year',ylab='Proportion Mature',
                 auto.key =   list(space = "right", points = FALSE, lines = TRUE)))
    cleanup()
  })

  ftable(round(tapply(NewCat$PROP_CAT,list(NewCat$species,NewCat$year,formatC(NewCat$age,digits=3,width=2,flag="0")),mean),2))


  b<-subset(NewBio,quarter==1 &  ICES_WSEA>0, select=c(species, year, age, quarter,WSEA,ICES_WSEA))
  b1<-b
  b1$weight<-b1$WSEA; b1$oldNew<-'SMS'
  b2<-b
  b2$weight<-b2$ICES_WSEA; b2$oldNew<-'ICES'
  b<-rbind(b1,b2)

  by(b,b$species,function(x) {
    trellis.device(device='png',file=file.path(OutputDir,paste0('Weight_Sea_Q1_',x[1,'species'],'.png')),width = 1000, height = 800)

    print(xyplot(weight~as.numeric(year)|
                   paste(species,formatC(age, digits = 0, width = 2, format = "f", flag = "0")),groups=oldNew,data=x,type='a',lwd=2,
                 scales = list(x=list(relation="same"),y=list(relation='free')),xlab='year',ylab='Mean weight in the sea, ICES annual and SMS quarter 1',
                 auto.key =   list(space = "right", points = FALSE, lines = TRUE)))
    cleanup()
  })


  b<-subset(NewBio,WSEA>0,select=c(species, year, age, quarter,WSEA))
  by(b,b$species,function(x) {
    trellis.device(device='png',file=file.path(OutputDir,paste0('Weight_SMS_Sea_',x[1,'species'],'.png')),width = 1000, height = 800)

    print(xyplot(WSEA~as.numeric(year)|
                   paste(species,formatC(age, digits = 0, width = 2, format = "f", flag = "0")),groups=quarter,data=x,type='a',lwd=2,
                 scales = list(x=list(relation="same"),y=list(relation='free')),xlab='year',ylab='Mean weight in the sea, SMS',
                 auto.key =   list(space = "right", points = FALSE, lines = TRUE)))
    cleanup()
  })

}

# tst<-subset(NewCat,species=='POK'); tapply(tst$CATCHN*tst$WCATCH,list(tst$year),sum)
if (FALSE) {
  sort(unique(NewCat$oldNew))
}

b<-expand.grid(species=unique(NewCat$species),year=unique(NewCat$year),age=unique(NewCat$age),quarter=unique(NewCat$quarter))
b<-merge(x=NewCat,y=b,all.y=T)
b[is.na(b$CATCHN),'CATCHN']<-0
b[is.na(b$WCATCH),'WCATCH']<-0
b<-b[order(b$species,b$year,as.numeric(b$age)+10,b$quarter),]
p<-1
for (i in (2:dim(b)[[1]])) {
  if (b[i,'species'] %in% c('COD','WHG','HAD','POK','MAC','PLE','SOL')) {
    if (b[i,'quarter']==1) p<-b[i,'PROP_CAT'] else b[i,'PROP_CAT']<-p
  }
}

if (FALSE) {
  tst<-filter(b,species=='SPR')
  head(tst,50)
  filter(tst,WCATCH<0)
  filter(tst,CATCHN==0)

  tapply(tst$CATCHN*tst$WCATCH,list(tst$year),sum)
  tapply(tst$CATCHN,list(tst$age),sum)
  tapply(tst$WCATCH,list(tst$age),sum)
  tapply(tst$CATCHN*tst$WCATCH,list(tst$year,tst$quarter),sum)
}


NewCat<-b

#head(subset(NewCat,species=='COD' & year==1989 & age==7),20)
# tst<-subset(NewCat,species=='POK'); tapply(tst$CATCHN*tst$WCATCH,list(tst$year),sum)
# tapply(b$CATCHN*b$WCATCH,list(b$species,b$year),sum)



vpa.code<-       c('COD','WHG','HAD','POK','MAC','HER','NSA','SSA','NOP','SPR','PLE','SOL')
la<-SMS@species.info[SMS@species.info[,'last-age-selec']>0,'last-age']
spOrder<-data.frame(nr=1:length(vpa.code),species=vpa.code,la=la)

fy<-SMS@first.year
ly<-SMS@last.year

a<-merge(NewCat,spOrder)

qPresent<-1   # 1:4

a<-subset(a,quarter %in% qPresent)
x<-a$CATCHN
la<-SMS@species.info[SMS@species.info[,'last-age-selec']>0,'last-age']
b<-tapply(x,list(a$species,a$year,a$quarter,a$age),sum,na.rm=T)
# b
dim(b)

# tst<-subset(a,species=='POK'); tapply(tst$CATCHN*tst$WCATCH,list(tst$year),sum)


### change mean weight using the old MSVPA data for the youngest ages, and scale the old MSVPA mean weight to the used single species mean weights


old<-read.table(file=file.path(root,"Data_NorthSea","Old_data",'NS_1998_all_sp','VPA_B01.DAT'),header=T)
a<-subset(old, year==1974 & age <=10 & ((age>=1) | (age==0 & quarter %in% c(3,4))),select=c(-year))

ftable(round(tapply(a$WSEA,list(a$species,a$quarter,a$age),sum),3))

# mean w
b<-aggregate(WSEA~species+age,mean,data=a)
head(b)
b$meanWSEA<-b$WSEA; b$WSEA<-NULL
head(b)
subset(b,species=='COD' & age==3 )

# minimum W
b2<-aggregate(WSEA~species+age,min,data=a)
head(b2)
b2$minWSEA<-b2$WSEA; b2$WSEA<-NULL
head(b2)
subset(b2,species=='COD' & age==3 )

b<-merge(b,b2)
subset(b,species=='COD' & age==3 )

a<-subset(a,species %in% c('COD','WHG','HAD','POK','MAC','HER','PLE','SOL'),select=c(species,quarter,age,WSEA))
head(a)
subset(a,species=='WHG' & age==1 )

ab<-merge(x=a,y=b,all.x=T)
subset(ab,species=='COD' & age==3 )  %>% arrange(age,quarter)
ab$wFac<-ab$WSEA/ab$meanWSEA
ab$wFacMin<-ab$WSEA/ab$minWSEA

subset(ab,species=='COD' & age==3 ) %>% arrange(age,quarter)

oldwsea<-subset(ab, select=c(species, quarter,age, WSEA))
oldwsea$oldWsea<-oldwsea$WSEA; oldwsea$WSEA<-NULL
head(oldwsea)
subset(oldwsea,species=='COD' & age==3 )

NewBio2<-merge(x=NewBio,y=oldwsea,all.x=T)
head(NewBio2)


ab<-subset(ab, species=='COD' & age>=3 | species=='WHG' & age >=1 | species=='HAD' & age>2 | species=='POK' & age>3 | species=='MAC' & age>2 |
             species=='HER' & age>2 | species=='PLE' & age >0 | species=='SOL' & age >0,
           select=c(species,quarter,age,wFac,wFacMin))
head(ab)
subset(ab,species=='WHG')

#ftable(round(tapply(ab$wFac,list(ab$species,ab$quarter,ab$age),sum),3))
#ftable(round(tapply(NewBio2$ICES_WSEA,list(NewBio2$species,NewBio2$quarter,NewBio2$age),mean),3))


b<-merge(x=NewBio2,y=ab,all.x=TRUE)
sort(unique(b$species))

bmean<-subset(b,!(species %in% c('WHG')))
sort(unique(bmean$species))


bmin<-subset(b,species %in% c('WHG'))
sort(unique(bmin$species))

#subset(bmean,year==1978 & species=='COD' & age==4)
bmean$WSEA<-ifelse(is.na(bmean$wFac) | is.na(bmean$ICES_WSEA),bmean$WSEA, ifelse(bmean$ICES_WSEA<=0,bmean$WSEA,bmean$wFac*bmean$ICES_WSEA))
bmean[bmean$WSEA<1.12 & bmean$species=='COD' & bmean$age=='3','WSEA']<-1.125
#subset(bmean,year==1978 & species=='COD' & age==4)

subset(bmin,year==1983 & species=='WHG' & age==4)
bmin$WSEA<-ifelse(is.na(bmin$wFacMin) | is.na(bmin$ICES_WSEA),bmin$WSEA, ifelse(bmin$ICES_WSEA<=0,bmin$WSEA,bmin$wFacMin*bmin$ICES_WSEA))
subset(bmin,year==1983 & species=='WHG' & age==4)
b78<-subset(bmin,year==1978 & species=='WHG')
bb78<-subset(b78,select=c(species,quarter,age,WSEA))
bb78$newWSEA<-bb78$WSEA; bb78$WSEA<-NULL

b7477<-subset(bmin,year<1978 & species=='WHG')
head(bb78)
bb<-merge(x=b7477,y=bb78,all.x=TRUE)

bb$WSEA<-bb$newWSEA; bb$newWSEA<-NULL;
bmin<-subset(bmin,!(year<1978 & species=='WHG'))
bmin<-rbind(bb,bmin)
sort(unique(bmin$species))
#ftable(round(tapply(bmin$WSEA,list(bmin$year,bmin$quarter,as.numeric(bmin$age)+10),sum),3))

b<-rbind(bmean,bmin)

b[b$WSEA<=0,'WSEA']<-NA
ftable(round(tapply(b$WSEA,list(b$species,b$quarter,as.numeric(b$age)+10),mean),3))
head(b)
head(NewBio)

NewBio<-subset(b,select=names(NewBio))

# b<-subset(NewBio,species=='WHG')
# ftable(round(tapply(b$WSEA,list(b$year,b$quarter,as.numeric(b$age)+10),mean),3))

if (makeGraphs) {
  b<-subset(NewBio,WSEA>0,select=c(species, year, age, quarter,WSEA))
  by(b,b$species,function(x) {
    trellis.device(device='png',file=file.path(OutputDir,paste0('Weight_Final_Sea_',x[1,'species'],'.png')),width = 1000, height = 800)

    print(xyplot(WSEA~as.numeric(year)|
                   paste(species,formatC(age, digits = 0, width = 2, format = "f", flag = "0")),groups=quarter,data=x,type='a',lwd=2,
                 scales = list(x=list(relation="same"),y=list(relation='free')),xlab='year',ylab='Mean weight in the sea, SMS',
                 auto.key =   list(space = "right", points = FALSE, lines = TRUE)))
    cleanup()
  })
}



if (makeGraphs) {
  NewBio$oldNew<-'key-2023';
  # head(subset(NewBio,quarter==1 & species=='COD' & age==10),200)
  OldBio$oldNew<-'key-2020'
  OldBio$PROP_CAT<-1
  OldBio$ICES_WSEA<- -1

  a<-rbind(NewBio,subset(OldBio,select=names(NewBio)))
  subset(a,species=='WHG' & year==2016 & age==4)
  b<-subset(a,quarter==1)
  b<-aggregate(WSEA~species+year+age+oldNew,sum,na.rm=T,data=b)
  b<-subset(b,WSEA>0)


  by(b,b$species,function(x) {
    trellis.device(device='png',file=file.path(OutputDir,paste0('Weights_sea_old_new_Q1',x[1,'species'],'.png')),width = 1000, height = 800)
    print(xyplot(WSEA~as.numeric(year)|
                   paste(species,formatC(age, digits = 0, width = 2, format = "f", flag = "0")),groups=oldNew,data=x,type='a',lwd=2,
                 scales = list(x=list(relation="same"),y=list(relation='free')),xlab='year',ylab='Weight in the Sea Quarter 1',
                 auto.key =   list(space = "right", points = FALSE, lines = TRUE)))
    cleanup()
  })


  NewBio$oldNew<-'key-2023';
  # head(subset(NewBio,quarter==1 & species=='COD' & age==10),200)
  OldBio$oldNew<-'key-2020'
  OldBio$PROP_CAT<-1
  OldBio$ICES_WSEA<- -1
  a<-rbind(NewBio,subset(OldBio,select=names(NewBio)))

  b<-aggregate(WSEA~species+year+age+quarter+oldNew,sum,na.rm=T,data=a)
  b<-subset(b,WSEA>0)


  by(b,b$species,function(x) {
    trellis.device(device='png',file=file.path(OutputDir,paste0('Weights_sea_old_new_Q1_4',x[1,'species'],'.png')),width = 1000, height = 800)
    print(xyplot(WSEA~as.numeric(year)|
                   paste(species,formatC(age, digits = 0, width = 2, format = "f", flag = "0"),"Q:",quarter),groups=oldNew,data=x,type='a',lwd=2,
                 scales = list(x=list(relation="same"),y=list(relation='free')),xlab='year',ylab='Weight in the',
                 auto.key =   list(space = "right", points = FALSE, lines = TRUE)))
    cleanup()
  })

  by(b,b$species,function(x) {
    trellis.device(device='png',file=file.path(OutputDir,paste0('Weights_sea_old_new_Q1_4_sameScale',x[1,'species'],'.png')),width = 1000, height = 800)
    print(xyplot(WSEA~as.numeric(year)|
                   paste(species,formatC(age, digits = 0, width = 2, format = "f", flag = "0"),"Q:",quarter),groups=oldNew,data=x,type='a',lwd=2,
                 scales = list(x=list(relation="same"),y=list(relation='same')),xlab='year',ylab='Weight in the Sea',
                 auto.key =   list(space = "right", points = FALSE, lines = TRUE)))
    cleanup()
  })

}


#####################################################################################
# make basis input catch file for 2023 keyrun
if (FALSE) {
  tst<-NewCat
  dim(tst)
  dups<- tst %>% select(year ,age ,species ,quarter )
  dim( dups %>% unique)
  t2<-tst %>% group_by(species,year,quarter,age) %>% mutate(n=dplyr::n()) %>%
      filter(n>1) %>% arrange(species,year,quarter,age)
  print(t2,n=50)
  print(filter(t2,CATCHN>0),n=50)
  sort(unique(t2$species))
  tst$dups<-duplicated(tst)
  tst
}
CAT.01<-NewCat
save(CAT.01,file=file.path(root,exchangeDir,'CAT_01.Rdata'))


# Biological data (Wsea,propmat, M and M1) and an new variable (proportion of the stock within area 4 (the North Sea))
BIO.01<-NewBio

# subset(NewBio,species=='NOP' & age==1 & year==1987)

save(BIO.01,file=file.path(root,exchangeDir,'BIO_01.Rdata'))


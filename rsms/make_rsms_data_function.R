make_rsms_data<-function(dir,sms.dat='sms.dat',adj_san0=TRUE,seasFfrom=c('F','catch')[1],sepSan=FALSE) {
 if (FALSE) {
   dir=data.path; sms.dat='rsms.dat'; seasFfrom<-'catch'; adj_san0=TRUE;seasFfrom=c('F','catch')[2]; 
 }

sms<-read.RSMS.control(dir,file=sms.dat)

multi<-!(sms@VPA.mode==0)

info<-sms@species.info
multi.environment<- any(info[,'predator']==2)

nSpecies<-sms@no.species-first.VPA+1L  # species with analytically assessment
nOthSpecies<-0L
ages<-sms@first.age:sms@max.age.all
nAges<-length(ages)
recSeason<-as.integer(sms@rec.season)
nSeasons<-as.integer(sms@last.season)
years<-sms@first.year.model:sms@last.year.model; nYears<-length(years)

off.age<- as.integer(1-sms@first.age)
off.year<- - as.integer(sms@first.year.model-1)
stockRecruitmentModelCode<-as.integer(sms@SSB.R)

info<-sms@species.info[first.VPA:sms@no.species,c("last-age", "first-age F>0",  "+group", "SpawningQ"),drop=FALSE]
info<-cbind(info,fModel=sms@fModel)
fSepar<-if_else(sms@fModel==2,sms@firstAgeYearEffect+off.age,99L)
info<-cbind(info,fSepar)
spNames<-dimnames(info)[[1]]
othspNames<-dimnames(info)[[1]][1:(first.VPA-1)]

names(stockRecruitmentModelCode)<-spNames
use_rho<-as.integer(sms@use_rho)
rec_loga <- rep(1,nSpecies)
rec_loga[stockRecruitmentModelCode==1] <- log(200)  # Ricker
rec_loga[stockRecruitmentModelCode==2] <- 4         # B&W

rec_logb <- rep(1,nSpecies)
rec_logb[stockRecruitmentModelCode==1] <- -12  # Ricker
rec_logb[stockRecruitmentModelCode==2] <- -12        # B&W

if (multi.environment) off.species<-as.integer(-first.VPA+1) else off.species<-0L

info<-cbind(info,s=(first.VPA:sms@no.species)+off.species)
off.oths<-0L #other species offset (for now)
off.season<-0L # not used 

fbarRange<-matrix(as.integer(sms@avg.F.ages),ncol=2,dimnames=dimnames(sms@avg.F.ages))+off.age;
seasonalCatches<-as.integer(sms@combined.catches==0)

## configuration of parameter keys

rho<-rep(0.7,nSpecies)

##### states at age for F random walk,
# use sms@keyLogFsta

x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nSpecies,dimnames=list(spNames,paste('age',ages)))
i<-1L
xx<-NULL
for (s in 1:nSpecies) {
  i<-1L
  la<-info[s,"last-age"]
  aa<-sort(unique(c(sms@keyLogFsta[[s]],la+1)))
  for (j in (1:(length(aa)-1))) {
    for (age in aa[j]:(aa[j+1]-1)) {
      x[s,age+off.age ]<-i;
      xx<-rbind(xx,data.frame(s=s,age=age,a=age+off.age,keyLogFsta=i))
    }
    i<-i+1L
  }
}

keyLogFsta<-x
keyLogFsta.df<-xx

keyLogFsta.list<-by(keyLogFsta.df,keyLogFsta.df$s,function(x) as.vector(x$keyLogFsta)-x[1,"keyLogFsta"]+1 )


nlogF<-unlist(lapply(sms@keyLogFsta,length))
nlogFto<-cumsum(nlogF)  
nlogFfrom<- c(1,head(cumsum(nlogF),-1)+1 )
nlogFfromTo<-matrix(c(as.integer(nlogFfrom),as.integer(nlogFto)),ncol=2)


keyLogFstaSd<-keyLogFsta
keyLogFstaSd[keyLogFstaSd<0]<- -9999L
keyLogFstaSd<-keyLogFstaSd+nlogFfrom-1
keyLogFstaSd[keyLogFstaSd<0]<- -1L
keyLogFstaSd
logSdLogFsta<-rep(-0.7,max(keyLogFstaSd))


##### year and Ages in separable F-Model
# use sms@catch.sep.year and sms@catch.sep.age

x<-map2(sms@catch.sep.age,sms@catch.sep.year,function(x,y) {
  xx<-expand.grid(age=x,year=y)
})

for (s in 1:nSpecies) x[[s]]<-data.frame(x[[s]],s=s,fModel=info[s,'fModel'])
ay<-do.call(rbind,x) %>% filter(fModel==2) %>%arrange (s,age,year) %>%
  mutate(index=1:dplyr::n(), y=year-sms@first.year.model+1,a=age-sms@first.age+1)

x<-lapply(1:nSpecies,function(x) {
  matrix(-1L,ncol=info[x,'last-age']-sms@first.age+1L,nrow=nYears,dimnames=list(years,paste('age',sms@first.age:info[x,'last-age'])))
})

for (i in (1:dim(ay)[[1]])) {
  x[[ay[i,'s']]][ay[i,'y']:nYears,ay[i,'a']]<-ay[i,'index'] 
}

x<-lapply(x,function(xx){
  for (a in 2:ncol(xx)) if (xx[1,a]== -1) xx[,a]<- xx[,a-1]
  xx
})

keyLogSeparF<- x
logSeparF<-rep(0,max(ay$index))

##### states at age for N random walk,

nlogN<-as.integer(info[,"last-age"]-sms@first.age+1L)

#### Process Noise N
# use sms@keyVarLogN

x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nSpecies,dimnames=list(spNames,paste('age',ages)))
i<-1L
xx<-NULL
for (s in 1:nSpecies) {
  #i<-1L
  la<-info[s,"last-age"]
  aa<-sort(unique(c(sms@keyVarLogN[[s]],la+1)))
  for (j in (1:(length(aa)-1))) {
    for (age in aa[j]:(aa[j+1]-1)) {
      x[s,age+off.age ]<-i;
      xx<-rbind(xx,data.frame(s=s,age=age,a=age+off.age,keyVarLogN=i))
    }
    i<-i+1L
  }
}

keyVarLogN<-x
keyVarLogN.df<-xx

x<-NULL
for (s in 1:nSpecies) {
  for (a in sms@keyVarLogN[[s]]) if (a==fa) x<-c(x,0.35) else x<-c(x,-0.35) 
}
logSdLogN<-x
# logSdLogN


nlogNto<-cumsum(nlogN)  
nlogNfrom<- c(1,head(cumsum(nlogN),-1)+1L )
nlogNfromTo<-matrix(c(as.integer(nlogNfrom),as.integer(nlogNto)),ncol=2)


#keyLogNsta[,1:2]<-c(seq(from=1,by=2,length.out=length(nlogN)),seq(from=2,by=2,length.out=length(nlogN)))   
#keyLogNsta

#### variance of catch and surveys


# sms@catch.s2.group
x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nSpecies,dimnames=list(spNames,paste('age',ages)))
i<-1L
xx<-NULL
for (s in 1:nSpecies) {
  la<-info[s,"last-age"]
  aa<-sort(unique(c(sms@catch.s2.group[[s]],la+1)))
  for (j in (1:(length(aa)-1))) {
    for (age in aa[j]:(aa[j+1]-1)) {
      x[s,age+off.age ]<-i;
      xx<-rbind(xx,data.frame(s=s,age=age,a=age+off.age,keyVarObsCatch=i))
    }
    i<-i+1L
  }
}
keyVarObsCatch<-x
keyVarObsCatch.df<-xx
logSdLogObsCatch<-rep(-0.35,max(keyVarObsCatch)) # parameter
 
###########  survey

a<-readNewFleetInfo(dir=dir,of='new_fleet_info.dat',off.age=off.age,off.year=off.year,verbose=FALSE) 
keySurvey<-a[['k']]
notUse<-keySurvey[keySurvey[,'useFleet']==0,'f']
catchability<-a[['catchability']]
catchability.var<-a[['catchability.var']]

if (length(notUse>0)) {
  catchability[notUse]<-NULL
  catchability.var[notUse]<-NULL
  keySurvey<-keySurvey[-notUse,]
  keySurvey[,'f']<-1L:dim(keySurvey)[[1]]
}

nFleets<-nrow(keySurvey)
fleetNames<-rownames(keySurvey)
 
x<-matrix(0.0,ncol=1L,nrow=nFleets,dimnames=list(paste(1:nFleets,fleetNames),c("sampleTimeWithin")))
x[,"sampleTimeWithin"]<- (as.numeric(keySurvey[,'endf'])-as.numeric(keySurvey[,'startf']))/2
sampleTimeWithin<-x

keySurvey.df<-as.data.frame(keySurvey) %>% tibble::rownames_to_column("fName")

## survey variance

x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nFleets,dimnames=list(fleetNames,paste('age',ages)))
i<-1L
xx<-NULL
for (f in (1:nFleets)) {
  jjj<-sort(unique(c(catchability.var[[f]],keySurvey[f,'maxage']+1)+off.age))
  s<-keySurvey[f,'s']
  for (jj in (1:(length(jjj)-1))) {
    for (j in jjj[jj]:(jjj[jj+1]-1)) {
      x[f,j]<-i;
      xx<-rbind(xx,data.frame(f=f,s=s,a=j,age=j-off.age,keyVarObsSurvey=i))
    }
    i<-i+1L
  }
}

keyVarObsSurvey<-x
keyVarObsSurvey.df<-xx

logSdLogObsSurvey<-rep(-0.35,max(keyVarObsSurvey))

## survey catchability
x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nFleets,dimnames=list(fleetNames,paste('age',ages)))
i<-1L
xx<-NULL
for (f in (1:nFleets)) {
  jjj<-sort(unique(c(catchability[[f]],keySurvey[f,'maxage']+1)+off.age))
  s<-keySurvey[f,'s']
  for (jj in (1:(length(jjj)-1))) {
    for (j in jjj[jj]:(jjj[jj+1]-1)) {
      x[f,j]<-i;
      xx<-rbind(xx,data.frame(f=f,s=s,a=j,age=j-off.age,keyCatchability=i))
    }
    i<-i+1L
  }
}
keyCatchability<-x
keyCatchability.df<-xx

logCatchability<-rep(0.0,max(keyCatchability)) 

## Catchability power age
x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nFleets,dimnames=list(fleetNames,paste('age',ages)))
i<- 0L
xx<-NULL
for (f in (1:nFleets)) {
  s<-keySurvey[f,'s']
  for (j in (keySurvey[f,'mina']:keySurvey[f,'maxa'])) {
    if (j==keySurvey[f,'PowerAge']) {i<-i+1L; x[f,j ]<-i} 
    xx<-rbind(xx,data.frame(f=f,s=s,species.n=s,age=j-off.age,a=j,keyPowerQ=i))
  }
}

keyQpow<-x
keyQpow.df<-xx
if (max(keyQpow) >0) logQpow<-rep(0.0,max(keyQpow)) else  logQpow<-rep(0.0,0)


### survey observations
cpue<-readFleetCatch(dir,of='new_fleet_info.dat',so='fleet_catch.in', condense=TRUE)

if (FALSE) {
  out<-by(cpue,cpue$f,function(x){
    (xtabs(obs~age+year,data=droplevels(x)))
  })
  out[[16]]
  out[[17]]
}

eff_fl<-keySurvey[keySurvey[,'type'] %in% c(4,5),,drop=FALSE]

if (dim(eff_fl)[[1]] > 0) { #inflate data with effort for each age group for type=effort
  ef<-  eff_fl[,'f']
  cpue_ok<-filter(cpue,!(f %in% ef)) 
  cpue_eff<-filter(cpue,f %in% ef) 
  a<-lapply(ef,function(ff){
    a<-filter(cpue_eff,f==ff)
    ages<-data.frame(age=-999,newage=keySurvey[ff,"minage"]:keySurvey[ff,"maxage"])
    left_join(a,ages,by = join_by(age),relationship = "many-to-many")  %>% mutate(age=newage,newage=NULL)
  })
  cpue_eff<-do.call(rbind,a)
  cpue<-rbind(cpue_ok,cpue_eff)
}

if (FALSE) {
  out<-by(cpue,cpue$f,function(x){
    (xtabs(obs~age+year,data=droplevels(x)))
  })
  out[[16]]
  out[[17]]
}
cpue<-filter(cpue,!(is.na(obs)))
cpue<-cpue %>% filter(year<=sms@last.year.model)   # SKAL RETTES


cpue<-cpue %>% mutate(y=year+off.year,a=age+off.age,species.n=s,season=q,fleet.no=f) %>%
  dplyr::select(year,y,season,q,age,a,species.n,s,fleet.no,f,obs) %>%
  filter(obs>0) %>% arrange(s,fleet.no,y,a) 
 cpue$obs.no<-1L:dim(cpue)[[1]]
#head(cpue)
 if (FALSE) {
   out<-by(cpue,cpue$f,function(x){
     (xtabs(obs~age+year,data=droplevels(x)))
   })
   out[[16]]
   out[[17]]
   out[[18]]
   out[[19]]
   out[[20]]
 }

fobs<-cpue %>% select(s,f,y,a) %>% group_by(s,f) %>%summarize(n=dplyr::n()) %>%  ungroup() %>% mutate(last=cumsum(n),first=lag(last)+1L)
fobs[1,"first"]<-1L


keySurvey.overview<-cbind(keySurvey,n=as.vector(fobs$n),first=as.vector(fobs$first),last=as.vector(fobs$last))

sumCreep<-sum(keySurvey.overview[,'techCreep'])
logTechCreep<-rep(0.0,sumCreep)


#head( keyVarObsSurvey.df); filter(keyVarObsSurvey.df,is.na(keyVarObsSurvey))
#head(keyCatchability.df)
#head(keyQpow.df)

k<-full_join(keyVarObsSurvey.df,keyCatchability.df) %>% dplyr::select(f, s, a, keyVarObsSurvey, keyCatchability) %>% arrange(f,s,a)
k<-full_join(k,select(keyQpow.df,s,f,a,keyPowerQ),by = join_by(f, s, a))

cpue<-left_join(cpue,k,by = join_by(a, s, f)) %>% arrange(f,s,a)
cpue$obs.no<-1:dim(cpue)[[1]]

if (FALSE) {
  out<-by(cpue,cpue$f,function(x){
    (xtabs(obs~age+year,data=droplevels(x)))
  })
  out[[16]]
  out[[17]]
}
keySurvey<-  cpue %>% dplyr::select(obs.no,f,s,y,a,q ,keyVarObsSurvey, keyCatchability,keyPowerQ)
summary(filter(keySurvey, is.na(keyVarObsSurvey) | is.na(keyCatchability)))
filter(keySurvey, is.na(keyVarObsSurvey) | is.na(keyCatchability))

stopifnot(dim(filter(keySurvey, is.na(keyVarObsSurvey) | is.na(keyCatchability)))[[1]]==0)
keySurvey<-keySurvey %>% mutate_if(is.numeric,as.integer)

keySurvey<-as.matrix(keySurvey)

logSurveyObs<-log(cpue$obs)


x<-scan(file=file.path(dir,'recruitment_years.in'),comment.char='#',quiet=TRUE)
stopifnot(tail(x,1)== -999)
rec_years<-matrix(head(x,-1),nrow=nSpecies,byrow=TRUE)==1

source(file.path(rsms.root.prog,"from_sms_format_to_rsms_data.R"))

d<-From_SMS_format_to_rsms(otherPredExist=multi.environment,catchMultiplier=1,dir=dir)

#merge data
b<-full_join(d$catch,d$bio,join_by(year, species.n, quarter, sub_area, age))
if (multi) {
  b2<-full_join(d$mean_l,d$consum,join_by(year, species.n, quarter, sub_area, age))
  b<-full_join(b,b2,join_by(year, species.n, quarter, sub_area, age))
}

a<-data.frame(species.n=first.VPA:nsp,annualC=sms@combined.catches)
b<-left_join(b,a,by = join_by(species.n))
#head(b)

x<-b$age==fa & b$quarter < recSeason & b$annualC==0
b[x,'CATCHN']<-0
b[x,'WCATCH']<-0


b0<-filter(b,annualC==0) %>% group_by(year,species.n, sub_area, age) %>% mutate(seasFprop=CATCHN/sum(CATCHN,na.rm=T)) %>% ungroup() %>%
  select(year, species.n, quarter, sub_area,   age, annualC, seasFprop)  
b1<-filter(b,annualC==1) %>%  select(year, species.n, quarter, sub_area,age, annualC) %>% mutate(seasFprop=0.25) %>%
     mutate(seasFprop=if_else(age==fa & quarter >= recSeason,0.5, seasFprop)) %>%
     mutate(seasFprop=if_else(age==fa & quarter < recSeason,0.0, seasFprop))
#filter(b1,year==1975 & age==0)
bb<-rbind(b0,b1) 

check<-bb %>% group_by(year,species.n, sub_area, age) %>%summarize(Fprop=sum(seasFprop)) %>% ungroup() 
# filter(check,!(Fprop>0.999 & Fprop<1.001) )
stopifnot(all(check$Fprop>0.999 & check$Fprop<1.001 | is.na(check$Fprop)))

b<-left_join(b,bb,by = join_by(year, species.n, quarter, sub_area, age,annualC))

b<-b %>% mutate(y=year+off.year,q=quarter,s=species.n+off.species,a=age+off.age)
b<-left_join(b,data.frame(s=info[,'s'],la=info[,"last-age"]+off.age),by = join_by(s)) %>% filter(a<=la)


propMat         <-by(b,b$s,function(x) {y<-tapply(x$PROPMAT,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
stockMeanWeight <-by(b,b$s,function(x) {y<-tapply(x$WSEA,list(x$y,x$q,x$a),sum) ; y[is.na(y)]<-0; y})
catchMeanWeight <-by(b,b$s,function(x) {y<-tapply(x$WCATCH,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
natMor          <-by(b,b$s,function(x) {y<-tapply(x$M,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
catchNumber     <-by(b,b$s,function(x) {y<-tapply(x$CATCHN,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})

# round(ftable(catchNumber[[1]][,,]/1000),1)



if (seasFfrom==c('F','catch')[2]) seasFprop <-by(b,b$s,function(x) {y<-tapply(x$seasFprop,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y}) else {
    seasFprop <-by(b,b$s,function(x) {y<-tapply(x$propF,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
}


# round(seasFprop[[1]][,,3],2) # cod age 3-1
# round(seasFprop[[6]][,,3],2) #her age 3-1
# round(seasFprop[[7]][,,3],2) #age 3-1=2

if (adj_san0 ) {
  sa<-grep('sandeel',spNames)
  if (length(sa) >0) {
    if (nSeasons==2) for (i in sa)  seasFprop[[i]][,,1]<-rep(c(0,1),each=nYears)  # age 0
    if (nSeasons==4) for (i in sa)  seasFprop[[i]][,,1]<-rep(c(0,0,1,0),each=nYears)  # age 0
  }
}
# round(ftable(seasFprop[[1]][,,]),2)

catchNoSeason   <-by(b,b$s,function(x) {y<-tapply(x$CATCHN,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})

zero<-propMat; for (s in 1:nSpecies) {zero[[s]][,,]<-0}


#### catch observations

#first zero-catch years
zy<-data.frame(s=1:nSpecies,y=0L)
if (sms@zero.catch.year.season==1) {
  fname<-'zero_catch_year_season.in'
  z<-head(scan(file.path(dir,fname),comment.char='#',quiet = TRUE),-1)
  bz<-expand.grid(s=(first.VPA:nsp) +off.species,y=1:nYears,q=1:nSeasons)
  bz<-bz[order(bz$s,bz$y,bz$q),]
  bz$z<-z 
  zeroCatchYear<-bz %>%group_by(s,y) %>% summarize(z=sum(z)) %>%  filter(z==0) %>% mutate(z=NULL)
} else zeroCatchYear<-NULL
zy<-rbind(zy,zeroCatchYear)

zeroCatchYear<-by(zy,list(zy$s),function(x) { x<-x$y[x$y>0]})
lapply(zeroCatchYear,function(x) length(x))
zeroCatchYearExistsSp<-unlist(lapply(zeroCatchYear,function(x) length(x)))
zeroCatchYearExistsSp[zeroCatchYearExistsSp>0] <-1L

zy<-zy %>% mutate(zero=TRUE) %>% filter(y>0)
if (dim(zy)[[1]]>0) zeroCatchYearExists<-1L else zeroCatchYearExists<-0L

catch<-aggregate(CATCHN~year+species.n+ sub_area+ age+y+ s+ a,FUN=sum,data=subset(b,CATCHN>=0 )) %>% filter(CATCHN>0)
catch<-left_join(catch,data.frame(s=info[,'s'],faf=info[,"first-age F>0"]+off.age),by = join_by(s)) %>% filter(a>=faf)
catch<-left_join(catch,zy,by = join_by(y, s)) %>% filter(is.na(zero)) %>% mutate(zero=NULL)

catch<-catch %>% arrange(s,y,a)
catch$obs.no<-1:dim(catch)[[1]]
logCatchObs<-log(catch$CATCHN)

keyCatch<-catch %>% mutate(CATCHN=NULL) %>% mutate_if(is.numeric,as.integer)
k<-full_join(keyLogFsta.df,keyVarObsCatch.df,by = join_by(s, age, a))
keyCatch<-left_join(keyCatch,k,by = join_by(s, age, a))
keyCatch<-keyCatch %>% arrange(obs.no)

stopifnot(dim(filter(keyCatch,is.na(keyLogFsta)))[[1]]==0)
keyCatch<-as.matrix(keyCatch)


info<-cbind(info,
    la=info[,"last-age"]+off.age,
    faf=info[,"first-age F>0"]+off.age)
 


Un<-matrix(0.0, nrow=sum(info[,"last-age"]-sms@first.age+1L),ncol=nYears )
Uf<-matrix(0.0, nrow=length(unlist(sms@keyLogFsta)),ncol=nYears)


#return values
list(
   data=list(                                          # Description of data, most of the text is copied from the ?sam.fit description
     sms.mode=sms@VPA.mode,                            #
     info=info,                                        # Various information for each species with analytically assessment
     nSpecies=nSpecies,                                # Number of species with analytically assessment
     #nOthSpecies=nOthSpecies,                          # Number of other predators (used for multispecies mode)   
     nYears=nYears,                                    # Number of years used in the model
     nAges=nAges,
     minAge=as.integer(sms@first.age),                 # A vector of integers defining the the lowest age class in the assessment for each species.
     #maxAge=as.integer(sms@max.age.all),               # Maximum age for all species. The actual age range by species is given in "info"
     recAge=as.integer(sms@first.age),                 # recruitmet age
     years=sms@first.year.model:sms@last.year.model,   # A vector of the years used in the model
     spNames=spNames,                                  # Species names of species with analytically assessment
     fleetNames=fleetNames,                            # names of survey fleets
     nSeasons=nSeasons,                                # number of seasons
     spawnSeason=1L,
     recSeason=recSeason,                              # recruitment season
     fbarRange=fbarRange,                              # Minimum and maximum age by species used to calculate Fbar 
     off.age=off.age,                                  # offset between age and the age index used in RSMS (where fist age index  have to be 1 )
     off.year=off.year,                                # offset between year and the year index used in RSMS (where fist year index  have to be 1 )
     #off.season=off.season,
     #off.species=off.species,
     #off.oths=off.oths,                          
     seasonalCatches=seasonalCatches,                  # are seasonal cataches available
     useRho=use_rho,                                   # use correlation in F at age
     stockRecruitmentModelCode=stockRecruitmentModelCode,
     zeroCatchYearExists=zeroCatchYearExists,
     zeroCatchYearExistsSp=zeroCatchYearExistsSp,
     zeroCatchYear=zeroCatchYear,
     keyLogFsta=keyLogFsta,                         # A matrix of integers. The number of rows is equal to the number of species fleets and the number of columns is equal to the number of age classes. The matrix describes the coupling of the fishing mortality states. '-1' is used for entries where no fishing mortality applies. For the valid entries consecutive integers starting 1 must be used, because they are used as indices in the corresponding state vector. If the same number is used for two fleet-age combinations, then the fishing mortality for those are assumed equal (linked to the same state).
     nlogF=nlogF,                                   # A vector with the sum of number of F state age groups by species. Each species number corresponds to the  number of unique (and not -1 value) in each (species) row of keyLogFsta
     keyLogFstaSd=keyLogFstaSd,                     # A matrix which links each F state age group by species to the logSdLogFsta=
     nlogFfromTo=nlogFfromTo,
     nlogN=nlogN,                                   # A vector with the sum of number of N state age groups by species
     nlogNfromTo=nlogNfromTo,
     
     keyVarLogN=keyVarLogN,
     keyVarObsCatch=keyVarObsCatch,
     keyVarObsSurvey=keyVarObsSurvey,
     keyCatchability=keyCatchability,
     propMat=  propMat,      
     stockMeanWeight=stockMeanWeight,
     catchMeanWeight=catchMeanWeight,
     catchNumber=catchNumber,
     seasFprop=seasFprop,
     natMor= natMor,
     propF=zero,
     propM=zero,
     keyCatch=keyCatch,
     catchNoSeason=catchNoSeason,
     logCatchObs=logCatchObs,
     recruitYears=rec_years,
     #nCatchObs=length(logCatchObs),
     keyLogSeparF=keyLogSeparF,
     keySurvey=keySurvey,
     keySurvey.overview=keySurvey.overview,  # not really used
     logSurveyObs=logSurveyObs,
     sampleTimeWithinSurvey=sampleTimeWithin
     #nSurveyObs=length(logSurveyObs)
   ),
   
  parameters=list(
      logSdLogObsCatch=logSdLogObsCatch,
      logCatchability=logCatchability,
      logSdLogObsSurvey=logSdLogObsSurvey,
      logQpow=logQpow,
      logSeparF=logSeparF,
      logTechCreep=logTechCreep,
      logSdLogFsta=logSdLogFsta,
      logSdLogN=logSdLogN,
      rho=rho,
      rec_loga=rec_loga, 
      rec_logb=rec_logb, 
      Un=Un,
      Uf=Uf
    ))
}



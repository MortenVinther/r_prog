make_rsms_data<-function(dir,annual=FALSE) {
# dir<-"S16_S21"; annual<-TRUE; dir="S19"
 
Init.function(dir=file.path(root,dir)) # initialize SMS environment
cat(SMS.control@species.names,'\n') # just cheking
  
rsms.root<-file.path("~","cod","RSMS");
sam.root<-file.path("~","cod");

sms<-SMS.control  # just shorter name
multi<-!(sms@VPA.mode==0)

info<-sms@species.info
multi.environment<- any(info[,'predator']==2)

nSpecies<-sms@no.species-first.VPA+1L  # species with analytical assessment
nOthSpecies<-0L
ages<-sms@first.age:sms@max.age.all
nAges<-length(ages)
recSeason<-as.integer(sms@rec.season)
nSeasons<-as.integer(sms@last.season)
years<-sms@first.year.model:sms@last.year.model; nYears<-length(years)
info<-sms@species.info[first.VPA:sms@no.species,c("last-age", "first-age F>0", "last-age-selec", "last-age-likelihood", "+group", "SSB/R"),drop=FALSE]
spNames<-dimnames(info)[[1]]
othspNames<-dimnames(info)[[1]][1:(first.VPA-1)]
off.age<- as.integer(1-sms@first.age)
off.year<- - as.integer(sms@first.year.model-1)
stockRecruitmentModelCode<-if_else (info[,"SSB/R"] >100,4L,info[,"SSB/R"])
# SKAL RETTES
stockRecruitmentModelCode[stockRecruitmentModelCode==4]<-2
stockRecruitmentModelCode[stockRecruitmentModelCode==3]<-1
stockRecruitmentModelCode[1]<-0
stockRecruitmentModelCode<-as.integer(stockRecruitmentModelCode)
names(stockRecruitmentModelCode)<-spNames

rec_loga <- rep(1,nSpecies)
rec_loga[stockRecruitmentModelCode==1] <- log(200)  # Ricker
rec_loga[stockRecruitmentModelCode==2] <- 1         # B&W

rec_logb <- rep(1,nSpecies)
rec_logb[stockRecruitmentModelCode==1] <- -12  # Ricker
rec_logb[stockRecruitmentModelCode==2] <- -12        # B&W

if (multi.environment) off.species<-as.integer(-first.VPA+1) else off.species<-0L

info<-cbind(info,s=(first.VPA:sms@no.species)+off.species)
off.oths<-0L #other species offset (for now)
off.season<-0L # not used 

fbarRange<-matrix(as.integer(sms@avg.F.ages),ncol=2,dimnames=dimnames(sms@avg.F.ages));

load(file=file.path(sam.root,"sam_par_dat.Rdata"),verbose=TRUE)
#str(sam_data)

# data for rsms



## configuration of parameter keys

##### states at age for F random walk,
# use sms@catch.season.age for now
#sam_data$keyLogFsta
#sms@catch.season.age

x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nSpecies,dimnames=list(spNames,paste('age',ages)))
i<-1L
xx<-NULL
for (s in 1:nSpecies) {
  i<-1L
  la<-info[s,"last-age-likelihood"]
  aa<-sort(unique(c(sms@catch.season.age[[s]],la+1)))
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

nlogF<-unlist(lapply(sms@catch.season.age,length))
nlogFto<-cumsum(nlogF)  
nlogFfrom<- c(1,head(cumsum(nlogF),-1)+1 )
nlogFfromTo<-matrix(c(as.integer(nlogFfrom),as.integer(nlogFto)),ncol=2)
  
##### states at age for N random walk,
sam_data$keyVarLogN
sam_parameters$logSdLogN

nlogN<-as.integer(info[,"last-age"]-sms@first.age+1L)

keyVarLogN<-keyLogFsta  # no suitable SMS structure available
keyVarLogN[,]<- -1
keyVarLogN[,1:2]<-c(seq(from=1,by=2,length.out=length(nlogN)),seq(from=2,by=2,length.out=length(nlogN)))   
for (s in (1:nSpecies)) keyVarLogN[s,3:(info[s,"last-age"]+off.age)] <-keyVarLogN[s,2] 
# keyVarLogN

logSdLogN<-rep( c(0.35,-0.35),nSpecies)
# logSdLogN

#keyVarLogN<-lapply(1:nSpecies,function(s) y<-keyVarLogN[s,keyVarLogN[s,]>0])
 


nlogNto<-cumsum(nlogN)  
nlogNfrom<- c(1,head(cumsum(nlogN),-1)+1L )
nlogNfromTo<-matrix(c(as.integer(nlogNfrom),as.integer(nlogNto)),ncol=2)


#keyLogNsta[,1:2]<-c(seq(from=1,by=2,length.out=length(nlogN)),seq(from=2,by=2,length.out=length(nlogN)))   
#keyLogNsta

#### variance of catch and surveys
# sam_data$keyVarObs

#split into two groups: catches and surveys
# sms@catch.s2.group
x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nSpecies,dimnames=list(spNames,paste('age',ages)))
i<-1L
xx<-NULL
for (s in 1:nSpecies) {
  la<-info[s,"last-age-likelihood"]
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
  
# read survey data and options into FLR objects
indices<-SMS2FLIndices(control=sms,path=file.path(root,dir))
indices[[1]]@range.SMS

nFleets<-length(indices)
spNo<-unlist(lapply(indices,function(x) x@range.SMS["species"]))
PowerAge<-unlist(lapply(indices,function(x) x@range.SMS["power.age"]))
season<-unlist(lapply(indices,function(x) x@range.SMS["season"]))
PowerAge[PowerAge<0]<- -9
fleetNames<-unlist(lapply(indices,function(x) x@name))
q.age<-unlist(lapply(indices,function(x) x@range.SMS['q.age']))

a<-do.call(data.frame,lapply(indices,function(x) x@range))
a<-rbind(a,s=spNo,q.age=q.age,f=1:nFleets,PowerAge=PowerAge,q=season)
a["plusgroup",]<-0L
a['mina',]<-a['min',]+off.age
a['maxa',]<-a['max',]+off.age
a['miny',]<-a['minyear',]+off.year
a['maxy',]<-a['maxyear',]+off.year
ra<-rownames(a)
ra[1]<-'minage'
ra[2]<-'maxage'
rownames(a)<-ra

x<-matrix(0.0,ncol=1L,nrow=nFleets,dimnames=list(paste(1:nFleets,fleetNames),c("sampleTimeWithin")))
x[,"sampleTimeWithin"]<- (as.numeric(a['endf',])-as.numeric(a['startf',]))/2
sampleTimeWithin<-x

a<- a%>% mutate_if(is.numeric,as.integer)
a<-t(a)[,c('f','s','minyear',"miny","maxyear","maxy","minage","mina","maxage","maxa","q.age","plusgroup","PowerAge","q","startf","endf")]
if (is.vector(a)) {lena<-length(a); nama=names(a); a<-matrix(a,nrow=1,ncol=lena); colnames(a)<-nama}
rownames(a)<-paste(1:nFleets,fleetNames)
keySurvey<-a
keySurvey.df<-as.data.frame(a) %>% tibble::rownames_to_column("fName")

keySurvey

## survey variance
# sam_data$keyVarObs

x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nFleets,dimnames=list(fleetNames,paste('age',ages)))
i<-1L
xx<-NULL
for (f in (1:nFleets)) {
  jjj<-sort(unique(c(indices[[f]]@range.SMS$var.age.group,keySurvey[f,'maxage']+1)+off.age))
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

sam_parameters$logSdLogObs
## survey catchability
# sam_data$keyLogFpar

x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nFleets,dimnames=list(fleetNames,paste('age',ages)))
i<-0L
xx<-NULL
for (f in (1:nFleets)) {
  s<-keySurvey[f,'s']
  for (j in (keySurvey[f,'mina']:keySurvey[f,'maxa'])) {
     if (j<=keySurvey[f,'q.age']+off.age+1L) i<-i+1L
     x[f,j ]<-i
     xx<-rbind(xx,data.frame(f=f,species.n=s,age=j-off.age,a=j,keyCatchability=i))
  }
}

keyCatchability<-x
keyCatchability.df<-xx
logCatchability<-rep(0.0,max(keyCatchability)) 

## Catchability power age
sam_data$keyQpow
x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nFleets,dimnames=list(fleetNames,paste('age',ages)))
i<- -1L
xx<-NULL
for (f in (1:nFleets)) {
  s<-keySurvey[f,'s']
  for (j in (keySurvey[f,'mina']:keySurvey[f,'maxa'])) {
    if (j==keySurvey[f,'PowerAge']+off.age) i<-i+1L
    x[f,j ]<-i
    xx<-rbind(xx,data.frame(f=f,s=s,species.n=s,age=j-off.age,a=j,keyPowerQ=i))
  }
}
x
keyQpow<-x
keyQpow.df<-xx
if (max(keyQpow) >0) logQpow<-rep(0.0,max(keyQpow)) else  logQpow<-rep(0.0,0)

### survey observations

cpue<-lapply(indices,function(x){
  nage<-x@range["max"]-x@range["min"]+1L
  a<-as.data.frame(x@catch.n /  rep(as.vector(x@effort),each=nage))
  a$species.n<-unlist(x@range.SMS["species"])
  a$fleet.no<-unlist(x@range.SMS["fleet.no"])
  a$f<-unlist(x@range.SMS["sp.fl"])
  a
})

cpue<-do.call(rbind,cpue)
head(cpue)
str(cpue)


cpue<-cpue %>% mutate(y=year+off.year,q=as.integer(as.character(cpue$season)),a=age+off.age,s=species.n)%>%
  mutate(quarter=season,obs=data) %>%
  dplyr::select(year,y,season,q,age,a,species.n,s,fleet.no,f,obs) %>%
  filter(obs>0) %>% arrange(s,fleet.no,y,a) 
 cpue$obs.no<-1L:dim(cpue)[[1]]
head(cpue)

cpue<-cpue %>% filter(year<=sms@last.year.model)   # SKAL RETTES
fobs<-cpue %>% select(s,f,y,a) %>% group_by(s,f) %>%summarize(n=dplyr::n()) %>%  ungroup() %>% mutate(last=cumsum(n),first=lag(last)+1L)
fobs[1,"first"]<-1L

tail(fobs)
keySurvey.overview<-cbind(keySurvey,n=as.vector(fobs$n),first=as.vector(fobs$first),last=as.vector(fobs$last))
str(keySurvey.overview)

# keyVarObsSurvey.df
# keyCatchability.df
# keyQpow.df

k<-full_join(keyVarObsSurvey.df,keyCatchability.df) %>% dplyr::select(f, s, a, keyVarObsSurvey, keyCatchability) %>% arrange(f,s,a)
k<-full_join(k,select(keyQpow.df,s,f,a,keyPowerQ),by = join_by(f, s, a))

cpue<-left_join(cpue,k,by = join_by(a, s, f)) %>% arrange(f,s,a)
cpue$obs.no<-1:dim(cpue)[[1]]

keySurvey<-  cpue %>% dplyr::select(obs.no,f,s,y,a,q ,keyVarObsSurvey, keyCatchability,keyPowerQ)
stopifnot(dim(filter(keySurvey, is.na(keyVarObsSurvey) | is.na(keyCatchability)))[[1]]==0)
keySurvey<-keySurvey %>% mutate_if(is.numeric,as.integer)

keySurvey<-as.matrix(keySurvey)

logSurveyObs<-log(cpue$obs)


source(file.path(rsms,"from_sms_format_to_rsms_data.R"))

d<-From_SMS_format_to_rsms(otherPredExist=multi.environment,catchMultiplier=1,dir=file.path(root,dir))

# str(d,2)

#merge data
b<-full_join(d$catch,d$bio,join_by(year, species.n, quarter, sub_area, age))
if (multi) {
  b2<-full_join(d$mean_l,d$consum,join_by(year, species.n, quarter, sub_area, age))
  b<-full_join(b,b2,join_by(year, species.n, quarter, sub_area, age))
}
filter(b,age==0 & CATCHN>0)



b$seasFprop<-0.25
b[b$age==0,"seasFprop"]<-0
b[b$age==0 & b$q %in% c(3,4),"seasFprop"]<-0.5

b<-b %>% mutate(y=year+off.year,q=quarter,s=species.n+off.species,a=age+off.age)
b<-left_join(b,data.frame(s=info[,'s'],la=info[,"last-age"]+off.age),by = join_by(s)) %>% filter(a<=la)


propMat         <-by(b,b$s,function(x) {y<-tapply(x$PROPMAT,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
stockMeanWeight <-by(b,b$s,function(x) {y<-tapply(x$WSEA,list(x$y,x$q,x$a),sum) ; y[is.na(y)]<-0; y})
catchMeanWeight <-by(b,b$s,function(x) {y<-tapply(x$WCATCH,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
natMor          <-by(b,b$s,function(x) {y<-tapply(x$M,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
seasFprop       <-by(b,b$s,function(x) {y<-tapply(x$seasFprop,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})


#zero<-lapply(propMat,function(x) x[,,]<-0) virker ikke ?
zero<-propMat
for (s in 1:nSpecies) {zero[[s]][,,]<-0}


#### catch observations

catch<-aggregate(CATCHN~year+species.n+ sub_area+ age+y+ s+ a,FUN=sum,data=subset(b,CATCHN>=0 )) %>% filter(CATCHN>0)
catch<-left_join(catch,data.frame(s=info[,'s'],faf=info[,"first-age F>0"]+off.age),by = join_by(s)) %>% filter(a>=faf)

catch<-catch %>% arrange(s,y,a)
catch$obs.no<-1:dim(catch)[[1]]
head(catch)
logCatchObs<-log(catch$CATCHN)

keyCatch<-catch %>% mutate(CATCHN=NULL) %>% mutate_if(is.numeric,as.integer)
k<-full_join(keyLogFsta.df,keyVarObsCatch.df,by = join_by(s, age, a))
keyCatch<-left_join(keyCatch,k,by = join_by(s, age, a))
keyCatch<-keyCatch %>% arrange(obs.no)
#str(keyCatch,2)

stopifnot(dim(filter(keyCatch,is.na(keyLogFsta)))[[1]]==0)
keyCatch<-as.matrix(keyCatch)


info<-cbind(info,
    la=info[,"last-age"]+off.age,
    faf=info[,"first-age F>0"]+off.age,
    lasel=info[,"last-age-selec"]+off.age,
    lalike=info[,"last-age-likelihood"]+off.age)

dat<-list(                                          # Description of data, most of the text is copied from the ?sam.fit description
  sms.mode=sms@VPA.mode,                            #
  info=info,                                        # Various information for each species with analytical assessment
  nSpecies=nSpecies,                                # Number of species with analytical assessment
  #nOthSpecies=nOthSpecies,                          # Number of other predators (used for multispecies mode)   
  nYears=nYears,                                    # Number of years used in the model
  nAges=nAges,
  minAge=as.integer(sms@first.age),                 # A vector of integers defining the the lowest age class in the assessment for each species.
  maxAge=as.integer(sms@max.age.all),               # Maximum age for all species. The actual age range by species is given in "info"
  recAge=as.integer(sms@first.age)+off.age,          # recruitmet age (index)
  maxAgePlusGroup=as.integer(info[,'+group']),                  #  Is last age group considered a plus group (1 yes, or 0 no).
  years=sms@first.year.model:sms@last.year.model,   # A vector of the years used in the model
  spNames=spNames,                                  # Species names of species with analytical assessment
  nSeasons=nSeasons,                                # number of seasons
  spawnSeason=1L,
  recSeason=recSeason,                              # recruitment season
  fbarRange=fbarRange,                              # Minimum and maximum age by species used to calculate Fbar 
  off.age=off.age,                                  # offset between age and the age index used in RSMS (where fist age index  have to be 1 )
  off.year=off.year,                                # offset between year and the year index used in RSMS (where fist year index  have to be 1 )
  #off.season=off.season,
  #off.species=off.species,
  #off.oths=off.oths,
  stockRecruitmentModelCode=stockRecruitmentModelCode,
  keyLogFsta=keyLogFsta,                         # A matrix of integers. The number of rows is equal to the number of species fleets and the number of columns is equal to the number of age classes. The matrix describes the coupling of the fishing mortality states. '-1' is used for entries where no fishing mortality applies. For the valid entries consecutive integers starting 1 must be used, because they are used as indices in the corresponding state vector. If the same number is used for two fleet-age combinations, then the fishing mortality for those are assumed equal (linked to the same state).
  nlogF=nlogF,                                   # A vector with the sum of number of F state age groups by species. Each species number corresponds to the  number of unique (and not -1 value) in each (species) row of keyLogFsta
  #keyLogFsta.list=keyLogFsta.list,
   # keyLogSdLogFsta  keyLogFstaSdSp=keyLogFstaSdSp,                 # A matrix which links each F state age group by species to the logSdLogFsta=
  nlogFfromTo=nlogFfromTo,
  nlogN=nlogN,                                   # A vector with the sum of number of N state age groups by species
  nlogNfromTo=nlogNfromTo,
  
  #keyLogNsta=keyLogNsta,                         # ?????????????????species/age key for sd of N random 
  #keyLogNstaSdSp=keyLogNstaSdSp,                 #

  #keyLogFstaSdSp=keyLogFstaSdSp,
  keyVarLogN=keyVarLogN,
  keyVarObsCatch=keyVarObsCatch,
  keyVarObsSurvey=keyVarObsSurvey,
  keyCatchability=keyCatchability,
  propMat=  propMat,      
  stockMeanWeight=stockMeanWeight,
  catchMeanWeight=catchMeanWeight,
  seasFprop=seasFprop,
  natMor= natMor,
  propF=zero,
  propM=zero,
  keyCatch=keyCatch,
  logCatchObs=logCatchObs,
  #nCatchObs=length(logCatchObs),
  keySurvey=keySurvey,
  keySurvey.overview=keySurvey.overview,  # not really used
  logSurveyObs=logSurveyObs,
  sampleTimeWithinSurvey=sampleTimeWithin
  #nSurveyObs=length(logSurveyObs)
)


#str(sam_data)
str(dat,1)

#str(sam_parameters)

logSdLogFsta<-rep(-0.7,max(keyLogFsta))
rho<-rep(0.7,nSpecies)
Un<-matrix(0.0, nrow=sum(info[,"last-age"]-sms@first.age+1L),ncol=nYears )
Uf<-matrix(0.0, nrow=length(unlist(sms@catch.season.age)),ncol=nYears)

parameters<-list(
  logSdLogObsCatch=logSdLogObsCatch,
  logCatchability=logCatchability,
  logSdLogObsSurvey=logSdLogObsSurvey,
  logQpow=logQpow,
  logSdLogFsta=logSdLogFsta,
  logSdLogN=logSdLogN,
  rho=rho,
  rec_loga=rec_loga, 
  rec_logb=rec_logb, 
  Un=Un,
  Uf=Uf
)
data<-dat

save(data,parameters,file=file.path(rsms.root,"rsms_input.Rdata")) 



if (annual) {
  
  #load(file=file.path(rsms.root,"rsms_input.Rdata"),verbose=TRUE)
  
  data$nSeasons<-1L
  data$recSeason<-1L
  
  data$propMat <-lapply(data$propMat,function(x) x[,1,,drop=FALSE])
  data$stockMeanWeight<- lapply(data$stockMeanWeight,function(x) {
    x<-apply(x,c(1,3),mean) 
    x<-array2DF(x) %>% mutate(q=1) %>% mutate_if(is.character,as.integer)
    x<-tapply(x$Value,list(x$Var1,x$q,x$Var2),sum)
  })
  data$catchMeanWeight<- lapply(data$catchMeanWeight,function(x) x[,1,,drop=FALSE])  # to be changed
  data$seasFprop<- lapply(data$seasFprop,function(x) {
    x<-apply(x,c(1,3),sum) 
    x<-array2DF(x) %>% mutate(q=1) %>% mutate_if(is.character,as.integer)
    x<-tapply(x$Value,list(x$Var1,x$q,x$Var2),sum)
  })
  data$propF<- lapply(data$propF,function(x) x[,1,,drop=FALSE])   
  data$propM<- lapply(data$propM,function(x) x[,1,,drop=FALSE])
  data$natMor<- lapply(data$natMor,function(x) {
    x<-apply(x,c(1,3),sum) 
    x<-array2DF(x) %>% mutate(q=1) %>% mutate_if(is.character,as.integer)
    x<-tapply(x$Value,list(x$Var1,x$q,x$Var2),sum)
  })
  data$seasFprop<- lapply(data$seasFprop,function(x) {
    x<-apply(x,c(1,3),sum) 
    x<-array2DF(x) %>% mutate(q=1) %>% mutate_if(is.character,as.integer)
    x<-tapply(x$Value,list(x$Var1,x$q,x$Var2),sum)
  })
  
  x<-cbind(data$sampleTimeWithinSurvey,q=keySurvey.overview[,'q'])
  x[,"sampleTimeWithin"]<-x[,"sampleTimeWithin"]/4+(x[,"q"]-1)*0.25
  data$sampleTimeWithinSurvey<-x[,1]

  
  data$keySurvey[,'q']<-1L
  
  save(data,parameters,file=file.path(rsms.root,"rsms_input.Rdata"))
}
 list(data=data,parameters=parameters)
}


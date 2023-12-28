sms<-SMS.control  # just shorter name
multi<-!(sms@VPA.mode==0)

info<-sms@species.info
multi.environment<- any(info[,'predator']==2)

nSpecies<-sms@no.species-first.VPA+1L  # species with analytical assessment
nOthSpecies<-0L
ages<-sms@first.age:sms@max.age.all
years<-sms@first.year.model:sms@last.year.model; nYears<-length(years)
info<-info[first.VPA:nsp,c(1:3,5:6,9),drop=FALSE]
spNames<-dimnames(info)[[1]]
off.age<- as.integer(1-sms@first.age)
off.year<- - as.integer(sms@first.year.model-1)
if (multi.environment) off.species<as.integer(-first.VPA+1) else off.species<-0L
info<-cbind(info,s=(first.VPA:sms@no.species)+off.species)
off.oths<-0L #other species offset (for now)
off.season=0L # not used 
info<-cbind(info,s=(first.VPA:sms@no.species)+off.species)

load(file=file.path(sam.root,"sam_par_dat.Rdata"),verbose=TRUE)
str(sam_data)

# data for rsms
dat<-list(
  sms.mode=sms@VPA.mode,
  info=info,
  nSpecies=nSpecies,  # species with analytical assessment
  nOthSpecies=nOthSpecies,
  nYears=nYears,
  minAge=as.integer(sms@first.age),
  maxAge=as.integer(sms@max.age.all),
  years=sms@first.year.model:sms@last.year.model,
  spNames=spNames,
  off.age=off.age,
  off.year=off.year,
  off.season=off.season,
  off.species=off.species,
  off.oths=off.oths
)
str(dat,2)



## configuration of parameter keys

##### states at age for F random walk,
# use sms@catch.season.age for now
data$keyLogFsta
sms@catch.season.age
x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nSpecies,dimnames=list(spNames,paste('age',ages)))

i<-1L
xx<-NULL
for (s in 1:nSpecies) {
  la<-info[s,"last-age-likelihood"]
  aa<-sort(unique(c(sms@catch.season.age[[s]],la+1)))
  for (j in (1:(length(aa)-1))) {
    for (a in aa[j]:(aa[j+1]-1)) {
      x[s,a+off.age ]<-i;
      xx<-rbind(xx,data.frame(s=s,age=a,j=a+off.age,keyLogFsta=i))
    }
    i<-i+1L
  }
}

x
xx

keyLogFsta<-x
keyLogFsta.df<-xx

#### variance of catch and surveys
data$keyVarObs
#split into two groups: catches and surveys
sms@catch.s2.group
x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nSpecies,dimnames=list(spNames,paste('age',ages)))
i<-1L
xx<-NULL
for (s in 1:nSpecies) {
  la<-info[s,"last-age-likelihood"]
  aa<-sort(unique(c(sms@catch.s2.group[[s]],la+1)))
  for (j in (1:(length(aa)-1))) {
    for (a in aa[j]:(aa[j+1]-1)) {
      x[s,a+off.age ]<-i;
      xx<-rbind(xx,data.frame(s=s,age=a,j=a+off.age,keyVarObsCatch=i))
    }
    i<-i+1L
  }
}
keyVarObsCatch<-x
keyVarObsCatch.df<-xx

# read survey data and options into FLR objects
indices<-SMS2FLIndices(sms)

nFleets<-length(indices)
spNo<-unlist(lapply(indices,function(x) x@range.SMS["species"]))
fleetNames<-unlist(lapply(indices,function(x) x@name))
q.age<-unlist(lapply(indices,function(x) x@range.SMS['q.age']))
a<-do.call(data.frame,lapply(indices,function(x) x@range))
a

a<-rbind(a,s=spNo,q.age=q.age,f=1:nFleets)
a["plusgroup",]<-0L
a['minj',]<-a['min',]+off.age
a['maxj',]<-a['max',]+off.age
a['mini',]<-a['minyear',]+off.year
a['maxi',]<-a['maxyear',]+off.year
ra<-rownames(a)
ra[1]<-'minage'
ra[2]<-'maxage'
rownames(a)<-ra
a

x<-matrix(0.0,ncol=1L,nrow=nFleets,dimnames=list(paste(1:nFleets,fleetNames),c("sampleTimeWithin")))
x[,"sampleTimeWithin"]<- (as.numeric(a['endf',])-as.numeric(a['startf',]))/2
sampleTimeWithin<-x

a<- a%>% mutate_if(is.numeric,as.integer)
t(a)
a<-t(a)[,c(10,8,4,5,13,14,6,7,1,2,11,12,3,9)]

rownames(a)<-paste(1:nFleets,fleetNames)
keyFleet<-a
keyFleet.df<-as.data.frame(a) %>% tibble::rownames_to_column("fName")

keyFleet
keyFleet.df

## survey variance
data$keyVarObs

x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nFleets,dimnames=list(fleetNames,paste('age',ages)))
i<-1L
xx<-NULL
for (f in (1:nFleets)) {
  jjj<-sort(unique(c(indices[[f]]@range.SMS$var.age.group,keyFleet[f,'maxage']+1)+off.age))
  s<-keyFleet[s,'s']
  for (jj in (1:(length(jjj)-1))) {
    for (j in jjj[jj]:(jjj[jj+1]-1)) {
      x[f,j]<-i;
      xx<-rbind(xx,data.frame(f=f,s=s,j=j,age=j-off.age,keyVarObsSurvey=i))
    }
    i<-i+1L
  }
}

keyVarObsSurvey<-x
keyVarObsSurvey.df<-xx

## survey catchability
data$keyLogFpar
x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nFleets,dimnames=list(fleetNames,paste('age',ages)))
i<-0L
xx<-NULL
for (f in (1:nFleets)) {
  s<-keyFleet[s,'s']
  for (j in (keyFleet[f,'minj']:keyFleet[f,'maxj'])) {
     if (j<=keyFleet[f,'q.age']+off.age+1L) i<-i+1L
     x[f,j ]<-i
     xx<-rbind(xx,data.frame(f=f,species.n=s,age=j-off.age,j=j,keyCatchability=i))
  }
}
keyCatchability<-x
keyCatchability.df<-xx


source(file.path(rsms,"from_sms_format_to_rsms_data.R"))

d<-From_SMS_format_to_rsms(otherPredExist=multi.environment,catchMultiplier=1)

# str(d,2)
summary(d$catch)
summary(d$bio)
#merge data
b<-full_join(d$catch,d$bio,join_by(year, species.n, quarter, sub_area, age))
if (multi) {
  b2<-full_join(d$mean_l,d$consum,join_by(year, species.n, quarter, sub_area, age))
  b<-full_join(b,b2,join_by(year, species.n, quarter, sub_area, age))
}

b<-b %>% mutate(i=year+off.year,q=quarter,s=species.n+off.species,j=age+off.age)
head(b)

zero<-tapply(b$WSEA,list(b$s,b$i,b$q,b$j),sum)
zero[,,,]<-0.0;
dat<-list(
  propMat=        tapply(b$PROPMAT,list(b$s,b$i,b$q,b$j),sum),
  stockMeanWeight=tapply(b$WSEA,list(b$s,b$i,b$q,b$j),sum),
  catchMeanWeight=tapply(b$WCATCH,list(b$s,b$i,b$q,b$j),sum),
  natMor=         tapply(b$M,list(b$s,b$i,b$q,b$j),sum),
  propF=zero,
  propM=zero
)


#### catch observations

catch<-aggregate(CATCHN~year+species.n+ sub_area+ age+i+ s+ j,FUN=sum,data=subset(b,CATCHN>=0 )) %>% filter(CATCHN>0)
head(catch)

catch<-left_join(catch,data.frame(s=info[,'s'],faf=info[,"first-age F>0"]+off.age),by = join_by(s)) %>% filter(j>=faf)

catch<-catch %>% arrange(s,i,j)
catch$obs.no<-1:dim(catch)[[1]]
head(catch)
logCatchObs<-log(catch$CATCHN)

keyCatch<-catch %>% mutate(CATCHN=NULL) %>% mutate_if(is.numeric,as.integer)
k<-full_join(keyLogFsta.df,keyVarObsCatch.df,by = join_by(s, age, j))
keyCatch<-left_join(keyCatch,k,by = join_by(s, age, j))
keyCatch<-keyCatch %>% arrange(obs.no)
str(keyCatch,2)

stopifnot(dim(filter(keyCatch,is.na(keyLogFsta)))[[1]]==0)
keyCatch<-as.matrix(keyCatch)


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

cpue<-cpue %>% mutate(i=year+off.year,q=season,j=age+off.age,s=species.n)%>%
  dplyr::select(year,i,quarter=season,q,age,j,species.n,s,fleet.no,f,obs=data)%>%
  filter(obs>0) %>% arrange(s,fleet.no,i,j)
head(cpue)

# keyVarObsSurvey.df
# keyCatchability.df

k<-full_join(keyVarObsSurvey.df,keyCatchability.df) %>% dplyr::select(f, s, j, keyVarObsSurvey, keyCatchability)
cpue<-left_join(cpue,k,by = join_by(j, s, f)) %>% arrange(f,s,j)
cpue$obs.no<-1:dim(cpue)[[1]]

keySurvey<-  cpue %>% dplyr::select( obs.no,f,s,i,j,q,  fleet.no ,  keyVarObsSurvey, keyCatchability)
logSurveyObs<-log(cpue$obs)


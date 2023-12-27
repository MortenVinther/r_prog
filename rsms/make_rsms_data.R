multi<-SMS.control@VPA.mode==1
sms<-SMS.control  # just shorter versions
nSpecies<-sms@no.species-first.VPA+1L
ages<-sms@first.age:sms@max.age.all
years<-sms@first.year.model:sms@last.year.model
info<-sms@species.info[first.VPA:nsp,c(1:3,5:6,9),drop=FALSE]
spNames<-dimnames(info)[[1]]
off.age<- as.integer(1-sms@first.age)
off.year<- - as.integer(sms@first.year.model-1)

## configuration of parameters

##### states at age for F random walk,
# use sms@catch.season.age for now
data$keyLogFsta
sms@catch.season.age
x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nSpecies,dimnames=list(spNames,paste('age',ages)))

i<-1L
xx<-NULL
for (s in 1:nSpecies) {
  la<-info[s,"last-age-likelihood"]
  aa<-sort(unique(c(sms@catch.season.age[[s]],la)))
  for (j in (1:(length(aa)-1))) {
    for (a in aa[j]:(aa[j+1]-1)) {
      x[s,a+off.age ]<-i;
      xx<-rbind(xx,data.frame(species.n=s,age=a,j=a+off.age,keyLogFsta=i))
    }
    i<-i+1L
  }
  x[s,la+off.age ]<-i-1L;
  xx<-rbind(xx,data.frame(species.n=s,age=la,j=la+off.age,keyLogFsta=i-1L))
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
  aa<-c(sms@catch.s2.group[[s]],la)
  for (j in (1:(length(aa)-1))) {
    for (a in aa[j]:(aa[j+1]-1)) {
      x[s,a+off.age ]<-i;
      xx<-rbind(xx,data.frame(species.n=s,age=a,j=a+off.age,keyVarObsCatch=i))
    }
    i<-i+1L
  }
  x[s,la+off.age ]<-i-1L;
  xx<-rbind(xx,data.frame(species.n=s,age=la,j=la+off.age,keyVarObsCatch=i-1L))

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

a<-rbind(a,species.n=spNo,q.age=q.age,f=1:nFleets)
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

indices[[1]]@range.SMS
unlist(lapply(indices,function(x) x@range.SMS["species"]))

x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nFleets,dimnames=list(fleetNames,paste('age',ages)))
i<-1L
xx<-NULL
for (f in (1:nFleets)) {
  jjj<-sort(unique(c(indices[[f]]@range.SMS$var.age.group,keyFleet[s,'maxage'])+off.age))
  s<-keyFleet[s,'species.n']
  for (jj in (1:(length(jjj)-1))) {
    for (j in jjj[jj]:(jjj[jj+1]-1)) {
      x[f,j]<-i;
      xx<-rbind(xx,data.frame(f=f,s=s,j=j,age=j-off.age,keyVarObsSurvey=i))
    }
    i<-i+1L
  }
  x[f,j+1]<-i-1;
  xx<-rbind(xx,data.frame(f=f,s=s,j=j+1,age=j+1-off.age,keyVarObsSurvey=i-1))
}
x
xx
keyVarObsSurvey<-x
keyVarObsSurvey.df<-xx

## survey catchability
data$keyLogFpar
x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nFleets,dimnames=list(fleetNames,paste('age',ages)))
i<-0L
xx<-NULL
for (f in (1:nFleets)) {
  s<-keyFleet[s,'species.n']
  for (j in (keyFleet[f,'minj']:keyFleet[f,'maxj'])) {
     if (j<=keyFleet[f,'q.age']+off.age+1L) i<-i+1L
     x[f,j ]<-i
     xx<-rbind(xx,data.frame(f=f,species.n=s,age=j-off.age,j=j,keyCatchability=i))
  }
}
x
xx

keyCatchability<-x
keyCatchability.df<-xx

dat<-list(
  noYears=length(years),
  nSeasons=sms@last.season,
  nSpecies=nSpecies,
  minAge=sms@first.age,
  maxAge=sms@max.age.all,
  species.info<-sms@species.info[first.VPA:nsp,c(1:3,5:6,9)]
)


source(file.path(rsms,"from_sms_format_to_rsms_data.R"))

d<-From_SMS_format_to_rsms(otherPredExist=FALSE,catchMultiplier=1)
str(d,2)

ftable(xtabs(CATCHN~year+quarter+age,data=subset(d$catch,year==2001 & species.n==1)))

#merge data
b<-full_join(d$catch,d$bio,join_by(year, species.n, quarter, sub_area, age))
if (multi) {
  b2<-full_join(d$mean_l,d$consum,join_by(year, species.n, quarter, sub_area, age))
  b<-full_join(b,b2,join_by(year, species.n, quarter, sub_area, age))
}

b<-b %>% mutate(i=year+off.year,q=quarter,s=species.n,j=age+off.age)
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

head(b)
catch<-aggregate(CATCHN~year+species.n+ sub_area+ age+i+ s+ j,FUN=sum,data=subset(b,CATCHN>=0 )) %>% filter(CATCHN>0)
head(catch)

h<-split(catch,catch$species.n)
h<-lapply(h,function(x) filter(x,age>=info[x[1,'species.n'],"first-age F>0"]))
catch<-do.call(rbind,h)
catch<-catch %>% arrange(s,i,j)
catch$obs.no<-1:dim(catch)[[1]]
head(catch)
logCatchObs<-log(catch$CATCHN)
keyCatch<-catch %>% mutate(CATCHN=NULL) %>% mutate_if(is.numeric,as.integer)
k<-full_join(keyLogFsta.df,keyVarObsCatch.df,by = join_by(species.n, age, j))
keyCatch<-left_join(keyCatch,k,by = join_by(species.n, age, j))
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

keyVarObsSurvey.df

keyCatchability.df

k<-full_join(keyVarObsSurvey.df,keyCatchability.df) %>% dplyr::select(f, s, j, keyVarObsSurvey, keyCatchability)
k

str(cpue)
cpue<-left_join(cpue,k,by = join_by(j, s, f)) %>% arrange(f,s,j)
cpue$obs.no<-1:dim(cpue)[[1]]
head(cpue)
keySurvey<-  cpue %>% dplyr::select( obs.no,f,s,i,j,q,  fleet.no ,  keyVarObsSurvey, keyCatchability)
logSurveyObs<-log(cpue$obs)

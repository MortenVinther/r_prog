multi<-FALSE

sms<-SMS.control  # just shorter versions
nSpecies<-sms@no.species-first.VPA+1L
ages<-sms@first.age:sms@max.age.all
years<-sms@first.year.model:sms@last.year.model
info<-sms@species.info[first.VPA:nsp,c(1:3,5:6,9),drop=FALSE]
spNames<-dimnames(info)[[1]]
off.age<-1L-sms@first.age

## configuration of parameters

data$keyLogFsta
x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nSpecies,dimnames=list(spNames,paste('age',ages)))

xx<-NULL
i<-1
for (s in 1:nSpecies) {
  for (a in info[s,"first-age F>0"]:info[s,"last-age-selec"]) {
    x[s,a+off.age ]<-i;
    xx<-rbind(xx,data.frame(species=s,age=a,ai=a+off.age,keyLogFsta=i))
    i<-i+1L
  }
  for (a in info[s,"last-age-selec"]:info[s,"last-age-likelihood"]) {
     x[s,a+off.age ]<-i;
     xx<-rbind(xx,data.frame(species=s,age=a,ai=a+off.age,keyLogFsta=i))
  }

}
x
keyLogFsta<-x
keyLogFsta.df<-xx

## variance of catch and surveys
data$keyVarObs
#split into two groups: catches and surveys
sms@catch.s2.group
x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nSpecies,dimnames=list(spNames,paste('age',ages)))
i<-1
xx<-NULL
for (s in 1:nSpecies) {
  aa<-c(sms@catch.s2.group[[s]],info[s,"last-age-likelihood"])
  for (j in (1:(length(aa)-1))) {
    for (a in aa[j]:aa[j+1]) {
      x[s,a+off.age ]<-i;
      xx<-rbind(xx,data.frame(species=s,age=a,ai=a+off.age,keyVarObsCatch=i))
    }
    i<-i+1L
  }
}
keyVarObsCatch<-x
keyVarObsCatch.df<-xx

# read survey data and options into FLR objects
indices<-SMS2FLIndices(sms)
indices[[1]]@range.SMS
indices[[1]]@range
indices[[1]]@catch.n
indices[[1]]@effort
indices[[1]]@name

str(indices[[1]])

nFleets<-length(indices)

spNo<-unlist(lapply(indices,function(x) x@range.SMS["species"]))
fleetNames<-unlist(lapply(indices,function(x) x@name))

a<-do.call(data.frame,lapply(indices,function(x) x@range))
a
a<-rbind(a,sp=spNo)
a["plusgroup",]<-0L

x<-matrix(0.0,ncol=1L,nrow=nFleets,dimnames=list(paste(1:nFleets,fleetNames),c("sampleTimeWithin")))
x[,"sampleTimeWithin"]<- (as.numeric(a['endf',])-as.numeric(a['startf',]))/2
sampleTimeWithin<-x

a<- a%>% mutate_if(is.numeric,as.integer)
a<-t(a)[,c(1:5,8)]
str(a)
rownames(a)<-paste(1:nFleets,fleetNames)
keyFleet<-a

keyFleet.df<-as.data.frame(a) %>% tibble::rownames_to_column("fName")

#$ sampleTimesSeason        : num [1:3] 1 1
#$ sampleTimeWithin         : num [1:3] 00.5 0.5


## survey catchability
data$keyLogFpar



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

#merge data to get consistent age and yrar index
b<-full_join(d$catch,d$bio,join_by(year, species.n, quarter, sub_area, age))
if (multi) {
  b2<-full_join(d$mean_l,d$consum,join_by(year, species.n, quarter, sub_area, age))
  b<-full_join(b,b2,join_by(year, species.n, quarter, sub_area, age))
}

b <-b %>% mutate(species=unclass(factor(species.n)),
                      yi=unclass(factor(year)),
                      qi=unclass(factor(quarter)),
                      ai=unclass(factor(age)),
                   areai=unclass(factor(sub_area)))

str(b)

table_dat<-function(x)  x<-tapply(x,list(species,yi,qi,ai),sum)
zero<-table_dat(b$WSEA); zero[,,,]<-0.0;
dat<-list(
  propMat=table_dat(b$PROPMAT),
  stockMeanWeight=table_dat(b$WSEA),
  catchMeanWeight=table_dat(b$WCATCH),
  natMor=table_dat(b$M),
  propF=zero,
  propM=zero
)
str(data)




obs<-do.call(cbind,matrix(ff))
colnames(obs)<-c("year","fleet","age","species","obs")
str(obs)
head(obs)

keys<-obs[,1:4]
obs<-log(obs[,5])

dat<-list(dat,keys=keys,obs=obs)

str(dat)

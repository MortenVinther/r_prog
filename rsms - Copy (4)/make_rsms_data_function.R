make_rsms_data<-function(dir,sms.dat='sms.dat',seasFfrom=c('F','catch')[1],multi=FALSE) {
 
 # TEST single species and one species
#  dir=data.path; sms.dat='rsms.dat';seasFfrom=c('F','catch')[2];multi=FALSE
  
  
 #TEST    dir=data.path; sms.dat='rsms.dat';seasFfrom=c('F','catch')[2]; multi=TRUE


sms<-read.RSMS.control(dir,file=sms.dat)

info<-sms@species.info
oldAllSpeciesNames<-rownames(info)
oldPredNames<-rownames(info[info[,'predator']>=1,])
oldPredNo<-match(oldPredNames,rownames(info))

incl_other<- any(info[,'predator']==2)

if (incl_other) {
  off.species<-as.integer(-first.VPA+1) 
  off.oths<-as.integer(first.VPA+1)  #other species offset # in rsms, other species come after VPA (analytical)  species
  info<-rbind(info[first.VPA:sms@no.species,],info[1:(first.VPA-1),])
  nOthSpecies<-as.integer(first.VPA-1)
} else {
  off.species<-0L
  off.oths<- -9L
  info<-info[first.VPA:sms@no.species,,drop=FALSE]
  nOthSpecies<-0L
}

info<-cbind(s=1L:dim(info)[[1]],info)

off.season<-0L # not used 


nSpecies<-sms@no.species-first.VPA+1L  # species with analytically assessment
otherFoodn<-nSpecies+1L
nAllSp<-nSpecies+nOthSpecies
nPredator<-sum(info[,'predator']>0)

predators<-info[info[,'predator']>0,'s']
preys<-info[info[,'prey']>0,'s']
oldNewSp_n<-c((1:nOthSpecies)+nSpecies,1:nSpecies)

ages<-sms@first.age:sms@max.age.all
nAges<-length(ages)
recSeason<-as.integer(sms@rec.season)
nSeasons<-as.integer(sms@last.season)
years<-sms@first.year.model:sms@last.year.model; nYears<-length(years)

yq_idx<-matrix(1:(nYears*nSeasons),nrow=nYears,ncol=nSeasons,byrow=TRUE)
off.age<- as.integer(1-sms@first.age)
off.year<- - as.integer(sms@first.year.model-1)


fqa<-rep(1L,max(info[,'last-age'])+off.age); fqa[1]<-recSeason

info<-cbind(info,fModel=c(sms@fModel,rep(-9,nOthSpecies)))
fSepar<-if_else(sms@fModel==2,sms@firstAgeYearEffect+off.age,99L)
info<-cbind(info,fSepar=c(fSepar,rep(-9,nOthSpecies)))

spNames<-dimnames(info)[[1]][1:nSpecies]
if (nOthSpecies>0) {
  othspNames<-dimnames(info)[[1]][(nSpecies+1):nAllSp] 
  allSpNames<-c(spNames,othspNames)
} else {
  allSpNames<-spNames
  othspNames<-as.character(NULL)
}
allSpNamesLong<-sms@species.names.long[ match(allSpNames,oldAllSpeciesNames) ]

predNames<- dimnames(info[info[,'predator']>0,])[[1]]
preyNames<- dimnames(info[info[,'prey']>0,])[[1]]
otherFoodName<-'OTH'
stockRecruitmentModelCode<-as.integer(sms@SSB.R)
names(stockRecruitmentModelCode)<-spNames
use_rho<-as.integer(sms@use_rho)
#use_rho<-factor(use_rho, levels = c(0L,1L,2L),labels = c("no_correlation", "compound_symmetry", "AR(1)"))
rec_loga <- rep(1,nSpecies)
rec_loga[stockRecruitmentModelCode==1] <- log(200)  # Ricker
rec_loga[stockRecruitmentModelCode==2] <- 4         # B&W

rec_logb <- rep(1,nSpecies)
rec_logb[stockRecruitmentModelCode==1] <- -12  # Ricker
rec_logb[stockRecruitmentModelCode==2] <- -12        # B&W

fbarRange<-matrix(as.integer(sms@avg.F.ages),ncol=2,dimnames=dimnames(sms@avg.F.ages))+off.age;
colnames(fbarRange)<-c('first-a','last-a')

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

# not used  keyLogFstaSeason<-lapply(1:nSeasons,function(q) keyLogFsta)


nlogF<-unlist(lapply(sms@keyLogFsta,length))
nlogFto<-cumsum(nlogF)  
nlogFfrom<- c(1,head(cumsum(nlogF),-1)+1 )
nlogFfromTo<-matrix(c(as.integer(nlogFfrom),as.integer(nlogFto)),ncol=2)


keyLogFstaSd<-keyLogFsta
keyLogFstaSd[keyLogFstaSd<0]<- -9999L
keyLogFstaSd<-keyLogFstaSd+nlogFfrom-1
keyLogFstaSd[keyLogFstaSd<0]<- -1L
#keyLogFstaSd
logSdLogFsta<-rep(-0.7,max(keyLogFstaSd))


##### year and Ages in separable F-Model
# use sms@catch.sep.year and sms@catch.sep.age

x<-map2(sms@catch.sep.age,sms@catch.sep.year,function(x,y) {
  xx<-expand.grid(age=x,year=y)
})

for (s in 1:nSpecies) x[[s]]<-data.frame(x[[s]],s=s,fModel=info[s,'fModel'])

ay<-do.call(rbind,x) %>% filter(fModel==2) %>%arrange (s,age,year) 

if (dim(ay)[[1]]>0) ay<-  ay %>%  mutate(index=1:dplyr::n(), y=year-sms@first.year.model+1,a=age-sms@first.age+1)

x<-lapply(1:nSpecies,function(x) {
  matrix(-1L,ncol=info[x,'last-age']-sms@first.age+1L,nrow=nYears,dimnames=list(years,paste('age',sms@first.age:info[x,'last-age'])))
})

if (dim(ay)[[1]]>0) for (i in (1:dim(ay)[[1]])) {
  x[[ay[i,'s']]][ay[i,'y']:nYears,ay[i,'a']]<-ay[i,'index'] 
}

x<-lapply(x,function(xx){
  for (a in 2:ncol(xx)) if (xx[1,a]== -1) xx[,a]<- xx[,a-1]
  xx
})

keyLogSeparF<- x
if (dim(ay)[[1]]>0) logSeparF<-rep(0,max(ay$index)) else logSeparF<-numeric(0)

##### states at age for N ,

nlogN<-as.integer(info[1:nSpecies,"last-age"]-sms@first.age+off.age)

#### Process Noise N
# use sms@keyVarLogN
inclProcess<-sms@incl.process.noise; names(inclProcess)<-spNames
# inclProcess[c(2:3,8)]<-2 ; inclProcess[c(5,7)]<-3; inclProcess 
#    1= Both recruits and older ages have process noise estimated
#    2= Recruits only (N at age older than recruits in the first year become parameters)
#    3= Do not include process noise at all (N in the first year and recruits become parameters)

x<-matrix(-1L,ncol=sms@max.age.all-sms@first.age+1L,nrow=nSpecies,dimnames=list(spNames,paste('age',ages)))
i<-1L
xx<-NULL
for (s in 1:nSpecies) {
  la<-info[s,"last-age"]
  aa<-sort(unique(c(sms@keyVarLogN[[s]],la+1)))
  if (inclProcess[s]==2) aa<-aa[1] else if (inclProcess[s]==3)  aa<-numeric(0)
  if (inclProcess[s]==1) {
    for (j in (1:(length(aa)-1))) {
      for (age in aa[j]:(aa[j+1]-1)) {
        x[s,age+off.age ]<-i;
        # skal rettes xx<-rbind(xx,data.frame(s=s,age=age,a=age+off.age,keyVarLogN=i))
      }
      i<-i+1L
    }} else if (inclProcess[s]==2) {
      x[s,off.age ]<-i;
      i<-i+1L  
   }
}
#inclProcess;x
keyVarLogN<-x
# keyVarLogN.df<-xx

aa<-sort(unique(as.vector(keyVarLogN))) ;aa<-aa[aa>0]
logSdLogN<-rep(log(0.2),length(aa))  



nlogN[inclProcess==2]<-1L # recruits only
nlogN[inclProcess==3]<-0L

nlogNto<-cumsum(nlogN)  
nlogNfrom<- c(1,head(cumsum(nlogN),-1)+1L )
nlogNfromTo<-matrix(c(as.integer(nlogNfrom),as.integer(nlogNto)),ncol=2)
nlogNfromTo[nlogNfromTo[,1]>nlogNfromTo[,2]]<-c(0L,0L);
rownames(nlogNfromTo)<-spNames

#### N in first year or recruits parameter, if needed
nlogN2<-as.integer(info[1:nSpecies,"last-age"]-sms@first.age+off.age)
logNfirstYparam<- numeric(0)
logNrecruitParam<-numeric(0)
logNfirstYparamfromTo<-matrix(0L,nrow=nSpecies,ncol=2); rownames(logNfirstYparamfromTo)<-spNames
logNrecruitParamfromTo<-logNfirstYparamfromTo
iniN<-log(1E7)

for (s in 1:nSpecies) {
  if (inclProcess[s] %in% c(2,3))  {
     for (a in (2L:nlogN2[s])) logNfirstYparam<-c(logNfirstYparam,iniN-a)
     logNfirstYparamfromTo[s,1]<-max(logNfirstYparamfromTo)+1L
     logNfirstYparamfromTo[s,2]<-max(logNfirstYparamfromTo)+nlogN2[s]-2L
  }
  if (inclProcess[s] %in% c(3))  {
    logNrecruitParam<-c(logNrecruitParam,rep(iniN,nYears))
    logNrecruitParamfromTo[s,1]<-max(logNrecruitParamfromTo)+1
    logNrecruitParamfromTo[s,2]<-max(logNrecruitParamfromTo)+nYears-1

  }  
}
#inclProcess
#logNfirstYparamfromTo
#logNrecruitParamfromTo

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
minCVsurvey<-a$cv
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

logSdLogObsSurvey<-rep(log(0.4^2),max(keyVarObsSurvey))

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
    a<-a %>% mutate(obs=obs/mean(obs)) 
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

fobs<-cpue %>% dplyr::select(s,f,y,a) %>% group_by(s,f) %>%summarize(n=dplyr::n(),.groups = "drop") %>% mutate(last=cumsum(n),first=lag(last)+1L)
fobs[1,"first"]<-1L


keySurvey.overview<-cbind(keySurvey,n=as.vector(fobs$n),first=as.vector(fobs$first),last=as.vector(fobs$last))

sumCreep<-sum(keySurvey.overview[,'techCreep'])
logTechCreep<-rep(0.0,sumCreep)


#head( keyVarObsSurvey.df); filter(keyVarObsSurvey.df,is.na(keyVarObsSurvey))
#head(keyCatchability.df)
#head(keyQpow.df)

k<-full_join(keyVarObsSurvey.df,keyCatchability.df,by = join_by(f, s, a, age)) %>% dplyr::select(f, s, a, keyVarObsSurvey, keyCatchability) %>% arrange(f,s,a)
k<-full_join(k,dplyr::select(keyQpow.df,s,f,a,keyPowerQ),by = join_by(f, s, a))

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
#summary(filter(keySurvey, is.na(keyVarObsSurvey) | is.na(keyCatchability)))
filter(keySurvey, is.na(keyVarObsSurvey) | is.na(keyCatchability))

stopifnot(dim(filter(keySurvey, is.na(keyVarObsSurvey) | is.na(keyCatchability)))[[1]]==0)
keySurvey<-keySurvey %>% mutate_if(is.numeric,as.integer)

keySurvey<-as.matrix(keySurvey)

logSurveyObs<-log(cpue$obs)


x<-scan(file=file.path(dir,'recruitment_years.in'),comment.char='#',quiet=TRUE)
stopifnot(tail(x,1)== -999)
rec_years<-matrix(head(x,-1),nrow=nSpecies,byrow=TRUE)==1

source(file.path(rsms.root.prog,"from_sms_format_to_rsms_data.R"))


d<-From_SMS_format_to_rsms(otherPredExist=multi,catchMultiplier=1,dir=dir)


#merge data
b<-full_join(d$catch,d$bio,join_by(year, species.n, quarter, sub_area, age))
if (multi) {
  b2<-full_join(d$mean_l,d$consum,join_by(year, species.n, quarter, sub_area, age))
  d$other$WSEAoth<-d$other$WSEA; d$other$WSEA<-NULL; d$other$species.n.1<-NULL
  b2<-full_join(b2,d$other,join_by(year, species.n, quarter, sub_area, age))
  b<-full_join(b,b2,join_by(year, species.n, quarter, sub_area, age))
}

a<-data.frame(species.n=first.VPA:nsp,annualC=sms@combined.catches)
b<-left_join(b,a,by = join_by(species.n))
#head(b)

x<-b$age==fa & b$quarter < recSeason & b$annualC==0
b[x,'CATCHN']<-0
b[x,'WCATCH']<-0


b0<-filter(b,annualC==0) %>% group_by(year,species.n, sub_area, age) %>% mutate(seasFprop=CATCHN/sum(CATCHN,na.rm=T)) %>% ungroup() %>%
  dplyr::select(year, species.n, quarter, sub_area,   age, annualC, seasFprop)  
b1<-filter(b,annualC==1) %>%  dplyr::select(year, species.n, quarter, sub_area,age, annualC) %>% mutate(seasFprop=1/nSeasons) %>%
     mutate(seasFprop=if_else(age==fa & quarter >= recSeason,0.5, seasFprop)) %>%
     mutate(seasFprop=if_else(age==fa & quarter < recSeason,0.0, seasFprop))
#filter(b1,year==1975 & age==0)
bb<-rbind(b0,b1) 

check<-bb %>% group_by(year,species.n, sub_area, age) %>%summarize(Fprop=sum(seasFprop),.groups = "drop") 
# filter(check,!(Fprop>0.999 & Fprop<1.001) )
stopifnot(all(check$Fprop>0.999 & check$Fprop<1.001 | is.na(check$Fprop)))

b<-left_join(b,bb,by = join_by(year, species.n, quarter, sub_area, age,annualC))

b<-b %>% mutate(y=year+off.year,q=quarter,a=age+off.age)  

if (nSpecies>1 & incl_other) b<-b %>% mutate(s=oldNewSp_n[species.n]) else b<-b %>% mutate(s=species.n)  # change species numbering
b<-left_join(b,data.frame(s=info[,'s'],la=info[,"last-age"]+off.age),by = join_by(s)) %>% filter(a<=la)


#sort(unique(paste(b$species.n,b$s,sep='_')))

if (multi) {b[b$s>nSpecies,'WSEA']<-b[b$s>nSpecies,'WSEAoth']; b$WSEAoth<-NULL}

propMat         <-by(b,b$s,function(x) {y<-tapply(x$PROPMAT,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
stockMeanWeight <-by(b,b$s,function(x) {y<-tapply(x$WSEA,list(x$y,x$q,x$a),sum) ; y[is.na(y)]<-0; y})
catchMeanWeight <-by(b,b$s,function(x) {y<-tapply(x$WCATCH,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
natMor          <-by(b,b$s,function(x) {y<-tapply(x$M,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
catchNumber     <-by(b,b$s,function(x) {y<-tapply(x$CATCHN,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})

if (multi) {
  consum          <-by(b,b$s,function(x) {y<-tapply(x$CONSUM,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
  meanL           <-by(b,b$s,function(x) {y<-tapply(x$mean_l,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
  propM2          <-by(b,b$s,function(x) {y<-tapply(x$PROP_M2,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
  natMor1         <-by(b,b$s,function(x) {y<-tapply(x$M1,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
  otherN          <-by(b,b$s,function(x) {y<-tapply(x$N,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
} else {consum<-NA;meanL<-NA;propM2=NA; natMor1=NA;  otherN<-NULL}


trimIt<-function(x,used=c(1,2),spnames=spNames){
  x<-x[used]
  names(x)<-spnames
  return(x)
}

propMat         <-trimIt(x=propMat,used=1:nSpecies)
stockMeanWeight <-trimIt(x=stockMeanWeight,used=1:nSpecies)
catchMeanWeight <-trimIt(x=catchMeanWeight,used=1:nSpecies)
natMor          <-trimIt(x=natMor,used=1:nSpecies)
catchMeanWeight <-trimIt(x=catchMeanWeight,used=1:nSpecies)
catchNumber     <-trimIt(x=catchNumber,used=1:nSpecies)

if (multi) {
  consum          <-trimIt(x=consum,used=1:nAllSp,spnames=allSpNames)
  meanL           <-trimIt(x=meanL,used=1:nAllSp,spnames=allSpNames)
  propM2          <-trimIt(x=propM2,used=1:nSpecies)
  natMor1         <-trimIt(x=natMor1,used=1:nSpecies)
  otherN          <-trimIt(x=otherN,used=1:nAllSp,spnames=allSpNames)
  
  NN<-filter(b,!is.na(N)) %>% transmute(y,q,predNo=s,predAge=a,N) %>%filter(N>0)%>%mutate(N=NULL)
  
}
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
  zeroCatchYear<-bz %>%group_by(s,y) %>% summarize(z=sum(z),.groups = "drop") %>%  filter(z==0) %>% mutate(z=NULL)
} else zeroCatchYear<-NULL
zy<-rbind(zy,zeroCatchYear)

zeroCatchYear<-by(zy,list(zy$s),function(x) { x<-x$y[x$y>0]})
#lapply(zeroCatchYear,function(x) length(x))
zeroCatchYearExistsSp<-unlist(lapply(zeroCatchYear,function(x) length(x)))
zeroCatchYearExistsSp[zeroCatchYearExistsSp>0] <-1L

zy<-zy %>% mutate(zero=TRUE) %>% filter(y>0)
if (dim(zy)[[1]]>0) zeroCatchYearExists<-1L else zeroCatchYearExists<-0L



catch<-left_join(b,
                 data.frame(s=1:nAllSp,fa=sms@first.age,combCat=c(sms@combined.catches,rep(-9,nOthSpecies)),recSeason=sms@rec.season,faf=info[,"first-age F>0"]+off.age),
                 by = join_by(s)) %>%
       filter(combCat==1 | (combCat==0 & (age>fa | (age==fa & q>=recSeason )))) 
#round(ftable(xtabs(CATCHN~year+q+age,data=catch)))

catch<-left_join(catch,bz,by = join_by(s,y, q)) %>% mutate(CATCHN=if_else(z==1,CATCHN,0),z=NULL)
catch2<-catch

catch<-left_join(catch,zy,by = join_by(y, s)) %>% filter(is.na(zero)) %>% mutate(zero=NULL)

#  round(ftable(xtabs(CATCHN~year+q+age,data=filter(catch2,s==7))))

catch<-aggregate(CATCHN~year+species.n+ sub_area+ age+y+ s+ a,FUN=sum,data=filter(catch,CATCHN>0 & a>=faf)) %>% arrange(s,y,a)

#  round(ftable(xtabs(CATCHN~year+age,data=filter(catch,s==7))))


catch$obs.no<-1:dim(catch)[[1]]
logCatchObs<-log(catch$CATCHN)



keyCatch<-catch %>% mutate(CATCHN=NULL) %>% mutate_if(is.numeric,as.integer)
k<-full_join(keyLogFsta.df,keyVarObsCatch.df,by = join_by(s, age, a))
keyCatch<-left_join(keyCatch,k,by = join_by(s, age, a)) %>% mutate(species.n=NULL)
keyCatch<-keyCatch %>% arrange(obs.no)

stopifnot(dim(filter(keyCatch,is.na(keyLogFsta)))[[1]]==0)
keyCatch<-as.matrix(keyCatch)

b0<-filter(catch2,annualC==0) %>% group_by(year,species.n, sub_area, age) %>% mutate(seasFprop=CATCHN/sum(CATCHN,na.rm=T)) %>% ungroup() 
b1<-filter(catch2,annualC==1) %>%   mutate(seasFprop=0.25) %>%
  mutate(seasFprop=if_else(age==fa & quarter >= recSeason,0.5, seasFprop)) %>%
  mutate(seasFprop=if_else(age==fa & quarter < recSeason,0.0, seasFprop))
#filter(b1,year==1975 & age==0)
bb<-rbind(b0,b1) 

if (seasFfrom==c('F','catch')[2]) seasFprop <-by(bb,bb$s,function(x) {y<-tapply(x$seasFprop,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y}) else {
  seasFprop <-by(bb,bb$s,function(x) {y<-tapply(x$propF,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
}
ageQCatchMis<-lapply(seasFprop,function(x) apply(x,c(1,3),function(x){sum(x)<=0}))

logSeasFprop<-seasFprop
for (s in (1:nSpecies)) logSeasFprop[[s]][seasFprop[[s]]>0] <-log(logSeasFprop[[s]][seasFprop[[s]]>0])

# round(seasFprop[[1]][,,3],2) # cod age 3-1
# round(seasFprop[[6]][,,3],2) #her age 3-1
# round(seasFprop[[7]][,,3],2) #age 3 -1=2
#round(seasFprop[[7]][,,1],2) #age 1 -1=0
# 
# if (adj_san0 ) {
#   sa<-grep('sandeel',spNames)
#   if (length(sa) >0) {
#     if (nSeasons==2) for (i in sa)  seasFprop[[i]][,,1]<-rep(c(0,1),each=nYears)  # age 0
#     if (nSeasons==4) for (i in sa)  seasFprop[[i]][,,1]<-rep(c(0,0,1,0),each=nYears)  # age 0
#   }
# }
# # round(ftable(seasFprop[[1]][,,]),2)

catchNoSeason   <-by(b,b$s,function(x) {y<-tapply(x$CATCHN,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})
catchNoSeasonUsed   <-by(catch2,catch2$s,function(x) {y<-tapply(x$CATCHN,list(x$y,x$q,x$a),sum); y[is.na(y)]<-0; y})



info<-cbind(info,
    la=info[,"last-age"]+off.age,
    faf=info[,"first-age F>0"]+off.age)
 
###############   stomach data

if (multi) {
  
  oldNewPred <-function(opt,label,info1=info) {
    x<-  slot(sms,opt)
    names(x)<-oldPredNames
    x<-x[predNames]
    info1<-cbind(info1,xx= -9)
    info1[rownames(info1) %in% predNames,'xx']<-x
    colnames(info1)<-c(head(colnames(info1),-1),label)
    return(info1)
  }

  info<-oldNewPred(opt='size.selection',label='sizeSlct')
  info<-oldNewPred(opt='sum.stom.like',label='sumStomLike')
  info<-oldNewPred(opt='stom.max.sumP',label='stomMaxSumP')
  info<-oldNewPred(opt='max.stom.sampl',label='stomMaxSamp')
  info<-oldNewPred(opt='size.other.food.suit',label='sizeOtherFoodSuit')
  info<-oldNewPred(opt='stomach.variance',label='stomachVariance')
  
s<-read_delim(file.path(stom.input,"stomcon_list_Boots_0500_haul_as_observed.dat"),show_col_types = FALSE) # the used one

if (FALSE) {
  chk<-s %>% group_by(SMS_area,year,quarter,pred,pred.no,pred.size,pred.size.class) %>% summarize(tot=sum(stomcon),.groups = "drop") %>% filter(tot>1.05)
  summary(chk)
  sort(unique(paste (chk$pred.no,chk$pred)))
  
  chk<-filter(s,pred=='GNT') %>% arrange(year,quarter,pred,prey,prey.size)
  print(chk,n=100)
  chk  # it is OK , even though OTH seems much larger than the rest of preys (for FUL at least)
}

s[s$pred=='W_H','pred']<-'WHM'
s[s$pred=='N_H','pred']<-'NHM'

preds<-info[info[,'predator']>0,'s']
minStom<-data.frame(pred=names(preds),minStom=sms@min.stom.cont)
s<-left_join(s,minStom,by = join_by(pred)) %>% mutate(stomcon=if_else(minStom>0 & minStom>stomcon,minStom,stomcon)) %>% filter(minStom>0) %>% mutate(minStom=NULL)

s<-s%>%mutate(tpred=paste(formatC(pred.no,flag='0',w=2),pred),tprey=paste(formatC(prey.no,flag='0',w=2),prey),prey.mean.length.ALK=NULL,pred.no=as.integer(pred.no),prey.no=as.integer(prey.no))
#xtabs(~tpred+tprey,data=s)  #OK but old species order

#s<-s %>% mutate(pred.no=oldNewSp_n[pred.no])
s$pred.no<-oldNewSp_n[s$pred.no]
excl<-!(s$prey %in% c(preyNames,otherFoodName))
#summary(excl)
soth<-as.data.frame(head(filter(s,prey==otherFoodName),1))

s[excl,'prey']<-otherFoodName
s[excl,'prey.no']<-soth$prey.no
s[excl,'prey.size']<-soth$prey.size
s[excl,'prey.size.class']<-soth$prey.size.class
s[excl,'prey.mean.length']<-soth$prey.mean.length
s[excl,'mean.weight']<-soth$mean.weight
s[excl,'calc.prey.number']<-soth$calc.prey.number
s[excl,'used.prey.number']<-soth$used.prey.number
s[excl,'haul.prey.no']<-soth$haul.prey.no
#names(s)
#summary(excl)

#sort(unique(s$type))

a<-left_join(data.frame(pred=oldPredNames,typ=sms@stom.type.include), 
  rbind(
  data.frame(typ=1,type=c('obs')),
  data.frame(typ=2,type=c('obs','mid')),
  data.frame(typ=3,type=c('obs','mid','tail'))
 ),relationship = "many-to-many",by = join_by(typ))

s<-inner_join(a,s,by = join_by(pred, type))


s<-s %>% mutate(noSampl=if_else(haul.no>info[pred,'stomMaxSamp'],info[pred,'stomMaxSamp'],haul.no),
                    phi=if_else(haul.no>info[pred,'stomMaxSumP'],info[pred,'stomMaxSumP'],phi))


s<-s %>% group_by(SMS_area,year,quarter,pred,pred.no,pred.size,pred.size.class,pred.mean.length,noSampl,phi,
               prey,prey.no,prey.size,prey.size.class,prey.mean.length,type,mean.weight,haul.prey.no,calc.prey.number,used.prey.number) %>%
               summarize(stomcon=sum(stomcon),.groups = "drop")
                         
                         
x2<-filter(s,prey=='OTH')
s<-filter(s,prey!='OTH') %>% mutate(prey.no=oldNewSp_n[prey.no])
s<-rbind(s,x2)
s<-s%>%mutate(tpred=paste(formatC(pred.no,flag='0',w=2),pred),tprey=paste(formatC(prey.no,flag='0',w=2),prey))
# xtabs(~tpred+tprey,data=s) #OK


# vulnerability
pp<-matrix(0L,nrow=nSpecies+1,ncol=nSpecies+nOthSpecies,dimnames=list(c(spNames,'OTH'),c(spNames,othspNames)))
pp2<- s%>% dplyr::select(pred,pred.no,prey,prey.no) %>% mutate(prey.no=as.integer(ifelse(prey.no==0,nSpecies+1,prey.no))) %>%unique() 
#xtabs(~prey.no+pred.no,data=pp2)
for (i in (1:dim(pp2)[[1]]))  pp[unlist(pp2[i,'prey.no']),unlist(pp2[i,'pred.no'])]<-1L
predPrey<-pp
# predPrey

pp2$no<-1:dim(pp2)[[1]]

vulneraIdx<-predPrey
vulneraIdx[vulneraIdx>0]<-1:sum(vulneraIdx>0) 
# vulneraIdx
vulnera<-rep(0,max(vulneraIdx)) #log vulnerability coefficient, parameter 

vulneraIdx.df<-array2DF(vulneraIdx) %>% rename(prey=Var1,pred=Var2) %>% filter(Value>0)

# pred prey overlap
ov<-scan(file.path(dir,"season_overlap.in") ,comment.char = "#",quiet=TRUE)
ov<-head(ov,-1)
ov<-matrix(ov,ncol=nSpecies+1,nrow=nAllSp*4,dimnames=list(paste0(rep(c(spNames,othspNames),each=4),' Q',1:4), c(spNames,otherFoodName)),byrow=TRUE )
ov<-array(ov,dim=c(4,nAllSp,nSpecies+1),dimnames=list(paste0('Q',1:4),c(spNames,othspNames), c(spNames,otherFoodName)))                                                        
ovd<-array2DF(ov)
names1<-c( spNames,othspNames)
names2<-c(spNames,otherFoodName)
ovdt<-ovd%>%mutate(q=as.integer(substr(Var1,2,2)),predNo=match(Var2,names1),preyNo=match(Var3,names2)) %>% 
   transmute(q,predNo,preyNo,overlap=Value)
overlap<-ov
# ftable(overlap)
overlap[overlap>0] =log(overlap[overlap>0])
overlapIdx<-filter(ovdt,overlap<0) %>% unique() %>% mutate(overlap= -overlap)
overlapP<-rep(0,length(unique(overlapIdx$overlap))) # log overlap parameter

# put weight on from length
lw<-scan(file.path(dir,"length_weight_relations.in") ,comment.char = "#",quiet=TRUE)
lw<-head(lw,-1)
lw<-matrix(lw,ncol=2,nrow=nSpecies+nOthSpecies,dimnames=list(c(spNames,othspNames),c('a','b')),byrow=TRUE )
lwdf<-data.frame(lw) %>% mutate(sp=1L:dim(lw)[[1]])

stomRange<-data.frame(pred.no=oldPredNo,srL=sms@size.range.lower/100,srU=sms@size.range.upper/100)
stomRange$pred.no<-oldNewSp_n[stomRange$pred.no]


s<-left_join(s,stomRange,by = join_by(pred.no))


#observed predator prey size ratio
ss<-s%>% left_join(x=s,y=lwdf,by = join_by(pred.no == sp)) %>% 
  mutate(pred.size.w=a*pred.mean.length**b,a=NULL,b=NULL,prey.size.w=mean.weight,pred.prey.size.w=pred.size.w/prey.size.w)  

a<-filter(ss,type=='obs' & prey !='OTH' ) %>% 
  mutate(log.pred.prey.size.w=log(pred.size.w/prey.size.w),log.pred.size.w=log(pred.size.w)) %>%
  select(SMS_area,pred.no, pred,tpred,prey.no,prey,tprey,log.pred.size.w,log.pred.prey.size.w,srL, srU ) 

#summary(exp(a$log.pred.size.w))

# constraint uniform size selection
# first the predator size independent prey range 

#aa<-a %>% group_by(pred.no,pred,tpred,prey.no,prey,tprey) %>% summarize(min_s=min(log.pred.prey.size.w),max_s=max(log.pred.prey.size.w),.groups = "drop") 

k<-a %>% group_by(pred.no,pred,prey.no,prey) %>% 
  summarize(min_s=quantile(log.pred.prey.size.w,probs=srL[1]),max_s=quantile(log.pred.prey.size.w,probs=srU[1]),n=dplyr::n(),.groups = "drop") 

if (FALSE) { # show the effects of quantile percentage
  t1<-a %>% group_by(pred.no,pred,prey.no,prey) %>% summarize(min_s=min(log.pred.prey.size.w),max_s=max(log.pred.prey.size.w),.groups = "drop") 
  t2<-a %>% group_by(pred.no,pred,prey.no,prey,srL,srU) %>% 
    summarize(min_sq=quantile(log.pred.prey.size.w,probs=srL[1]),max_sq=quantile(log.pred.prey.size.w,probs=srU[1]),n=dplyr::n(),.groups = "drop") 
  print(left_join(t1,t2,by = join_by(pred.no, pred, prey.no, prey)),n=200)
 
}
# round(xtabs(min_s~tpred+tprey,data=aa),3) 
# 
#round(xtabs(max_s~tpred+tprey,data=aa),3) 
#  round(xtabs(exp(min_s)~tpred+tprey,data=aa),2) ;round(xtabs(exp(max_s)~tpred+tprey,data=aa),2)

if (FALSE) { #quantile regression, not used
  pp<-array(0,dim=c(4,nSpecies,nSpecies+nOthSpecies),dimnames=list(c('lower_intc','lower_slope','upper_intc','upper_slope'),spNames,c(spNames,othspNames)) )                                                         
  for (i in (1:dim(aa)[[1]])) {
    pred<-unlist(aa[i,'pred'])
    prey<-unlist(aa[i,'prey'])
    pp[1,prey,pred]<-unlist(aa[i,'min_s'])
    pp[3,prey,pred]<-unlist(aa[i,'max_s'])
  }
  #  round(ftable(pp),2)
  
  
  
  # add parameters from quantile regressions
  
  library(quantreg)
  q1<-0.025 #lower quantile
  q2<-0.975 #higher quantile
  
  if(FALSE) by(a,list(a$pred.no),function(x) {
    
    pp<-ggplot(data=x, aes(x=log.pred.size.w, y=log.pred.prey.size.w)) +
      geom_point()+
      labs(title=x[1,'pred'],y='log(predator size/prey size)',x='log(predator size)')+
      geom_smooth(method = "lm", se =F)+
      facet_wrap(vars(prey),ncol=3)
  })
    
  
  if (FALSE) {
    a1<-filter(a,pred.no %in% c(1:9) )
    aa<-by(a1,list(a1$prey.no,a1$pred.no),function(x,minNobs=5) {
      titl<-paste(x[1,'pred'],'eating',x[1,'prey'])
      print(titl)
      
      if (dim(x)[[1]]<minNobs) {
        print('too few observations')
      } else {
        ru<-rq(log.pred.prey.size.w~log.pred.size.w ,data=x, tau=q1)
        ru<-ru[["coefficients"]]
        rl<-rq(log.pred.prey.size.w~log.pred.size.w ,data=x, tau=q2)
        rl<-rl[["coefficients"]]
        
        
        pp<-ggplot(data=x, aes(x=log.pred.size.w, y=log.pred.prey.size.w)) +
          geom_point()+
          labs(title=titl,y='log(predator size/prey size)',x='log(predator size)')+
          geom_smooth(method = "lm", se =F)+
          geom_abline(slope=ru[2],intercept=ru[1])+
          geom_abline(slope=rl[2],intercept=rl[1])
        print(pp)
      } 
      if (dim(x)[[1]]<minNobs) {
        return(list(pred=unlist(x[1,'pred']),prey=unlist(x[1,'prey'])))
      } else return(list(pred=unlist(x[1,'pred']),prey=unlist(x[1,'prey']),ru=ru,rl=rl,plt=pp))
      
    })
  }
  
  # sort(unique(a$tpred))
  a1<-filter(a,pred.no %in% c(1:5) )
  n<-a1 %>% group_by(pred.no,prey.no) %>% summarize(n=dplyr::n(),.groups = "drop")
  an<-left_join(a1,n,by = join_by(pred.no, prey.no))
  
  aa<-filter(an,n>5) %>% 
    nest_by(pred.no,pred,prey.no,prey) %>%  
    mutate(rl=list(rq(log.pred.prey.size.w~log.pred.size.w, tau=q1,data=data)),
           ru=list(rq(log.pred.prey.size.w~log.pred.size.w, tau=q2,data=data)))
  
  
  for (i in (1:dim(aa)[[1]])) {
    pred<-unlist(aa[i,'pred'])
    prey<-unlist(aa[i,'prey'])
    cof<-aa[i,'rl'][[1]][[1]][['coefficients']]
    pp[1,prey,pred]<-cof[1]
    pp[2,prey,pred]<-cof[2]
    cof<-aa[i,'ru'][[1]][[1]][['coefficients']]
    pp[3,prey,pred]<-cof[1]
    pp[4,prey,pred]<-cof[2]
  }  
  #round(ftable(pp),2)  
  
  predPreySizeParm<-pp
  
  names1<-c( spNames,othspNames)
  names2<-c(spNames,otherFoodName)
  
  k<-array2DF(pp) %>% transmute(preyNo=match(Var2,names2),predNo=match(Var3,names1),Var1,Value) %>% 
    pivot_wider(names_from = Var1, values_from = Value) 
  
  
  if (FALSE){
    u<-scan(file.path(stom.input,"pred_prey_size_range_param.csv"),comment.char = "#",sep=',') 
    u<-u[!is.na(u)]
    u<-head(u,-1)
    u<-matrix(u,nrow=(nSpecies)*4,ncol=nSpecies++nOthSpecies,byrow=TRUE)
    round(ftable(u),1)
    
    uu<-array(0,dim=c(4,nSpecies,nSpecies+nOthSpecies),dimnames=list(c('lower_intc','lower_slope','upper_intc','upper_slope'),spNames,c(spNames,othspNames)) )                                                         
    round(ftable(uu),1)
    for (i in (1:4)) {
      ul<-nSpecies*(i-1)+1
      up<-nSpecies*(i) 
      uu[i,,]<-u[ul:up,]
    }
    
    pp<-uu
    round(ftable(pp),1)
  }
}  # end not used qualtile regression


#sort(unique(pp2$prey))
#sort(unique(pp2$prey.no))
# the full set of possible combinations of pred, pred age, prey and prey age
pp3<-by(pp2,list(pp2$no),function(x){
  predAge<-1L:info[x$pred.no,'la']
  preyAge<-1L:info[x$prey.no,'la']
  expand.grid(y=1:nYears,q=1L:nSeasons,predNo=x$pred.no,predAge=predAge,preyNo=x$prey.no,preyAge=preyAge)
})
pp3<-do.call(rbind,pp3)
dim(pp3)

#filter(pp3,preyAge==1 & preyNo==1 & q==3 & y==1)
#xtabs(~preyNo+preyAge+q,data=pp3)

pp3<-filter(pp3,(predAge > 1) | (q>=recSeason & predAge==1)) %>% 
     filter((    preyAge>1)   | (q>=recSeason &  preyAge==1)) %>% 
  filter(preyNo<otherFoodn | (preyNo==otherFoodn & preyAge==2))  %>%
  arrange(y,q,preyNo,preyAge,predNo,predAge)

dim(pp3)
#filter(pp3,preyAge==1 & preyNo==1 & q==3 & y==1)

#xtabs(~preyNo+preyAge+q,data=pp3)

#sort(unique(pp3$preyNo))
# all combinations
b7<-left_join(b,lwdf%>%transmute(s=sp,l_a=a,l_b=b),by = join_by(s))

aa2<-left_join(pp3%>% transmute(y,q,predNo,predAge,preyNo,preyAge),b7%>%  transmute(y,q,predNo=s,predAge=a,predW=WSEA,predL=mean_l,predWl=l_a*predL**l_b,N),by = join_by(y,q, predNo, predAge))
aa2<-left_join(aa2,b7%>%  transmute(y,q,preyNo=s,preyAge=a,preyW=WSEA, preyL=mean_l,preyWl=l_a*preyL**l_b),by = join_by(y,q, preyNo, preyAge))
# sort(unique(aa2$preyNo))
#xtabs(~preyNo+preyAge+q,data=aa2)

dim(aa2)
aa2<- aa2 %>% filter(predNo<=nSpecies | predNo>nSpecies & N>0)  %>% mutate(N=NULL)

# filter(aa2,preyNo==13)
aa2<- aa2%>%mutate(preyW=if_else(preyNo==otherFoodn,1,preyW))


aa2<-filter(aa2,predW>0 & preyW>0 & preyNo <=nSpecies+1)
#filter(aa2,preyNo==13)

# actual comb
aa3<-left_join(aa2,k,by = join_by(predNo==pred.no, preyNo==prey.no))
#filter(aa3,preyAge==1 & preyNo==1 & q==1 & y==1)
#filter(aa3,preyAge==1 & preyNo==1 & q==3 & y==1)

#aa3<-aa3 %>% mutate(predW=predWl,preyW=preyWl)
aaa<-aa3 %>%  mutate(logRatio=log(predW/preyW),predW=log(predW),preyW=log(preyW),predNo=as.integer(predNo),preyNo=as.integer(preyNo)) %>% 
       mutate(incl=if_else(logRatio >=min_s & logRatio <=max_s,TRUE,FALSE )) %>%
       mutate(incl=if_else(preyNo==otherFoodn,TRUE,incl)) %>% mutate( min_s=NULL,max_s=NULL,n=NULL)
xtabs(~preyNo+preyAge+q,data=aaa)
#filter(aaa,preyAge==1 & preyNo==1 & q==3 & y==1)
#summary(aaa)
#filter(aaa,preyNo==otherFoodn)
#ftable(xtabs(~incl+predNo+preyNo,data=aaa))

if (FALSE) ggplot(subset(aaa,predNo<=4), aes(x=predW, y=logRatio,col=factor(incl))) +
  geom_point() + 
    facet_wrap(vars(predNo),ncol=5,nrow=5)+
  labs(y='log(predator size/prey size)',x='log(predator size)')

if (FALSE) ggplot(subset(aaa,predNo==1 & preyNo<=6), aes(x=predW, y=logRatio,col=factor(incl))) +
  geom_point() + 
  facet_wrap(vars(predNo,preyNo),ncol=5,nrow=5)+
  labs(y='log(predator size/prey size)',x='log(predator size)')

if (FALSE) {
  t1<-xtabs(~preyNo+preyAge+q,data=aaa)
  t2<-xtabs(~preyNo+preyAge+q,data=filter(aaa,incl))
  t3<-xtabs(~preyNo+preyAge+q,data=filter(aaa,!incl))
  t2 
  q<-3
  cat('all\n');t1[,,q];cat('included\n');t2[,,q];cat('excluded\n');t3[,,q]
  # filter(aaa,preyNo==3 & q==3 &y==1)
}


#subset(aaa,predNo==23)
aaa<-aaa %>% filter(incl) %>% mutate(lower_intc=NULL,lower_slope=NULL, upper_intc=NULL, upper_slope=NULL,incl=NULL) %>%
        arrange(y,q,preyNo,preyAge,predNo,predAge)

#filter(aaa,preyNo==otherFoodn)
#dim(aaa)
#subset(aaa,predNo==13 )
#head(subset(aaa,predNo==1 & predAge<=3 ))
#aaa$n<-1L:dim(aaa)[[1]]

#sort(unique(aaa$preyNo))
#add vulnerability param index
aaa<-left_join(aaa, pp2 %>% transmute(area=1L,predNo=pred.no,preyNo=prey.no,vulneraIdx=no),by = join_by(predNo, preyNo))

#filter(aaa,preyNo==otherFoodn)
#filter(aaa,predNo==1 & predAge==2 &y==1 & q==3)

if (FALSE) {
  one<-aaa %>% group_by(y,q,predNo,predAge) %>% summarize(n=dplyr::n(),.groups = "drop") # %>% ungroup()
  filter(one,n<=1)
  othErr<-left_join(one,aaa,by = join_by(y, q, predNo, predAge))
  sort(unique(othErr$predNo))
  filter(othErr,preyNo != otherFoodn)
  aaa<-left_join(aaa,one,by = join_by(y, q, predNo, predAge)) %>% filter(!(n==1 & preyNo==otherFoodn)) %>% mutate(n=NULL)
}
#summary(aaa); dim(a)


check2<-function(sp=1,yy=1,qq=1) {
  cod<-filter(aaa,predNo==sp & y==yy & q==qq)
  ftable(xtabs(~predAge+preyNo+preyAge,data=cod))
}

#check2(1)
#check2(2)
names(aaa)
filter(aaa, preyNo==1 & q==3 & y==1 &predNo==1)
a4<-aaa %>% select(y,q,predNo,predAge,predW,preyNo,preyAge,preyW,logRatio, vulneraIdx,preyNo,preyAge,preyW,logRatio,vulneraIdx) %>% mutate(area=1L)
filter(a4,preyAge==1 & preyNo==1 & q==3 & y==1)


suitIdx<-a4 %>% nest(data=c(predNo,predAge,predW,preyNo,preyAge,preyW,logRatio, vulneraIdx)) %>% 
  rowwise() %>% mutate(data=list(data %>% nest(data=c(preyNo,preyAge,preyW,logRatio,vulneraIdx))))


partM2Idx<-a4 %>% mutate(partM2=0,suit=0, availFood=0,predW=NULL,preyW=NULL,logRatio=NULL, vulneraIdx=NULL) %>%
  nest(data=c(predNo,predAge,preyNo,preyAge,partM2,suit,availFood)) %>% 
  rowwise() %>% mutate(data=list(data %>% nest(data=c(preyNo,preyAge,partM2,suit,availFood))))

#partM2Idx
#partM2Idx$data[[1]]
#partM2Idx$data[[1]]$data[[1]]


# suitIdx
# suitIdx[,'y']
# suitIdx[1,'data'][[1]][[1]]
# #suitIdx[1,'data'][[1]][[1]][['data']]
# suitIdx[1,'data'][[1]][[1]][['data']][[2]]



## stomach contents

oth<-ss$prey==otherFoodName

ss[oth,'prey.no']<-otherFoodn
ss[oth,'prey.size.class']<-1
ss[oth,'prey.size.w']<-1
ss[oth,'pred.prey.size.w']<-1 # reset later on

# glimpse(filter(ss,oth))
aBasis<-ss%>% transmute(area=as.integer(SMS_area),
                        year=as.integer(year),
                        y=as.integer(year+off.year),
                        q=as.integer(quarter+off.season),
                        predC=pred,
                        pred=pred.no,
                        predSize=pred.size,
                        predSizeClass=as.integer(pred.size.class),
                        predMeanLength=as.integer(pred.mean.length),  
                        predSizeW=pred.size.w,
                        noSampl=as.integer(noSampl),
                        phi=phi,
                        preyC=prey,
                        prey=as.integer(prey.no),
                        preySize=prey.size,
                        preySizeClass=as.integer(prey.size.class),
                        preyMeanLength=as.integer(prey.mean.length),
                        preySizeW=prey.size.w,
                        logPPsize=log(pred.prey.size.w),
                        type=type,
                        stomcon=stomcon )
aBasis[aBasis$preyC=='OTH','logPPsize']<-NA

# filter(aBasis,y==8,q==1,pred==1,predSizeClass==6) %>% arrange(prey,preySize) %>% select(-preySizeClass,-preyMeanLength,-predMeanLength)
# sum up to 1 (into proportions)
aBasis<-aBasis %>% group_by(area,year,y,q,predC,pred,predSize,predSizeClass,predMeanLength,predSizeW,noSampl,
                       phi,preyC,prey,preySize,preySizeClass,preyMeanLength,preySizeW,logPPsize,type) %>% summarize(stomcon=sum(stomcon),.groups = "drop") %>%  
  group_by(area,year,y,q,predC,pred,predSize,predSizeClass,predMeanLength,predSizeW,noSampl,phi) %>% mutate(stomcon=stomcon/sum(stomcon)) %>% ungroup()
   
# filter(aBasis,y==8,q==1,pred==1,predSizeClass==6) %>% arrange(prey,preySize) %>% select(-preySizeClass,-preyMeanLength,-predMeanLength)



#vulneraIdx.df
aBasis<-left_join(aBasis,select(vulneraIdx.df,pred,prey,vulneraIdx=Value),by=join_by(predC==pred,preyC==prey))

aBasis<-aBasis %>% group_by(area,year,y,q, pred,  predSizeClass,prey) %>% 
  mutate(stomconTot=sum(stomcon),minL=min(preySizeClass),nPreyGroups=dplyr::n(),firstL=if_else(minL==preySizeClass,TRUE,FALSE)) %>%
  ungroup() %>% arrange(area,year,y,q, pred,  predSizeClass,prey,preySizeClass) %>%
  group_by(area,year,y,q, pred,  predSizeClass) %>% mutate(nPreyGroups=sum(firstL),preyIdx=cumsum(firstL), minL=NULL) %>% ungroup()

# filter(aBasis,y==8,q==1,pred==1,predSizeClass==6) %>% arrange(prey,preySize)%>% select(-preySizeClass,-preyMeanLength,-predMeanLength)


#tmp<-filter(aBasis2,y==8&q==1&pred==1 & predSizeClass >5) %>% select(predSizeClass,prey,preySizeClass,stomcon,stomconTot,nPreyGroups,firstL,preyIdx) %>% arrange (predSizeClass , prey ,preySizeClass )
#print(tmp,n=20)


logStomObsVar<-1L:length(predNames)
names(logStomObsVar)<-predNames 
aBasis$stomObsVarIdx<-logStomObsVar[aBasis$predC]


#sort(unique(paste(aBasis$predC,aBasis$stomObsVarIdx)))
logStomObsVar[]<-log(0.5^2)

usestomObsVar<- !info[, 'stomachVariance'][predNames] %in% c(4)
#glimpse(aBasis)

# filter(aBasis,y==8,q==1,pred==1,predSizeClass==6) %>% arrange(prey,preySize)%>% select(-preySizeClass,-preyMeanLength,-predMeanLength)


if (FALSE) {
  chk<- aBasis %>% select(area, y,q,pred, predSizeClass, predSizeW ,noSampl, phi,nPreyGroups,prey,preySizeClass, preySizeW, logPPsize, type, firstL, stomconTot ,stomcon,preyIdx,vulneraIdx,stomObsVarIdx) %>% 
    nest(data=c(pred,predSizeClass,predSizeW,stomObsVarIdx,noSampl,   phi,nPreyGroups,prey,preySizeClass,preySizeW,logPPsize,type,stomcon, stomconTot ,preyIdx,firstL,vulneraIdx))
  chk[1,'data'][[1]]
  
  chk2<-  chk%>%rowwise() %>% mutate(data=list(data %>% nest(data=c(prey,preySizeClass,preySizeW,logPPsize,type,firstL,preyIdx,stomconTot ,stomcon,vulneraIdx))))
  chk2[1,'data'][[1]]
  chk2[1,'data'][[1]][[1]]$data
}


###### stom ALK
s<-read_delim(file.path(stom.input,"ALK_stom_list_Simple_0001_sample_id_as_observed.dat"),show_col_types = FALSE)
# sort(unique(s$prey))
d2010<-filter(s,year==2013) %>%mutate(year=2010)  #snyd
s<- rbind(filter(s,year !=2010),d2010)
s<-s%>%mutate(prey.no=as.integer(prey.no),prey.age=as.integer(prey.age),quarter=as.integer(quarter),year=as.integer(year),prey.size.class=as.integer(prey.size.class))


# new number order
s$prey.no<-oldNewSp_n[s$prey.no]
s<-filter(s,prey %in% union(preyNames,setdiff(predNames,othspNames)))  #just preys and analytical predators


alkBasis<-s%>% transmute(area=as.integer(SMS_area),
                        year,
                        y=as.integer(year+off.year),
                        q=as.integer(quarter+off.season),
                        sp=prey,
                        s=prey.no,
                        age=as.integer(prey.age),
                        a=as.integer(prey.age+off.age),
                        size=prey.size,
                        sizeClass=as.integer(prey.size.class),
                        meanLength=ALK.MeanL,  
                        alk=ALK/100
)

#a        <-by(alkBasis,list(alkBasis$y,alkBasis$q),function(x) {y<-tapply(x$alk,list(x$s,x$a,x$sizeClass),sum); y[is.na(y)]<-0; y})
#str(a)
#round(ftable(a[[1]]),3)

alkBasis<-left_join(alkBasis,data.frame(s=info[,'s'],la=info[,'la']),by = join_by(s))
plusG<-alkBasis$a>alkBasis$la
# summary(plusG)
alkBasis$wf<-1
alkBasis[plusG,'wf']<-exp(-(alkBasis[plusG,'a']-alkBasis[plusG,'la'])*0.5) #quick and dirty plusgroup
alkBasis[plusG,'a']<-alkBasis[plusG,'la']
tst1<-filter(alkBasis,s==2&a>6)
#tst1


alkBasis<-alkBasis %>% group_by(area,year,y ,q, sp,s, a ,size,sizeClass) %>% 
  summarize(meanLength=weighted.mean(x=meanLength, w=wf),alk=weighted.mean(x=alk, w=wf),.groups = "drop") %>%
  group_by(area,year,y ,q, sp,s, a ) %>%
  mutate(alk=alk/sum(alk))  %>% ungroup()
#tst1
#filter(alkBasis,s==2&a>6)
  
fl<- alkBasis %>% group_by(area,y,q,s,a) %>% summarize(minSize=min(sizeClass), maxSize=max(sizeClass),.groups = "drop")

#tapply(fl$maxSize,list(fl$s,fl$y),max)
#tapply(fl$minSize,list(fl$s,fl$y),min)
maxSizeCl<-max(fl$maxSize)
alk<- left_join(alkBasis %>% select(area,y,q,s,a,sizeClass,meanLength,alk),fl,by = join_by(area, y, q, s, a)) %>% mutate(s=as.integer(s))
if (FALSE) {
  chk<-filter(alk,y==40,q==3,s==2)
  chk<-xtabs(alk~a+sizeClass,data=chk)
  round(cbind(chk,rowSums(chk)),3)
}
#alk$alkObsNo<-1L:dim(alk)[[1]]
#aBasis$stomObsNo=1L:dim(aBsis)[[1]]

#check that ALK and STOM match
alk <- alk %>% mutate(yqALK=paste(y,q,sep='-'))

alkyq<-alk %>% select(area,y,q,yqALK) %>% unique() %>% arrange(area,y,q) 
aBasis<-left_join(aBasis,alkyq,by = join_by(area, y, q))
#glimpse(aBasis)


a1<-sort(unique(alk$yqALK))
a2<-sort(unique(aBasis$yqALK))
stopifnot(length(setdiff(a2,a1))==0) # missing ALKs
a1a2<-setdiff(a1,a2) #There are ALK for y and q combinations that are not used
#a1a2
if (length(a1a2)>0) {
  alk<-alk %>% filter(yqALK %in% a2)
}


# check there is at least on ALK observation for a predator size class
stomChkPred<-aBasis %>% filter(pred<=nSpecies) %>%select("area","year","y","q","predC","pred","predSizeClass") %>% unique() %>% mutate(fromStom=TRUE)
alkChkPred<-alk %>%  select("area","y","q","s","sizeClass") %>% unique() %>% mutate(fromALK=TRUE)
chkPred<-full_join(x=stomChkPred,y=alkChkPred,by=join_by(area==area,y==y,q==q,pred==s, predSizeClass==sizeClass))
# summary(stomChkPred);summary(alkChkPred); summary(chkPred)
#print(filter(chkPred,is.na(fromStom)),n=20) #  ALK but no STOM, not a problem
problem<-filter(chkPred,is.na(fromALK))
print(problem,n=20) #  STOM but no ALK, that is a problem

aBasis<-anti_join(aBasis,problem%>%mutate(fromStom,fromALK)) 
#stopifnot(dim(problem)[[1]]==0)

# check there is at least on ALK observation for a prey size class
stomChkPrey<-aBasis %>% select("area","year","y","q","preyC","prey","preySizeClass") %>% filter(preyC !=otherFoodName) %>%unique()%>% mutate(fromStom=TRUE)
alkChkPrey<-alk %>%  select("area","y","q","s","sizeClass") %>% unique() %>% mutate(fromALK=TRUE)
chkPrey<-full_join(x=stomChkPrey,y=alkChkPrey,by=join_by(area==area,y==y,q==q,prey==s, preySizeClass==sizeClass))
problem<-filter(chkPrey,is.na(fromALK))
#summary(stomChkPrey);summary(alkChkPrey);summary(chkPrey)
print(problem,n=200) #  STOM but no ALK, that is a problem
#dim(problem)
stopifnot(dim(problem)[[1]]==0)

alkyq<-alk %>% select(area,y,q,yqALK) %>% unique() %>% arrange(area,y,q) %>% mutate(yqALK=1L:dplyr::n())
aBasis$yqALK<-NULL
aBasis<-left_join(aBasis,alkyq,by = join_by(area, y, q))
#summary(aBasis)
filter(aBasis,is.na(yqALK) & pred <=nSpecies ) %>% select(area,y,year,q,pred,predC) %>% unique()

stom<- aBasis %>% select(area, y,q,yqALK,pred, predSizeClass, predSizeW ,noSampl, phi,nPreyGroups,prey,preySizeClass, preySizeW, logPPsize, type, firstL, stomconTot ,stomcon,preyIdx,vulneraIdx,stomObsVarIdx) %>% 
  nest(data=c(pred,predSizeClass,predSizeW,stomObsVarIdx,noSampl,   phi,nPreyGroups,prey,preySizeClass,preySizeW,logPPsize,type,stomcon, stomconTot ,preyIdx,firstL,vulneraIdx)) %>% 
  rowwise() %>% mutate(data=list(data %>% nest(data=c(prey,preySizeClass,preySizeW,logPPsize,type,firstL,preyIdx,stomconTot ,stomcon,vulneraIdx))))

if (FALSE){
  str(stom$data[[1]],2)
  str(stom$data[[1]][['data']],2)
  stom[,'y']
  stom[1,'data'][[1]][[1]]
  stom[1,'data'][[1]][[1]][['data']]
  stom[1,'data'][[1]][[1]][['data']][[1]]
  stom[1,'data'][[1]][[1]][['data']][[2]]
  stom[1,'data'][[1]][[1]]$data[[3]]
  

  stst<-unnest(stom,cols = c(data))
  stst<-unnest(stst,cols = c(data))
  names(stst)
  tst<-select(stst,area,y,q,pred,predSizeClass,prey,preySizeClass)
  dim(stst);dim(tst); dim(unique(tst))
  
}


alk<-alk%>%  nest(data=c(a,s,minSize,maxSize,sizeClass,meanLength,alk)) %>% 
  rowwise() %>% mutate(data=list(data %>% nest(data=c(sizeClass,meanLength,alk))))



#str(alk$data[[1]],2)
#alk[1,'data'][[1]][[1]]
#alk[1,'data'][[1]][[1]][['data']]




#other food
of<-scan(file.path(dir,"other_food.in"),comment.char = "#",quiet=TRUE)
of<-head(of,-1)
names(of)<-c(spNames,othspNames)
of
} #end if multi

#Un<-matrix(0.0, nrow=sum(info[,"last-age"]-sms@first.age+1L),ncol=nYears )
Un<-matrix(0.0, nrow=sum(nlogN),ncol=nYears )

Uf<-matrix(0.0, nrow=length(unlist(sms@keyLogFsta)),ncol=nYears)


#return values
rl<-list(
   data=list(                                          
     options=sms,
     sms.mode=sms@VPA.mode,                            #
     info=info,                                        # Various information for each species with analytically assessment
     nSpecies=nSpecies,                                # Number of species with analytically assessment
     preds=predators,                                  # predator species
     preys=preys,                                      # prey species
     nOthSpecies=nOthSpecies,                          # Number of other predators (used for multispecies mode) 
     otherFoodn=otherFoodn,                            # Species number for other food
     nYears=nYears,                                    # Number of years used in the model
     nAges=nAges,
     nSeasons=nSeasons,                                # number of seasons
     fqa=fqa,                                          # fist quater for each age
     yqIdx=yq_idx,                                     # index for year-season combinations
     minAge=as.integer(sms@first.age),                 # A vector of integers defining the the lowest age class in the assessment for each species.
     #maxAge=as.integer(sms@max.age.all),              # Maximum age for all species. The actual age range by species is given in "info"
     recAge=as.integer(sms@first.age),                 # recruitment age
     years=sms@first.year.model:sms@last.year.model,   # A vector of the years used in the model
     spNames=spNames,                                  # Species names of species with analytically assessment
     othspNames=othspNames,                            # Species names of "other" species (without analytically assessment)
     allSpNames=allSpNames, 
     allSpNamesLong=allSpNamesLong,
     predNames=predNames,
     preyNames=preyNames,
     otherFoodName=otherFoodName,
     fleetNames=fleetNames,                            # names of survey fleets
     spawnSeason=1L,
     recSeason=recSeason,                              # recruitment season
     fbarRange=fbarRange,                              # Minimum and maximum age by species used to calculate Fbar 
     off.age=off.age,                                  # offset between age and the age index used in RSMS (where fist age index  have to be 1 )
     off.year=off.year,                                # offset between year and the year index used in RSMS (where fist year index  have to be 1 )
     #off.season=off.season,
     #off.species=off.species,
     #off.oths=off.oths,                          
     seasonalCatches=seasonalCatches,                  # are seasonal catches available
     doProcessNoise=inclProcess,
     doProcessN_any=inclProcess %in% c(1,2),
     doProcessN_old=inclProcess %in% c(1),
     doProcessN_none=inclProcess==3,
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
     logNfirstYparamfromTo=logNfirstYparamfromTo,
     logNrecruitParamfromTo=logNrecruitParamfromTo,
     
     keyVarLogN=keyVarLogN,
     keyVarObsCatch=keyVarObsCatch,
     minVarObsCatch=sms@min.catch.CV,
     minVarObsSurvey=minCVsurvey,
     keyVarObsSurvey=keyVarObsSurvey,
     keyCatchability=keyCatchability,
     propMat=  propMat,      
     stockMeanWeight=stockMeanWeight,
     catchMeanWeight=catchMeanWeight,
     catchNumber=catchNumber,
     seasFprop=seasFprop,
     logSeasFprop=logSeasFprop,
     ageQCatchMis=ageQCatchMis,
     natMor= natMor,
     propF=zero,
     propM=zero,
     
     
     keyCatch=keyCatch,
     catchNoSeason=catchNoSeason,
     catchNoSeasonUsed=catchNoSeasonUsed,
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
      logNrecruitParam=logNrecruitParam,
      logNfirstYparam=logNfirstYparam,
      rho=rho,
      rec_loga=rec_loga, 
      rec_logb=rec_logb,
       Un=Un,
      Uf=Uf
    ))

if (multi){
  rl[['data']]<-c(rl[['data']],
     list(
      #multi species
      consum=consum,
      meanL=meanL,
      propM2=propM2,          
      natMor1=natMor1,         
      otherN=otherN, 
      stom=stom,
      alk=alk,
      maxSizeCl=maxSizeCl,
      otherFood=of*1000,
     # predPreySizeParm=predPreySizeParm,
      vulneraIdx=vulneraIdx,
      suitIdx=suitIdx,
      partM2Idx=partM2Idx,
      overlap=overlap, 
      overlapIdx=overlapIdx,
      usestomObsVar=usestomObsVar)
  )
      
  rl[['parameters']]<-c(rl[['parameters']],
   list(  
    vulnera=vulnera,
    overlapP=overlapP,
    logStomObsVar=logStomObsVar
   ))
}
return(rl)
}



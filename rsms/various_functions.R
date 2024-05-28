calcCV<-function(sdrep,data) {
  cv_par<-function(par,key,type){
    species<-rownames(key)
    rownames(key)<-paste0('var_',1:dim(key)[[1]])
    s<-data.frame(species,Var1=rownames(key),sp.no=1:length(species))
    a<-array2DF(key)  %>%filter(Value>0) %>% 
      mutate(Var=as.numeric(substr(Var1,5,6)),age=as.numeric(substr(Var2,5,6)),param=par) %>%
      arrange(Var,age)  %>%group_by(Var) %>% mutate(dup=duplicated(Value)) %>%filter(!dup) %>%
      mutate(nexta=lead(age),group=if_else(is.na(nexta),paste0(age,'-'),paste0(age,'-',nexta-1)),type=type) %>% ungroup() 
    a<-   left_join(a,s,by = join_by(Var1)) %>%
      select(sp.no,param,Var,age,group,Value,type)
   if (type=='fleet')   a<-a %>% mutate(sp.no= data$keySurvey.overview[sp.no,'s'])
   a
  }
  
  k<-rbind(
    cv_par(par='logSdLogObsCatch',key=data$keyVarObsCatch,type='species'),
    cv_par(par='logSdLogFsta',key=data$keyLogFstaSd,type='species'),
    cv_par(par="logCatchability",key=data$keyCatchability,type='fleet'),
    cv_par(par="logSdLogObsSurvey",key=data$keyVarObsSurvey,type='fleet'),
    cv_par(par="logSdLogN",key=data$keyVarLogN,type='species'),
    data.frame(param="rho", sp.no=1:data$nSpecies, Var=1:data$nSpecies,age=-9,group='-9',Value=as.integer(data$useRho),type='species') %>% filter(Value>0) %>%
      mutate(Value=cumsum(Value))%>%as_tibble()
  )
  
  pf<-sdrep[['par.fixed']]
  pcv<-data.frame(param=names(pf),parValue=pf,sd=sqrt(diag(sdrep[['cov.fixed']])),gradient=sdrep$gradient.fixed) %>% 
    mutate(cv=round(abs(sd/parValue),3),mark=if_else(cv>0.9,'check','')) %>%
    group_by(param) %>% mutate(n=dplyr::row_number()) %>% ungroup() %>% data.frame()
  pcv
 
  cv<-left_join(pcv,k,by=join_by(param==param,n==Value)) %>% 
    mutate(vari=if_else(type=='species',data$spNames[Var],rownames(data$keySurvey.overview)[Var])) %>%
    select(sp.no,param,parValue,sd, vari,gradient, ages=group) %>%arrange(desc(sd)) %>% as_tibble()
  cv
}


# print(calcCV(sdrep),n=40)



readFleetCatch<-function(dir,of='new_fleet_info.dat', so='fleet_catch.in',condense=TRUE) {
  a<-readNewFleetInfo(dir=dir,of=of,off.age=1,off.year=1)
  
  k<-a[['k']]
  notUse<-k[k[,'useFleet']==0,'f']
  
  nFleets<-nrow(k)
  nSpecies<-length(unique(k[,'s']))
  fleetNames<-rownames(k)
  b<-scan(file=so,quiet=TRUE,comment.char = "#")
  i<-1
  out<-lapply(1:nFleets,function(x) {
    ny<-k[x,"maxy"]-k[x,'miny']+1L
    na<-k[x,"maxage"]-k[x,'minage']+1L
    
    if (k[x,'type'] %in% c(2L,3L)) na<-1L
    if (k[x,'type'] %in% c(4L,5L)) na<-0L
    ii<-ny*(na+1)
    #cat(x,i,ii,ny,na,k[x,'type'],'\n')
    bb<-b[i:(i+ii-1)]
    bb<-matrix(bb,nrow=ny,ncol=na+1,byrow=T)
    
    if (na>0) colnames(bb)<-c('effort',paste('age',k[x,'minage']:k[x,'maxage'],sep='_')) else colnames(bb)<-'effort'
    rownames(bb)<-paste('y',k[x,'minyear']:k[x,'maxyear'],sep='_')
    
    if (condense) {
      if (na>0) bb<- bb[,2:(na+1),drop=FALSE]  /bb[,1]  #cpue
      bb<-array2DF(bb) %>% 
        transmute(year=as.integer(substr(Var1,3,6)),q=k[x,'q'],age=substr(Var2,5,6),s=k[x,'s'],
                  f=k[x,'f'],type=k[x,'type'],obs=Value) %>% as_tibble()
    } 
    i<<-i+ii 
    return(bb)
  })
  if (length(notUse)>0) {
    out<-out[-notUse]
    for (f in (1:length(out))) {
      out[[f]]$f<-f
    }
  }
  if (condense) {
    out<-do.call(rbind,out) 
    out[out$age=='rt','age']<- "-999"
    out$age<-as.integer(out$age)
  }
  out
}

#readFleetCatch(condense=TRUE)
#readFleetCatch(condense=FALSE)



index_old_new<-function(dir,outDir=dir,sms.dat='sms.dat',of="new_fleet_info.dat",verbose=FALSE) {
  if (FALSE) {
    dir=data.path; sms.dat='rsms.dat'; seasFfrom<-'catch'; of="new_fleet_info.dat"
  }
  
  Init.function(dir=dir,sms.dat=sms.dat) # initialize SMS environment
  #cat(SMS.control@species.names,'\n') # just cheking
  
  sms<-SMS.control  # just shorter name
  info<-sms@species.info
  multi.environment<- any(info[,'predator']==2)
  
  nSpecies<-sms@no.species-first.VPA+1L  # species with analytical assessment
  nOthSpecies<-0L
  ages<-sms@first.age:sms@max.age.all
  nAges<-length(ages)
  recSeason<-as.integer(sms@rec.season)
  nSeasons<-as.integer(sms@last.season)
  years<-sms@first.year.model:sms@last.year.model; nYears<-length(years)
  info<-sms@species.info[first.VPA:sms@no.species,c("last-age", "first-age F>0", "last-age-selec", "last-age-likelihood", "+group", "SSB/R","RecAdd2"),drop=FALSE]
  spNames<-dimnames(info)[[1]]
  othspNames<-dimnames(info)[[1]][1:(first.VPA-1)]
  
  
  # read survey data and options into FLR objects
  indices<-SMS2FLIndices(control=sms,path=file.path(root,dir))
  #indices[[1]]@range.SMS
  
  nFleets<-length(indices)
  spNo<-unlist(lapply(indices,function(x) x@range.SMS["species"]))
  PowerAge<-unlist(lapply(indices,function(x) x@range.SMS["power.age"]))
  season<-unlist(lapply(indices,function(x) x@range.SMS["season"]))
  PowerAge[PowerAge<0]<- -9
  PowerAge[PowerAge>=0]<- PowerAge[PowerAge>=0] +1L
  fleetNames<-unlist(lapply(indices,function(x) x@name))
  q.age<-unlist(lapply(indices,function(x) x@range.SMS['q.age']))
  var.age<-lapply(indices,function(x) x@range.SMS['var.age.group'])
  a<-do.call(data.frame,lapply(indices,function(x) x@range))
  a<-rbind(a,s=spNo,q.age=q.age,f=1:nFleets,PowerAge=PowerAge,q=season)
  a["plusgroup",]<-0L
  ra<-rownames(a)
  ra[1]<-'minage'
  ra[2]<-'maxage'
  rownames(a)<-ra
  
  a<- a%>% mutate_if(is.numeric,as.integer)
  a<-t(a)[,c('f','s','minyear',"maxyear","minage","maxage","q.age","plusgroup","PowerAge","q","startf","endf")]
  a<-cbind(a,type=1L,techCreep=0L)
  
  if (is.vector(a)) {lena<-length(a); nama=names(a); a<-matrix(a,nrow=1,ncol=lena); colnames(a)<-nama}
  rownames(a)<-paste(1:nFleets,fleetNames)
  
  keySurvey<-a
  keySurvey.df<-as.data.frame(a) %>% tibble::rownames_to_column("fName") %>% mutate(species=spNames[s])
  
  ## fleets per species
  a<-keySurvey.df %>% group_by(s,species) %>% summarize(n=dplyr::n()) %>% data.frame()
  
  cat("# File for configuration of survey observations\n",file=of,append=F)
  cat(0.2, " # minimum CV of CPUE observations\n",file=of,append=T)
  cat("### number of fleets by species\n",file=of,append=T) 
  for (i in (1:dim(a)[[1]])) cat(a[i,'n']," # ",a[i,'species'],'\n',file=of,append=T)
  
  cat("####### survey options ###### \n",
      "#      Fleet name\n",
      "# 1:   Use the fleet asssessment, 1=yes, 0=no",
      "# 2:   Index type and input:\n",
      "#        1 = Effort and index catch in numbers by age;\n",
      "#        2 = Effort and index for total exploitable biomass, (including all ages with F>0);\n",
      "#        3 = Effort and index for toal spawning stock biomass;\n",
      "#        4 = Commercial effort used with the assumption that Effort is a proxy for F at age;\n",
      "#        5 = Commercial effort used with the assumption that Effort is a proxy for Fbar ;\n",
      "# 3-4:  First and last year\n",
      "# 5:    Season for survey, use 1 for annual data\n",
      "# 6.7:  Alpha and beta - the start and end of the survey period for the fleet given as fractions of the season\n",
      "#                  (or of the year if annual data are used)\n",
      "# 8-9:  First and last age, use -1 if not relevant\n",
      "# 10:    Last age is a plus group? (0=no, 1=yes)\n",
      "# 11:   Last age for stock size dependent catchability (power model), -1 shows no ages uses power model\n",
      "# 12:   Technical creep in effort (only used for type 4 index)\n",
      "# 13:   Number  of cathability groups\n",
      "#         First ages in group of ages  for estimates of catcahability at age\n",
      "#       Number of variance of cathability groups\n",
      "#         First ages in group of ages  for estimates of variance of catcahability at age\n",
      "#       Check value (must be -999)\n",
      "########################",file=of,append=T)
  
  txt<-FALSE 
  ff<-0
  for (sp in (1:dim(a)[[1]])) {
    cat("#########  ",a[sp,'species'],'\n',file=of,append=T)
    fls<-filter(keySurvey.df,sp==s) 
    for (i in (1:dim(fls)[[1]])) {
      ff<-ff+1
      cat("##  ",fls[i,'species'],"# Fleet:",i,fleetNames[fls[i,'f']], "  fleet",ff,'\n',file=of,append=T)
      cat(fls[i,'useFleet'],"\t\t# Use the fleet\n",file=of,append=T)
      cat(fls[i,'type'],"\t\t# Index type\n",file=of,append=T)
      cat(fls[i,'minyear'],fls[i,'maxyear'],"\t# First and last year\n",file=of,append=T)
      cat(fls[i,'q'],"\t\t# Season for survey\n",file=of,append=T)
      cat(fls[i,'startf'],fls[i,'endf'],"\t\t# Alpha and beta\n",file=of,append=T)
      cat(fls[i,'minage'],fls[i,'maxage'],"\t\t# First and last age\n",file=of,append=T)
      cat(0,"\t\t# Is last age a plus group (0=no, 1=yes)\n",file=of,append=T)
      cat(fls[i,'PowerAge'],"\t\t# Last age power model\n",file=of,append=T)
      cat("0 \t\t# Technical creep\n",file=of,append=T)
      nag<-fls[i,'minage']:fls[i,'q.age']
      if (fls[i,'q.age']<fls[i,'maxage']) nag<-c(nag,tail(nag,1)+1)
      cat(length(nag),"\t\t# Number of cathability groups\n",file=of,append=T)
      cat("  ",nag,if_else(length(nag)<3,"\t\t","\t"),"# First ages, catchability groups\n",file=of,append=T)
      va<-var.age[[ff]][[1]]
      cat(length(va),"\t\t# Number of variance of catchability groups\n",file=of,append=T)
      cat("  ",va,if_else(length(va)<3,"\t\t","\t"),"# First ages, variance of catchahability at age\n",file=of,append=T)
      cat("-999 \t\t# check value (must be -999)\n",file=of,append=T)
    }
  }
 cat("file: ",of,"has been written\n")
}

# index_old_new(dir,outDir=dir,sms.dat='sms.dat',of=file.path(root,dir,"new_fleet_info.dat")) 



readNewFleetInfo<-function(dir=dir,of='new_fleet_info.dat',off.age,off.year,verbose=FALSE){
  
  # dir=data.path; of='new_fleet_info.dat'; off.age=1;off.year=-1973; verbose=FALSE
  
  Init.function(dir=dir,sms.dat='sms.dat') # initialize SMS environment
  #cat(SMS.control@species.names,'\n') # just checking
  sms<-SMS.control  # just shorter name
  nSpecies<-sms@no.species-first.VPA+1L  # species with analytically assessment
  
  a<-scan(file=file.path(dir,of),comment.char = "#",quiet=TRUE)
  cv<-a[1]
  a<-a[-1]
  nf<-a[1:nSpecies]
  a<-a[-(1:nSpecies)]
  
  nFleets<-sum(nf)
  # Fleetnames
  s<-readLines(file.path(dir,'fleet_names.in'), n=1000)
  s<-gsub('_',' ',s)
  s<-sub('[[:space:]]+$', '', s)
  s<-s[substr(s,1,1)!='#']
  flNames<-s[1:nFleets]
  
  
  k<-matrix(-999L,nrow=nFleets,ncol=18)
  colnames(k)<-c("f","s","minyear","maxyear","miny","maxy","minage","maxage","mina","maxa","plusgroup","PowerAge", "q","startf","endf","type","techCreep","useFleet")
  rownames(k)<-flNames
  
  k[,'f']<-1:nFleets
  fl<-0
  catchability<-list()
  catchability.var<-list()
  
  for (s in (1:nSpecies)) {
    if (verbose) cat('reading species:',s,'\n')
    for (nfl in(1:nf[s]) ) {
      if (verbose) cat('\treadingfleet:',nfl)
      i<-1
      fl<-fl+1
      k[fl,'s']<-s
      k[fl,'useFleet']<-a[i]; i<-i+1
      k[fl,'type']<-a[i]; i<-i+1
      k[fl,'minyear']<-a[i]; i<-i+1
      k[fl,'maxyear']<-a[i]; i<-i+1
      k[fl,'q']<-a[i]; i<-i+1
      k[fl,'startf']<-a[i]; i<-i+1
      k[fl,'endf']<-a[i]; i<-i+1
      k[fl,'minage']<-a[i]; i<-i+1
      k[fl,'maxage']<-a[i]; i<-i+1
      k[fl,'plusgroup']<-a[i]; i<-i+1
      k[fl,'PowerAge']<-a[i]; i<-i+1
      k[fl,'techCreep']<-a[i]; i<-i+1
      n.catcha<-a[i];i<-i+1
      catchability<-c(catchability,list(a[i:(i+n.catcha-1)])); i<-i+n.catcha
      n.catcha<-a[i];i<-i+1
      catchability.var<-c(catchability.var,list(a[i:(i+n.catcha-1)])); i<-i+n.catcha
      chk<-a[i]; 
      if (verbose) cat(' done, check=:',chk,'\n')
      stopifnot(chk== -999)
      a<-a[-(1:i)]
      
    }
  } 
  
  k[,'mina']<-k[,'minage']+off.age
  k[,'maxa']<-k[,'maxage']+off.age
  k[,'miny']<-k[,'minyear']+off.year
  k[,'maxy']<-k[,'maxyear']+off.year
  #mode(k)<-'integer'
  list(k=k,cv=cv,catchability=catchability,catchability.var=catchability.var)
}

# j<-readNewFleetInfo(dir=data.path,of='new_fleet_info.dat',off.age=data$off.age,off.year=data$off.year,verbose=TRUE); str(j,1)


writeNewFleetInfo<-function(dir=dir,of='new_fleet_info.dat',key,speciesNames){
  
  a<-key[['k']]
  catcb<-key[['catchability']]
  catcbVar<-key[['catchability.var']]
  key.df<-as.data.frame(a) %>% tibble::rownames_to_column("fName") %>% mutate(species=speciesNames[s])
  
  ## fleets per species
  b<-key.df %>% group_by(s,species) %>% summarize(n=dplyr::n()) %>% ungroup() %>% data.frame()
  
  cat("# File for configuration of survey observations\n",file=of,append=F)
  cat(key$cv, " # minimum CV of CPUE observations\n",file=of,append=T)
  cat("### number of fleets by species\n",file=of,append=T) 
  for (i in (1:dim(b)[[1]])) cat(b[i,'n']," # ",b[i,'species'],'\n',file=of,append=T)
  
  cat("####### survey options ###### \n",
      "#      Fleet name\n",
      "# 1:   Use the fleet asssessment, 1=yes, 0=no\n",
      "# 2:   Index type and input:\n",
      "#        1 = Effort and index catch in numbers by age;\n",
      "#        2 = Effort and index for total exploitable biomass, (including all ages with F>0);\n",
      "#        3 = Effort and index for toal spawning stock biomass;\n",
      "#        4 = Commercial effort used with the assumption that Effort is a proxy for F at age;\n",
      "#        5 = Commercial effort used with the assumption that Effort is a proxy for Fbar ;\n",
      "# 3-4:  First and last year\n",
      "# 5:    Season for survey, use 1 for annual data\n",
      "# 6.7:  Alpha and beta - the start and end of the survey period for the fleet given as fractions of the season\n",
      "#                  (or of the year if annual data are used)\n",
      "# 8-9:  First and last age, use -1 if not relevant\n",
      "# 10:    Last age is a plus group? (0=no, 1=yes)\n",
      "# 11:   Last age for stock size dependent catchability (power model), -1 shows no ages uses power model\n",
      "# 12:   Technical creep in effort (only used for type 4 index)\n",
      "# 13:   Number  of cathability groups\n",
      "#         First ages in group of ages  for estimates of catcahability at age\n",
      "#       Number of variance of cathability groups\n",
      "#         First ages in group of ages  for estimates of variance of catcahability at age\n",
      "#       Check value (must be -999)\n",
      "########################",file=of,append=T)
  
  txt<-FALSE 
  ff<-0
  for (sp in (1:dim(b)[[1]])) {
    cat("#########  ",b[sp,'species'],'\n',file=of,append=T)
    fls<-filter(key.df,s==sp) 
    for (i in (1:dim(fls)[[1]])) {
      ff<-ff+1
      cat("##  ",fls[i,'species'],"# Fleet:",i,fls[i,'fName'], "  fleet",ff,'\n',file=of,append=T)
      cat(fls[i,'useFleet'],"\t\t# Use the fleet\n",file=of,append=T)
      cat(fls[i,'type'],"\t\t# Index type\n",file=of,append=T)
      cat(fls[i,'minyear'],fls[i,'maxyear'],"\t# First and last year\n",file=of,append=T)
      cat(fls[i,'q'],"\t\t# Season for survey\n",file=of,append=T)
      cat(fls[i,'startf'],fls[i,'endf'],"\t\t# Alpha and beta\n",file=of,append=T)
      cat(fls[i,'minage'],fls[i,'maxage'],"\t\t# First and last age\n",file=of,append=T)
      cat(0,"\t\t# Is last age a plus group (0=no, 1=yes)\n",file=of,append=T)
      cat(fls[i,'PowerAge'],"\t\t# Last age power model\n",file=of,append=T)
      cat("0 \t\t# Technical creep\n",file=of,append=T)
      nag<-length(catcb[[ff]])
      cat(nag,"\t\t# Number of cathability groups\n",file=of,append=T)
      cat("  ",catcb[[ff]],if_else(length(nag)<3,"\t\t","\t"),"# First ages, catchability groups\n",file=of,append=T)
      va<-catcbVar[[ff]]
      cat(length(va),"\t\t# Number of variance of catchability groups\n",file=of,append=T)
      cat("  ",va,if_else(length(va)<3,"\t\t","\t"),"# First ages, variance of catchahability at age\n",file=of,append=T)
      cat("-999 \t\t# check value (must be -999)\n",file=of,append=T)
    }
  }
  cat("file: ",of,"has been written\n")
}

# a<-readNewFleetInfo(dir=data.path,of='new_fleet_info.dat',off.age=data$off.age,off.year=data$off.year)
#a<-readNewFleetInfo(dir=data.path,of='new_fleet_info.dat',off.age=1L,off.year=-1973L)

#writeNewFleetInfo(dir=data.path,of='new_fleet_info3.dat',key=a,speciesNames=data$spNames)
#writeNewFleetInfo(dir=data.path,of='new_fleet_info3.dat',key=a,speciesNames=spNames)



announce<-function(x) {
  cat("\nobjective:",x$objective,"  convergence:",x$convergence, "  ", x$message, "  iterations:",x$iterations, "  evaluations:",x$evaluations)
}

inputToDF<-function(data) {
  a<-rbind(
    list_rbind(lapply(data$catchMeanWeight,array2DF),names_to='s') %>% mutate(var="weca"),
    list_rbind(lapply(data$catchNumber,array2DF),names_to='s') %>% mutate(var="canum"),
    list_rbind(lapply(data$stockMeanWeight,array2DF),names_to='s') %>% mutate(var="west"),
    list_rbind(lapply(data$propMat,array2DF),names_to='s') %>% mutate(var="propMat"),
    list_rbind(lapply(data$propM,array2DF),names_to='s') %>% mutate(var="propM"),
    list_rbind(lapply(data$propF,array2DF),names_to='s') %>% mutate(var="propF"),
    list_rbind(lapply(data$natMor,array2DF),names_to='s') %>% mutate(var="M"),
    list_rbind(lapply(data$seasFprop,array2DF),names_to='s') %>% mutate(var="seasFprop")

  )
  pivot_wider(a,names_from=var,values_from=Value) %>% mutate_if(is.character,as.integer) %>%
    mutate(species=data$spNames[s], year=Var1-data$off.year, q=Var2,age=Var3-data$off.age) %>% 
    mutate(Var1=NULL,Var2=NULL,Var3=NULL)
}


#a<-inputToDF(data)


outputAnnuToDF<-function(obj) {
  rep<-obj$report()
  a<-rbind(
    list_rbind(lapply(rep$Chat,array2DF),names_to='s') %>% mutate(var="Chat")
  )  
  pivot_wider(a,names_from=var,values_from=Value)  %>%
    mutate(species=s, year=as.integer(y),age=parse_number(a)-data$off.age) %>% 
    mutate(y=NULL,a=NULL,Var3=NULL,s=NULL,Chat=exp(Chat))
} 
#  b<-outputAnnuToDF(obj)

outputToDF<-function(obj) {
  rep<-obj$report()
  a<-rbind(
    list_rbind(lapply(rep$logNq,array2DF),names_to='s') %>% mutate(var="N",Value=exp(Value)),
    list_rbind(lapply(rep$Zq,array2DF),names_to='s') %>% mutate(var="Z")
   )  
     pivot_wider(a,names_from=var,values_from=Value)  %>%
    mutate(species=s, year=as.integer(y),q=as.integer(q),age=parse_number(a)-data$off.age) %>% 
    mutate(y=NULL,a=NULL,Var3=NULL,s=NULL,deadZ=N*(1-exp(-Z)))
} 
#  b<-outputToDF(obj)


FToDF<-function(obj,sdrep,data) {
  if (missing(sdrep)) sdrep <- sdreport(obj)
  x<-as.list(sdrep, what="Est")
  ff<-exp(x$Uf)
  colnames(ff)<-data$years
  FF<-NULL
  for (s in (1:data$nSpecies)) {
    i<-data$nlogFfromTo[s,]
    key<-data$keyLogFsta[s,][data$keyLogFsta[s,]>0]
    faf<-data$info[s,'faf']; la<-data$info[s,'la']
    fff<-ff[i[1]:i[2],,drop=FALSE][key,]
    if (data$info[s,'fModel']==2) {
      fSepar<-data$info[s,'fSepar']; la<-data$info[s,'la']
      key2<-t(data$keyLogSeparF[[s]][,fSepar:la] )
      dims<-dim(key2)
      seaF<-exp(x$logSeparF[key2])
      seaF<-matrix(seaF,ncol=dims[[2]])
      fff[fSepar:la,]<-fff[fSepar:la,]*seaF
    }
    rownames(fff)<-(faf:la)-data$off.age
    fff<-array2DF(fff) %>% mutate(s=s,year=as.integer(Var2),age=as.integer(Var1),Var1=NULL,Var2=NULL) %>%rename(FF=Value)
    
    FF<-rbind(FF,fff)
  }
  FF
}

seasonalF<-function(obj,sdrep,data) {
  FF<-FToDF(obj,sdrep,data) 
  a<-inputToDF(data) %>% select(s,seasFprop, species, year, q, age)
  left_join(FF,a,by = join_by(s, year, age)) %>% mutate(SF=FF*seasFprop)
}
#  FF<-FToDF(obj,sdrep=sdrep)

#b<-left_join(inputToDF(data),outputToDF(obj),by = join_by(species, year, q, age))

plotF<-function(obj,sdrep,data,combineAges=FALSE) {
  ff<- FToDF(obj,sdrep,data) %>% mutate(species=data$spNames[s])
  
  a<-t(data$keyLogFsta)
  a<-cbind(a,Age=data$minAge-1 +(1:data$nAges))
  a<-data.frame(a) %>% pivot_longer(cols=1:data$nSpecies) %>%rename(species=name,aGroup=value,age=Age)
  
  ff<-left_join(ff,a) %>% as_tibble()
  
  ag<-a%>%  group_by(species,aGroup) %>% summarise(mina=min(age),maxa=max(age)) %>% ungroup() %>%
    mutate(ages=paste(mina,maxa,sep='-'),mina=NULL,maxa=NULL)
  ag
  
  ff<-left_join(ff,ag,by = join_by(species, aGroup))%>% mutate(ages=factor(ages),age=factor(age))
  
  fplt<-by(ff,ff$s,function(x) ggplot(x,aes(x=year,y=FF,shape=age,col=age))+
             geom_point(size=2)  + geom_line()+labs(title=x[1,'species'],ylab='F')+
             facet_wrap(vars(data$spNames[s]),ncol=2,scales="free_y")
  )
  return(fplt)
  if (combineAges) {
    fff<-ff %>%mutate(age=NULL) %>% unique()
    
    by(fff,fff$s,function(x) ggplot(x,aes(x=year,y=FF,shape=ages,col=ages))+
         geom_point(size=2)  + geom_line()+labs(title=x[1,'species'],ylab='F')+
         facet_wrap(vars(data$spNames[s]),ncol=2,scales="free_y")
    )
  }
}


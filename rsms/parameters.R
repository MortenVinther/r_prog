extractParameters<-function(sdrep,myMap,data) {
  a<-data.frame(name=attr(sdrep$par.fixed,'names'),estimate=sdrep$par.fixed,estimate.sd=sqrt(diag(sdrep$cov.fixed)))
  if (length(dim(sdrep$gradient.fixed))==0)   a$gradient<-sdrep$gradient.fixed else  a$gradient<-sdrep$gradient.fixed[1,] 
  a<-a %>%mutate(order=1:dim(a)[[1]])
  a$key<-ave(seq_len(nrow(a)), a$name, FUN = seq_along)
  
  transKey<-function(key,pp,type) {
    array2DF(key) %>% filter(Value>0) %>%mutate(param=pp,type=type)
  }
  transKeyNoNA<-function(key,pp,type) {
    array2DF(key)%>%filter(!is.na(Value))%>% mutate(param=pp,type=type)
  }
  
  parName<-unique(a$name)
  keys<-NULL
  
  if ('logSdLogObsCatch' %in% parName) keys<- rbind(keys,
                                                    left_join(data$keyVarObsCatch,data.frame(s=1:data$nSpecies,seasonal= data$seasonalCatches),by = join_by(s))  %>% 
                     transmute(Var1=if_else(seasonal==1,paste0(data$spNames[s],' Q:',q),data$spNames[s]) ,Var2=paste('age',a-data$off.age),Value=varGroup,param='logSdLogObsCatch',type='species'))
  
  if ('logFSeasonal' %in% parName) {
    k<-do.call(rbind,lapply(1:data$nSpecies,function(s) if (data$FfromSeparableModel[s]==1) array2DF(data$keylogSeasonF[[s]]) %>% filter(Value>0) %>%
             transmute(s,y=as.integer(Var1),q=as.integer(Var2),age=as.integer(Var3)-data$off.age,grp=as.integer(Value)) %>% 
             group_by(s,q,grp) %>%summarize(fa=min(age),la=max(age),fy=min(y)-data$off.year,ly=max(y)-data$off.year,.groups='drop') %>%
             transmute(Var1=paste0(data$spNames[s],' Q:',q,' ',fy,'-',ly), Var2=paste(fa,la,sep='-'), Value=grp,param="logFSeasonal",type='qData') else NULL))

    k<-k %>% arrange(Value) %>%mutate(used= !is.na(myMap$logFSeasonal)) %>% filter(used) %>% mutate(used=NULL,Value=1:dplyr::n())
    
   keys<-rbind(keys,k)
  }                                                  
  
  if ('logCatchability' %in% parName) keys<-rbind(keys,transKey(key=data$keyCatchability,pp='logCatchability', type='fleet'))
  if ('logSdLogObsSurvey' %in% parName) keys<-rbind(keys,transKey(key=data$keyVarObsSurvey,pp='logSdLogObsSurvey', type='fleet'))
 
  if ('logSdLogFsta' %in% parName) {
    k<-data$keyLogFsta;k
    kk<-cumsum(apply(k,1,max))
    if (dim(k)[[1]]>=2) for (i in (2:dim(k)[[1]])) k[i,k[i,]>0]<-k[i,k[i,]>0]+kk[i-1]
    keys<-rbind(keys,transKey(key=k,pp='logSdLogFsta', type='species'))
  }

  if ('logYearEffectF' %in% parName) {
    x<-parameters$logYearEffectF
    x[!is.na( myMap$logYearEffectF)]<-1:sum(!is.na(myMap$logYearEffectF))   
    x[is.na( myMap$logYearEffectF)]<-NA
    transKeyNoNA(key=x,pp="logYearEffectF",type='species')
    keys<-rbind(keys,transKeyNoNA(key=x,pp="logYearEffectF",type='species'))
  }
  
  
 
  if ('logNrecruitParam' %in% parName) keys<-rbind(keys,do.call(rbind,lapply(data$spNames,function(x) 
    data.frame(Var1=x,Var2=seq(1,data$logNrecruitParamfromTo[x,2]-data$logNrecruitParamfromTo[x,1]+1)-data$off.year,
               Value=seq(data$logNrecruitParamfromTo[x,1],data$logNrecruitParamfromTo[x,2]),param='logNrecruitParam',type='species'))) %>% 
      filter(Value >0))
  
 
  if ('logNfirstYparam' %in% parName) keys<-rbind(keys,do.call(rbind,lapply(data$spNames,function(x) {
    if (data$logNfirstYparamfromTo[x,1] >0 ) {
      ages<-data$logNfirstYparamfromTo[x,2]-data$logNfirstYparamfromTo[x,1]+1
      k<-data.frame(Var1=x,Var2=paste('age',1:ages),Value=seq(data$logNfirstYparamfromTo[x,1],data$logNfirstYparamfromTo[x,2])) %>%
      mutate(param='logNfirstYparam',type='species')
    } else return(NULL)
  }))) 
    
  if ('logSeparAgeF' %in% parName) {
    k<-data$keylogSeparAgeF
    x<-lapply(1:data$nSpecies,function(x) {
     kk<-unique(k[[x]])
     if (any(kk!= -1)) {
       array2DF(kk) %>% filter(Value>0) %>%
         transmute(Var1=paste(data$spNames[x],Var1),Var2,Value,param='logSeparAgeF',type='fleet')
     } 
    })
    keys<-rbind(keys,  do.call(rbind,x))
  }
  
  if ('logSsbRsd' %in% parName) {
    k<-data$inclSsbR;
    keys<-rbind(keys,data.frame(Var1=names(k),Var2=-9,Value=1:length(k),param='logSsbRsd',type='species'))
  }
  
  
 
  if ('logSdLogN' %in% parName) {
    keys<-rbind(keys,transKey(key=data$keyVarLogN,pp='logSdLogN', type='species'))
  }
  
  
  if ('vulnera' %in% parName) keys<-rbind(keys,transKey(key=data$vulneraIdx,pp='vulnera', type='prey-pred'))

  if ('overlapP' %in% parName) {
    k<- data$overlapIdx %>%  arrange(overlap,predNo,preyNo) %>% 
      transmute(Var1=data$allSpNames[predNo],Var2=c(data$spNames,data$otherFoodName)[preyNo],Value=overlap,param='overlapP',type='prey-pred')
    keys<-rbind(keys,k)
  }
 
  parMap<-function(vari,type='species') {
    y<-filter(a,name==vari)$estimate; l=length(y)
    if (vari %in% names(myMap)) y<-  y[myMap[[vari]]]
    data.frame(Var1=data$spNames[!is.na(y)],Var2=-9,Value=1:l,param=vari,type=type)
  }
  
  parMapPred<-function(vari,type='predator') {
    y<-filter(a,name==vari)$estimate; l=length(y)
    if (vari %in% names(myMap)) y<-  y[myMap[[vari]]]
    data.frame(Var1=data$predNames[!is.na(y)],Var2=-9,Value=1:l,param=vari,type=type)
  }
  
  
  m1<-NULL
  if ('rec_loga' %in% parName) m1<-parMap(vari='rec_loga')
  if ('rec_logb' %in% parName) m1<-rbind(m1,parMap(vari='rec_logb'))
  if ('rho' %in% parName) m1<-rbind(m1,parMap(vari='rho'))
  
  if ('logStomObsVar' %in% parName) m1<-rbind(m1,parMapPred(vari='logStomObsVar',type='predator'))
 
  keys<-ungroup(keys) %>% dplyr::select(Var1,Var2,Value,param,type)
  
  mk<-rbind(m1,keys) %>% rename(key=Value,name=param)
  mk<-mk %>% rowwise() %>% mutate(s=if_else(type %in% c('species','predator','qData'),match(unlist(strsplit(Var1,' '))[1],data$allSpNames),-9))
  tmp<-filter(mk,type=='fleet') %>% transmute(Var1,Var2,name,s2= match(unlist(lapply(strsplit(Var1,' '),function(x) x[1])),data$allSpNames))
  mk<-left_join(mk,tmp,by = join_by(Var1, Var2, name)) %>% mutate(s=if_else(s== -9 & !is.na(s2),s2,s), s2=NULL) 
     
  all<-left_join(a,mk,by = join_by(name, key)) %>%as_tibble()

  tmp<-!(all$type %in% c('prey-pred','predator'))
  all$age<- -9
  all[tmp,'age']<-parse_number(unlist(all[tmp,'Var2']))
  tmp<-!(all$type %in% c('prey-pred'))
  all<-all %>% mutate(Var1=if_else(tmp,Var1,paste(Var1,Var2,sep='-'))) %>% mutate(expEst=if_else(estimate<5,round(exp(estimate),2),-9))
  all2<-all %>%group_by(name,estimate,estimate.sd,gradient,expEst,order,key,Var1,s) %>% 
    summarize(min_age=min(age),max_age=max(age),.groups = "drop")%>%arrange(order) 

    list(parameters=all,parmSummary=all2)
} 


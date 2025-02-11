
simCatchN<-function(inp=inp_all, maxVar=1.0 ^2) {
  # Load the MASS package
  #library(MASS)
  # Set the seed for reproducibility
  set.seed(123)
  
  # test inp=inp_all; maxVar=1.0 ^2
    
  spPlus<-inp[['data']]$nOthSpecies
  a<-Read.summary.data() %>% transmute(year=Year,Quarter,s=Species.n-spPlus,species=Species,age=Age,chat=C.hat) %>%
    filter(chat>=0)
  seasonalC<-inp$data$seasonalCatches
  names(seasonalC)<-inp[['data']]$spNames
  b<-data.frame(s=1:length(seasonalC),seasonal=seasonalC)
  
  a<-left_join(a,b,by = join_by(s)) %>% mutate(quarter=if_else(seasonal==1,Quarter,1)) %>%
    group_by(s,species,year,quarter,seasonal,age) %>% summarize(chat=sum(chat),.groups='drop') 
  
  #parameters from old key run           
  p<-Read.SMS.std(dir=data.path,stdFile="run_ms3.std",gfile='run_ms3.grd',excludeNotUsed=TRUE,rec_scale=TRUE,remove=NULL) %>% as_tibble()
  #sort(unique(p$name))
  p<-p[grep("catch_s2_ini",p$name),]%>%transmute(s=species-spPlus,quarter,age,variance=value) %>% arrange(s,quarter,age)
  aa<-a %>% transmute(s,species,quarter,age,seasonal) %>% unique()
  
  
  a3<-left_join(aa,p,by = join_by(s, quarter, age)) %>% arrange(s,quarter,age)
  
  for (i in 1:10) a3<<-a3%>%   mutate(variance=if_else(is.na(variance),lag(variance,1),variance))
  a4<-a3 %>% filter((seasonal==0 & quarter==1) | (seasonal==1))  %>% arrange(s,quarter,age)
  
  # round(xtabs(sqrt(variance)~paste(s,species)+age,data=filter(a4,quarter==1)),2) #sd
  
  a5<-left_join(a,a4,by = join_by(s, species, quarter, seasonal, age)) %>% mutate(variance=if_else(is.na(variance),2,variance))
  a5<-data.frame(a5)
  
  b<-by(a5,list(a$s,a$species,a$year,a$quarter),function (x) {
    # x<-filter(a,year==1974 & s==1)
    nage<-dim(x)[[1]]
    chat<-log(x[,'chat'])
    cov<-matrix(0,ncol=nage,nrow=nage)
    diag(cov)<-min(x[,'variance'],maxVar)
    x$Csim<-exp(MASS::mvrnorm(n=1, mu = chat, Sigma = cov))
    x
  })
  
  b1<-do.call(rbind,b)
  
  newC<-inp[['data']]$catchNumber
  out<-lapply(names(newC),function(x){
    #cat(x,'\n')
    a<-filter(b1,species==x)
    b<-tapply(a$Csim,list(a$year,a$quarter,a$age),sum)
    b[is.na(b)]<-0
    if (seasonalC[x])  newC[[x]][,,]<<-b else {
      for (q in as.numeric(dimnames(b)[[2]]) ) newC[[x]][,q,]<<-b
    }
    return(x)
  })
  inp[['data']]$catchNumber<-newC
  
  template<-data.frame(s=data$keyCatch[,'s'],year=data$keyCatch[,'year'],age=data$keyCatch[,'age'],obs.no=data$keyCatch[,'obs.no'])
  b2<-b1 %>% group_by(s,year,age) %>% summarize(Csim=sum(Csim),.groups='drop')
  b2<-left_join(template,b2,by = join_by(s, year, age)) %>% arrange(obs.no)
  inp[['data']]$logCatchObs<-log(b2$Csim)
  
  return(inp)
}

#newC<-simCatchN(inp=inp_all, maxVar=1.0 ^2)


sdrep <- sdreport(obj); 
cat('Hesssian:',sdrep$pdHess,'\n')

sdrep

a<-data.frame(name=attr(sdrep$par.fixed,'names'),value=sdrep$par.fixed,gradient=sdrep$gradient.fixed) 
a<-a %>%mutate(order=1:dim(a)[[1]])
a$key<-ave(seq_len(nrow(a)), a$name, FUN = seq_along)
head(a)



transKey<-function(key,pp,type) {
  array2DF(key) %>% filter(Value>0) %>%mutate(param=pp,type=type)
}

keys<-           transKey(key=data$keyVarObsCatch, pp='logSdLogObsCatch',type='species')
keys<-rbind(keys,transKey(key=data$keyCatchability,pp='logCatchability', type='fleet'))
keys<-rbind(keys,transKey(key=data$keyVarObsSurvey,pp='logSdLogObsSurvey', type='fleet'))
keys<-rbind(keys,transKey(key=data$keyLogFsta,pp='logSdLogFsta', type='species'))
keys<-rbind(keys,transKey(key=data$keyVarLogN,pp='logSdLogN', type='species'))


sdrep$value
opt$par

parMap<-function(vari,type='species') {
  y<-filter(a,name==vari)$value; l=length(y)
  if (vari %in% names(my.map)) y<-  y[my.map[[vari]]]
  data.frame(Var1=data$spNames[!is.na(y)],Var2=-9,Value=1:l,param=vari,type=type)
}

m1<-parMap(vari='rec_loga')
m1<-rbind(m1,parMap(vari='rec_logb'))
m1<-rbind(m1,parMap(vari='rho'))
m1
head(keys)
mk<-rbind(m1,keys) %>% rename(key=Value,name=param)
head(mk)
head(a)
all<-left_join(a,mk)

all%>% mutate(numVar=parse_number(Var2))  %>% arrange(name,key,numVar)

group_by(name,key) %>% 
  summarize(min_a=min(numVar),max_a=max(numVar)) %>%ungroup() %>% mutate(name=pp,type=type)
summary(all)

transKeyXXX<-function(key,pp,vars,type) {
  y<-array2DF(key); 
  y<-filter(y,Value>0) %>%rename(key=Value) %>% mutate(age=parse_number(Var2)) %>% group_by(Var1,key) %>% 
    summarize(min_a=min(age),max_a=max(age)) %>%ungroup() %>% arrange(key) %>% mutate(name=pp,type=type)
  left_join(vars,y,by = join_by(name, key))
}


keys

data$keyLogFsta

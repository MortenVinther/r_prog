# Root Mean square Error, RMSE
library(plotrix)

fleet.names<-Read.fleet.names()

file<-file.path(data.path,'fleet_info.dat') 
finfo<-scan(file,comment.char = "#") 

file<-file.path(data.path,'catch_survey_residuals.out')
res<-read.table(file,comment.char = "#",header=T)

res<-tibble(subset(res,res$residual!=-99.9))

par_exp1<-read.table(file.path(data.path,"par_exp.out"),header=TRUE) %>% filter(par=="qq_s2_ini") %>%
  rename(Species.n=species,groupAge=age) %>% mutate(Age=groupAge) %>% 
  dplyr::select(parNo,Species.n,fleet,Age,groupAge)

pars<-scan(file.path(data.path,"sms.par"),comment.char = "#")
par_exp1$value<-sqrt(pars[par_exp1$parNo])
ress<-filter(res,data=='survey')

surv<-left_join(ress,par_exp1)  %>% 
      fill(parNo) %>% fill(groupAge) %>% fill(value) %>%
      group_by(Species.n, fleet, Quarter,groupAge) %>% mutate(mina=min(Age),maxa=max(Age)) %>% ungroup() %>%
      mutate(groupAge=paste(mina,maxa,sep=':'),mina=NULL,maxa=NULL)
surv

# just to check that the original parameter values are the same as RMSE 
check<-surv %>% group_by(Species.n,fleet,Quarter,groupAge,value) %>% summarise(n=dplyr::n(),RMSE=sqrt(sum(residual^2)/n)) %>%
  ungroup()%>% mutate(difval=value-RMSE)
check

filter(check,difval>0.000001) # should be empty

tst<-subset(surv,Species.n==2 & groupAge=="1:3" &  fleet==1)

taylor.diagram(ref=tst$observed,model=tst$model,add=FALSE)
taylor.diagram(ref=log(tst$observed),model=log(tst$model),add=FALSE,col='red')

tst<-subset(surv,Species.n==2 & Age==1 & fleet==1)
taylor.diagram(ref=log(tst$observed),model=log(tst$model),add=TRUE,col='blue')
tst<-subset(surv,Species.n==2 & Age==2 & fleet==1)
taylor.diagram(ref=log(tst$observed),model=log(tst$model),add=TRUE,col='orange')
tst<-subset(surv,Species.n==2 & Age==3 & fleet==1)
taylor.diagram(ref=log(tst$observed),model=log(tst$model),add=TRUE,col='green')



#catch

par_exp1<-
  read.table(file.path(data.path,"par_exp.out"),header=TRUE) %>% filter(par=="catch_s2_ini") %>%
  rename(Species.n=species,groupAge=age,Quarter=quarter) %>% mutate(Age=groupAge) %>% 
  dplyr::select(parNo,Species.n,,Quarter,Age,groupAge)
par_exp1


pars<-scan(file.path(data.path,"sms.par"),comment.char = "#")
par_exp1$value<-sqrt(pars[par_exp1$parNo])

ress<-filter(res,data=='catch')
ress
tail(ress)
info<-data.frame(Species.n=first.VPA:nsp,combine_catches=SMS.control@combined.catches,combine_s2=SMS.control@seasonal.catch.s2)
info
ress<-left_join(ress,info)
ress

catch1<-left_join(filter(ress,combine_s2==1),par_exp1) %>% dplyr::arrange(Species.n,Age,Quarter,Year) %>%   mutate(data=NULL,fleet=NULL) %>%
  fill(parNo) %>% fill(groupAge) %>% fill(value) %>%
  group_by(Species.n,  Quarter,groupAge) %>%
  mutate(groupAge=paste(min(Age),max(Age),sep=':')) %>%
  group_by(Species.n) %>% mutate(minq=min(Quarter),maxq=max(Quarter)) %>% ungroup() %>%
  mutate(groupQuarter=if_else(combine_s2==1,paste(minq,maxq,sep=':'),as.character(Quarter))) %>%
           mutate(minq=NULL,maxq=NULL)

catch0<-left_join(filter(ress,combine_s2==0),par_exp1) %>% dplyr::arrange(Species.n,Quarter,Age,Year) %>%   mutate(data=NULL,fleet=NULL) %>%
  fill(parNo) %>% fill(groupAge) %>% fill(value) %>%
  group_by(Species.n,  Quarter,groupAge) %>%
  mutate(groupAge=paste(min(Age),max(Age),sep=':')) %>%
  group_by(Species.n) %>% mutate(minq=min(Quarter),maxq=max(Quarter)) %>% ungroup() %>%
  mutate(groupQuarter=if_else(combine_s2==1,paste(minq,maxq,sep=':'),as.character(Quarter))) %>%
  mutate(minq=NULL,maxq=NULL)

catch<-rbind(catch0,catch1)

tail(catch)
filter(catch,Species.n==2& groupAge=="1:8")
filter(catch,Species.n==2& groupAge=="2:8")

# just to check that the original parameter values are the same as RMSE 
check<-catch %>% group_by(Species.n,,groupQuarter,groupAge,value) %>% summarise(n=dplyr::n(),RMSE=sqrt(sum(residual^2)/n)) %>%
  ungroup()%>% mutate(difval=(value-RMSE)^2)
check

filter(check,difval>1E-8) # should be empty

tst<-subset(catch,Species.n==2 & groupAge=="2:8" & groupQuarter=="1:4")

taylor.diagram(ref=log(tst$observed),model=log(tst$model),add=FALSE,col='red')

tst<-subset(surv,Species.n==2 & Age==2 )
taylor.diagram(ref=log(tst$observed),model=log(tst$model),add=TRUE,col='blue')
tst<-subset(surv,Species.n==2 & Age==3 )
taylor.diagram(ref=log(tst$observed),model=log(tst$model),add=TRUE,col='orange')
tst<-subset(surv,Species.n==2 & Age==4 )
taylor.diagram(ref=log(tst$observed),model=log(tst$model),add=TRUE,col='green')



 
Init.function() # get SMS.contol object  including sp.names

nsp<-nsp-first.VPA+1

nox.noy<-nox*noy
plot.no<-0

years<-rep(0,2)
ages<-rep(0,2)
if (is.na(my.species)) my.species<-1:nsp

for (sp in 1:nsp) if (sp %in% my.species) {
  aa<-subset(res,Species.n==first.VPA+sp-1)
  nf<-finfo[sp+1] #no. of fleets
  sp.name<-sp.names[sp+first.VPA-1]
  
  print(paste(sp.name,sp,"Number of fleets:",nf))
  
  for (f in 1:nf) {
    
    print(paste("sp:",sp,"  fleet:",f))
    
    nyr<-years[2]-years[1]+1
    nag<-ages[2]-ages[1]+1
    plot.no<-plot.no+1
    
    if (plot.no%%nox.noy==0 || f==1){
      newplot(dev,nox,noy,Portrait=Portrait,filename=paste('Survey',sp.name,f,plot.no),dir=oDir,...)
      par(mar=c(3,4,3,2))
      if (dev=="wmf") par(mar=c(2,4,2,2))
      plot.no<-0
    }
    
    bb<-subset(aa,fleet==f)
    tmp<-tapply(bb$residual,list(age=bb$Age,year=bb$Year),sum,na.rm=T)
    tmp[tmp==-99.99]<-0
    
    
    xpos <- as.numeric(dimnames(tmp)[[2]]) # years
    ypos <- as.numeric(dimnames(tmp)[[1]]) #ages
    title<- paste(sp.name," fleet:",f,sep="")
    title<- fleet.names[sp,f]
    #title<- paste(sp.name,": ",fleet.names[f],sep="")
    
    if (length(ypos)==1) {
      
      tmp2<-tmp
      tmp2[]<-0
      if (ypos[1]>=1) {
        tmp<-rbind(tmp2,tmp,tmp2)
        ypos <- (ypos-1):(ypos+1)
      }
      else if (ypos[1]==0) {
        tmp<-rbind(tmp,tmp2)
        ypos <- ypos:(ypos+1)
      }      
    }
    #print(tmp)
    if (over.all.max>0) residplot(tmp,xpos,ypos,main=title,refdot=use.ref.dot,start.year=start.year,end.year=end.year,maxn=over.all.max)
    else residplot(tmp,xpos,ypos,main=title,refdot=use.ref.dot,start.year=start.year,end.year=end.year)
    
  }
}

cat("Max buble size=",max.buble,'\n')
if (reverse.colors) cat("log(Survey observed CPUE) - log(expected CPUE). 'Red' negative, 'White' positive\n")
if (! reverse.colors) cat("log(Survey observed CPUE) - log(expected CPUE). 'Red' positive, 'White' negative\n")



s<-read_delim(file.path(stom.input,"stomcon_list_Boots_0500_haul_as_observed.dat"))
s[s$pred=='W_H','pred']<-'WHM'
s[s$pred=='N_H','pred']<-'NHM'

s<-s%>%mutate(tpred=paste(formatC(pred.no,flag='0',w=2),pred),tprey=paste(formatC(prey.no,flag='0',w=2),prey),prey.mean.length.ALK=NULL,pred.no=as.integer(pred.no),prey.no=as.integer(prey.no))
#xtabs(~tpred+tprey,data=s)  #OK but old species order



#s<-s %>% mutate(pred.no=oldNewSp_n[pred.no])
s$pred.no<-oldNewSp_n[s$pred.no]
x2<-filter(s,prey=='OTH')
s<-filter(s,prey!='OTH') %>% mutate(prey.no=oldNewSp_n[prey.no])
s<-rbind(s,x2)
s<-s%>%mutate(tpred=paste(formatC(pred.no,flag='0',w=2),pred),tprey=paste(formatC(prey.no,flag='0',w=2),prey))
xtabs(~tpred+tprey,data=s) #OK



pp<-matrix(0L,nrow=nSpecies+1,ncol=nSpecies+nOthSpecies,dimnames=list(c(spNames,'OTH'),c(spNames,othspNames)))
pp2<- s%>% select(pred,pred.no,prey,prey.no) %>% mutate(prey.no=as.integer(ifelse(prey.no==0,nSpecies+1,prey.no))) %>%unique() 
xtabs(~prey.no+pred.no,data=pp2)
for (i in (1:dim(pp2)[[1]]))  pp[unlist(pp2[i,'prey.no']),unlist(pp2[i,'pred.no'])]<-1L
predPrey<-pp
predPrey

vulneraIdx<-predPrey
vulneraIdx[vulneraIdx>0]<-1:sum(vulneraIdx>0)
vulneraIdx
vulnera<-rep(1,max(vulneraIdx)) #vulnerability coefficient, parameter 


lw<-scan(file.path(dir,"length_weight_relations.in") ,comment.char = "#")
lw<-head(lw,-1)
lw<-matrix(lw,ncol=2,nrow=nSpecies+nOthSpecies,dimnames=list(c(spNames,othspNames),c('a','b')),byrow=TRUE )
lwdf<-data.frame(lw) %>% mutate(sp=1:dim(lw)[[1]])
#lwdf

ss<-s%>% left_join(x=s,y=lwdf,by = join_by(pred.no == sp)) %>% 
  mutate(pred.size.w=a*pred.mean.length**b,a=NULL,b=NULL,prey.size.w=mean.weight,pred.prey.size.w=pred.size.w/prey.size.w)  
 

#xtabs(~tpred+tprey,data=ss)  #OK

library(quantreg)
q1<-0.025 #lower quantile
q2<-0.975 #higher quantile


xtabs(~tpred+tprey,data=ss)  #OK

a<-filter(ss,type=='obs' & prey !='OTH' ) %>% 
  mutate(log.pred.prey.size.w=log(pred.size.w/prey.size.w),log.pred.size.w=log(pred.size.w)) %>%
  select(SMS_area,pred.no, pred,tpred,prey.no,prey,tprey,log.pred.size.w,log.pred.prey.size.w ) 
a
xtabs(~tpred+tprey,data=a)  #OK

# predator size independent prey range
aa<-a %>% group_by(pred.no,pred,tpred,prey.no,prey,tprey) %>% summarize(min_s=min(log.pred.size.w),max_s=max(log.pred.size.w))
round(xtabs(min_s~tpred+tprey,data=aa),3)  
round(xtabs(max_s~tpred+tprey,data=aa),3)  
aa

#
pp<-array(0,dim=c(4,nSpecies,nSpecies+nOthSpecies),dimnames=list(c('lower_intc','lower_slope','upper_intc','upper_slope'),spNames,c(spNames,othspNames)) )                                                         
for (i in (1:dim(aa)[[1]])) {
  pred<-unlist(aa[i,'pred'])
  prey<-unlist(aa[i,'prey'])
  pp[1,prey,pred]<-unlist(aa[i,'min_s'])
  pp[3,prey,pred]<-unlist(aa[i,'max_s'])
}
round(ftable(pp),2)


#quantile regressions

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

a1<-filter(a,pred.no %in% c(1:9) )
n<-a1 %>% group_by(pred.no,prey.no) %>% summarize(n=dplyr::n())
an<-left_join(a1,n,by = join_by(pred.no, prey.no))

aa<-filter(an,n>5) %>% 
  nest_by(pred.no,pred,prey.no,prey) %>%  
  mutate(ru=list(rq(log.pred.prey.size.w~log.pred.size.w, tau=q1,data=data)),
         rl=list(rq(log.pred.prey.size.w~log.pred.size.w, tau=q2,data=data)))


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
round(ftable(pp),2)  
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
  ftable(pp)
}

## stomach contents

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
                   noStom=as.integer(stom.no),
                   noHaul=as.integer(haul.no),
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

  
a<- aBasis %>% mutate(predSize=NULL,  predMeanLength=NULL, predSizeW=NULL,,preySize=NULL, preySizeClass=NULL, preyMeanLength=NULL ) %>% 
  group_by( area,year,y,q,predC,pred, predSizeClass,noHaul,noStom,phi) %>% 
  nest()

if (FALSE) {
  head(a)
  a[1,'data'][[1]]
  a[6,'data'][[1]]
  a[6,'data'][[1]][[1]]
  
  print(a[6,'data'][[1]][[1]],n=50)
}

stom<-a

#other food
of<-scan(file.path(dir,"other_food.in"),comment.char = "#")
of<-head(of,-1)
names(of)<-c(spNames,othspNames)
of

stomD<-list(stom=stom,of=of,predPredSize=pp,vulneraIdx=vulneraIdx,vulnera=vulnera)



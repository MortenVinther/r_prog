selPred=c(1L,2L)

suitIdx<-unnest(data$suitIdx,cols = c(data))
suitIdx
suitIdx$data[[1]]
s<-filter(suitIdx,predNo %in% selPred) %>% unnest(cols = c(data)) %>% filter(preyNo !=data$otherFoodn)  %>% 
  mutate(prey=factor(data$spNames[preyNo]))
s 
summary(s)
xtabs(~prey+preyAge+q,data=s)

stom<-unnest(data$stom,cols = c(data))
stom
stom<-filter(stom,pred %in% selPred) %>% unnest(cols = c(data)) %>% filter(prey !=data$otherFoodn)  %>% 
  transmute(area,y,q,preyNo=prey,prey=data$spNames[preyNo],predNo=pred,pred=data$predNames[predNo],predW=predSizeW,logRatio=logPPsize)
stom 
summary(stom)
sort(names(stom));sort(names(s)); 

ss<-rbind(stom %>% select(area,y,q,preyNo,predNo,predW,logRatio) %>% mutate(type='stom'),
            s %>% select(area,y,q,preyNo,predNo,predW,logRatio) %>% mutate(type='M2'))  %>%
     mutate(prey=data$spNames[preyNo],pred=data$predNames[predNo])
ss

#library(plyr)
sss<-filter(ss,predNo==1)
by(sss,sss$prey,function(x) by(x,x$type,function(xx) range(xx$logRatio)))
by(sss,sss$prey,function(x) by(x,x$type,function(xx) exp(range(xx$logRatio))))

mu <- ddply(sss, "type", summarise, grp.mean=mean(logRatio))
mu
ggplot(sss, aes(x=logRatio, fill=type, color=type)) +
  geom_histogram(position="identity", alpha=0.5)+
  facet_wrap(~prey,  ncol=2, strip.position = 'right')
ggplot(sss, aes(x=logRatio, fill=type, color=type)) +
  geom_histogram(position="identity", alpha=0.5)+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=type),
             linetype="dashed")+
  facet_wrap(~paste(prey,type),  ncol=2, strip.position = 'right',scales='free_y')



out<-by(s,list(s$predNo),function(x) {
  p=unlist(x[1,'predNo'])
  st<-filter(stom,predNo==p)
  tit<- data$predNames[(x$predNo)[1]] 
  a<-ggplot(x,aes(x=predW, y = logRatio, col=prey,shape=prey)) +
    facet_wrap(~paste0('Q:',q),  ncol=2, strip.position = 'right')+
    geom_point( )+
    geom_point(data=st,aes(x=log(predW), y = logRatio), col='black') +
    labs(x='log (Predator size)',y='log(predator/prey size ratio)',title=tit)
    #theme_minimal() +
    #theme( panel.grid.major = element_line(linetype = "blank"),
    #       panel.grid.minor = element_line(linetype = "blank")
    #)
    print(a)
})

cleanup()
out<-by(s,list(s$predNo,s$preyNo),function(x) {
  st<-filter(stom,predNo==unlist(x[1,'predNo']) & preyNo==unlist(x[1,'preyNo']))
  tit<- paste(data$predNames[(x$predNo)[1]], 'eating',data$spNames[(x$preyNo)[1]]) 
  a<-ggplot(x,aes(x=predW, y = logRatio, col=prey,shape=prey)) +
    facet_wrap(~paste0('Q:',q),  ncol=2, strip.position = 'right')+
    geom_point( )+
    geom_point(data=st,aes(x=log(predW), y = logRatio), col='black') +
    labs(x='log (Predator size)',y='log(predator/prey size ratio)',title=tit)
  #theme_minimal() +
  #theme( panel.grid.major = element_line(linetype = "blank"),
  #       panel.grid.minor = element_line(linetype = "blank")
  #)
  print(a)
})

cleanup()

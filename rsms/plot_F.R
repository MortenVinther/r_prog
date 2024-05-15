ff<- FToDF(obj,sdrep) %>% mutate(species=data$spNames[s])

a<-t(data$keyLogFsta)
a<-cbind(a,Age=data$minAge-1 +(1:data$nAges))
a<-data.frame(a) %>% pivot_longer(cols=1:data$nSpecies) %>%rename(species=name,aGroup=value,age=Age)

ff<-left_join(ff,a) %>% as_tibble()

ag<-a%>%  group_by(species,aGroup) %>% summarise(mina=min(age),maxa=max(age)) %>% ungroup() %>%
  mutate(ages=paste(mina,maxa,sep='-'),mina=NULL,maxa=NULL)
ag

ff<-left_join(ff,ag,by = join_by(species, aGroup))%>% mutate(ages=factor(ages),age=factor(age))

ff


by(ff,ff$s,function(x) ggplot(x,aes(x=year,y=FF,shape=age,col=age))+
  geom_point(size=2)  + geom_line()+labs(title=x[1,'species'],ylab='F')+
  facet_wrap(vars(data$spNames[s]),ncol=2,scales="free_y")
)

fff<-ff %>%mutate(age=NULL) %>% unique()

by(fff,fff$s,function(x) ggplot(x,aes(x=year,y=FF,shape=ages,col=ages))+
     geom_point(size=2)  + geom_line()+labs(title=x[1,'species'],ylab='F')+
     facet_wrap(vars(data$spNames[s]),ncol=2,scales="free_y")
)


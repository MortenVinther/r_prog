b<-left_join(inputToDF(data),outputToDF(obj),by = join_by(species, year, q, age))
b
summary(b)

ff<-FToDF(obj,sdrep=sdrep)

ab<-left_join(b,ff,by = join_by(s, year, age)) %>%
  mutate(PopeF=canum/(N*exp(-M/2)),FF=seasFprop*FF,annualC=data$combinedCatches[s]==1)

data$combinedCatches

san<-filter(ab,!annualC & q %in% c(2,3) & age>0) %>% group_by(s,species,year,age) %>%
  mutate(fProp=PopeF/sum(PopeF),fPropCanum=canum/sum(canum),q=factor(q))
san

x<-filter(san,s==2) %>% select(-weca,-west,-propMat,-propM,-propF,-Z,-M    )
x

ggplot(data=x, aes(x=fPropCanum, y=fProp, group=q)) +
  geom_point(aes(shape=,col=q))+
  facet_grid(vars(age), scales="free_y")+
  ggtitle(unlist(x[1,'species']))

ggplot(data=x, aes(x=seasFprop, y=fProp, group=q)) +
  geom_point(aes(shape=,col=q))+
  facet_grid(vars(age), scales="free_y")+
  ggtitle(unlist(x[1,'species']))

ggplot(data=x, aes(x=year, y=seasFprop, group=q)) +
  geom_point(aes(shape=,col=q))+
  facet_grid(vars(age), scales="free_y")+
  ggtitle(unlist(x[1,'species']))

ggplot(data=x, aes(x=year, y=fProp, group=q)) +
  geom_point(aes(shape=,col=q))+
  geom_line(aes(shape=,col=q))+
  facet_grid(vars(age), scales="free_y")+
  ggtitle(unlist(x[1,'species']))



#fleet with commercial effort (assumed to cover the full international fishery)
efl<-data$keySurvey.overview[data$keySurvey.overview[,'type']==4,]

eff<-lapply(1:dim(efl)[[1]],function(i) {
  x<-filter(data.frame(data$keySurvey),f==efl[i,'f'] & s==efl[i,'s']) %>% select(obs.no,f,s,y,q) %>% 
    mutate(effort=exp(data$logSurveyObs[obs.no])) %>% as_tibble()
})

eff<-do.call(rbind,eff) %>% mutate(obs.no=NULL,year=y-data$off.year)

yield<-inputToDF(data) %>%group_by(s,species,year,q) %>% summarize(yield=sum(canum*west))

eff
yield
cpue<-left_join(eff,yield,by = join_by(s, q, year)) %>% mutate(cpue=yield/effort/1000,fleet=factor(f))

ggplot(cpue,aes(x=year,y=effort,fill=fleet,shape=fleet,col=fleet))+
  geom_point(size=2)+
  geom_smooth(method=lm)+  facet_wrap(vars(data$spNames[s]),ncol=2,scales="free_y")

ggplot(cpue,aes(x=year,y=yield,fill=fleet,shape=fleet,col=fleet))+
  geom_point(size=2)+
  geom_smooth(method=lm)+  facet_wrap(vars(data$spNames[s]),ncol=2,scales="free_y")

ggplot(cpue,aes(x=year,y=cpue,fill=fleet,shape=fleet,col=fleet))+
  geom_point(size=2)+
  geom_smooth(method=lm)+  facet_wrap(vars(data$spNames[s]),ncol=2,scales="free_y")

  

catch_compare_hat<-function(
    indir=data.path,
    paper=TRUE,
    compare.dir=data.path,
    makeAllGraphs=FALSE # make plots for HTML output
    )
{
  if (makeAllGraphs) dev<-'png' else dev='screen'
  filename='Catch_compare_hat.png'
  if (makeAllGraphs) filename=file.path(compare.dir,filename) 
  png(file=filename,width = 800, height = 800)
  
  faF<-data.frame(Species.n=1:nsp,faF=SMS.control@species.info[,"first-age F>0"])
                  
  dat<-Read.summary.data(dir=indir,extend=FALSE,read.init.function=F) %>% filter(weca>=0) 
  dat<-left_join(dat,faF) %>% filter(Age>=faF) %>%
     mutate(observed=Yield,predicted=C.hat*weca) %>% group_by(Species.n,Species,Year) %>% 
    summarize(observed=sum(observed,na.rm=TRUE),predicted=sum(predicted,na.rm=T)) %>% ungroup() %>%
    pivot_longer( cols = c("observed", "predicted"),names_to="Yield")
  
  x<-ggplot(dat, aes(x=Year, y=value/1000,shape=Yield,color=Yield)) + 
     geom_point()+ geom_line()+facet_wrap(vars(Species),scales="free_y")+
     ylab("Yield (1000 tonnes)")+theme_minimal(base_size = 18)+  theme(legend.position="bottom")+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5 ))
  print(x)
 if (paper) cleanup() 
}

#catch_compare_hat(indir=data.path,makeAllGraphs=FALSE)



catch_compare_age_hat<-function(
    indir=data.path,
    paper=TRUE,
    compare.dir=data.path,
    makeAllGraphs=FALSE # make plots for HTML output
){
  if (makeAllGraphs) dev<-'png' else dev='screen'
  filename=file.path('Catch_compare_hat_age.png')
  if (makeAllGraphs) filename=file.path(compare.dir,filename)  
  ageFormat<-data.frame(a=(0:15),lab=paste('Age',(0:15))) 
 
  faF<-data.frame(Species.n=1:nsp,faF=SMS.control@species.info[,"first-age F>0"])
  
dat<-Read.summary.data(dir=indir,extend=FALSE,read.init.function=F) %>% filter(weca>=0) 
dat<-left_join(dat,faF) %>% filter(Age>=faF) %>%
  mutate(observed=Yield,predicted=C.hat*weca, Age=factor(Age,levels=ageFormat$a,labels=ageFormat$lab)) %>% 
  group_by(Species.n,Species,Year,Age) %>% 
  summarize(observed=sum(observed,na.rm=TRUE),predicted=sum(predicted,na.rm=T)) %>% ungroup() %>%
  pivot_longer( cols = c("observed", "predicted"),names_to="Yield")

pp<-unique(dat$Species.n)
for (p in pp) {
  dat3<-  dat %>% filter(Species.n==p)
  if (makeAllGraphs) dev<-'png' else dev<-'screen'
  
  filename=paste0('Catch_compare_hat_age_',dat3[1,'Species'],'.png')
  if (makeAllGraphs) filename=file.path(compare.dir,filename)  
  png(file=filename,width = 800, height = 800)
  
  x<-ggplot(dat3, aes(x=Year, y=value/1000,shape=Yield,color=Yield)) + 
    geom_point()+ geom_line()+facet_wrap(vars(Age),scales="free_y")+
    labs(title=dat3[1,'Species'])+
    ylab("Yield (1000 tonnes)")+theme_minimal(base_size = 18)+   theme(legend.position="bottom")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5 ))
  print(x)
  if (paper) cleanup() 
}
}


#catch_compare_age_hat(indir=data.path,makeAllGraphs=FALSE)



overw<-as.data.frame(data$keySurvey.overview) %>%  rownames_to_column(var = "fleet") %>% select(f,fleet,s,type)

survResid<-as.data.frame(data$keySurvey) %>% as_tibble()  %>% 
  mutate(predict=rep$predSurveyObs,obs=data$logSurveyObs,Residual=predict-obs,year=y-data$off.year,Age=factor(a-data$off.age),species=data$spNames[s]) %>% 
  left_join(.,overw,by = join_by(f, s)) %>% select(s,species,f,fleet,type,year,q,Age,s,predict,obs,Residual) %>%
  mutate(cols=if_else(Residual<0,'negativ','positiv')) 

plt<-by(survResid,survResid$s,function(x){
  plt<-x %>% ggplot(aes(year, Age,color= cols,fill=cols,size = sqrt(abs(Residual)))) +
    geom_point(shape = 21,alpha=0.75) +
    scale_color_manual(values = c("blue", "red")) +
    labs(x = "", y = "Age",title=paste('Survey residuals, ',x[1,'species'])) +
    # facet_grid(rows =vars(fleet), scales="free_y")
    facet_wrap(vars(fleet),ncol=1)+
    theme_bw() + 
    theme(strip.text = element_text(size = 10,margin = margin(0.01,0,0.01,0, "cm")))
  print(plt)
})
print(plt)

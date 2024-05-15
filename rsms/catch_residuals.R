b1<-outputAnnuToDF(obj)
b2<-inputToDF(data) %>% group_by(s,species,year,age) %>% summarize(canum=sum(canum)) %>% ungroup()
faf<-data.frame(species=data$spNames,faf=  data$info[,'first-age F>0'])
Resid<-left_join(b1,b2, by = join_by(species, year, age)) %>% left_join(.,faf,by = join_by(species)) %>% filter(age>=faf)  %>%
  mutate(Residual=log(Chat)-log(canum), cols=if_else(Residual<0,'negativ','positiv'),Age=factor(age))

plt<-by(Resid,Resid$s,function(x){
  plt<-x %>% ggplot(aes(year, Age,color= cols,fill=cols,size = sqrt(abs(Residual)))) +
    geom_point(shape = 21,alpha=0.75) +
    scale_color_manual(values = c("blue", "red")) +
    labs(x = "", y = "Age",title=paste('Catch residuals, ',x[1,'species'])) +
    theme_bw() + 
    theme(strip.text = element_text(size = 10,margin = margin(0.01,0,0.01,0, "cm")))
  print(plt)
})
print(plt)

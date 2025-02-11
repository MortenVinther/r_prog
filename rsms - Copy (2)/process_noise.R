
N<-lapply(rep$logNq,function(x) (cbind(a1=x[,data$recSeason,1],x[,1,-1]))) # N in recruiting season for the first age and age 1 jan for the rest
# Process N
resid<-lapply(data$spNames,function(x) {
  y<-t(N[[x]])-rep$predN[[x]]; 
  rownames(y)<-(1:dim(y)[[1]])-data$off.age
  y<-y[,-1]
  attr(y,'species')<-x; 
  y
})

plt<-lapply(resid,function(x) {
  as.data.frame(x) %>%
    rownames_to_column(var = "Age") %>%
    gather(key, Residual, -Age) %>% as_tibble() %>%
    mutate(Age=factor(as.integer(Age)),year=as.integer(key),cols=if_else(Residual<0,'negativ','positiv')) %>% #mutate(cols=factor(cols)) %>%
    ggplot(aes(year, Age,color= cols,fill=cols,size = sqrt(abs(Residual)))) +
    geom_point(shape = 21,alpha=0.75) +
    scale_color_manual(values = c("blue", "red")) +
    labs(x = "", y = "Age",title=paste('Process noise,  ',attr(x,'species'))) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5))
})

print(plt)


Standardized.residuals<-function(dir=data.path) {
 
  file<-file.path(dir,'catch_survey_residuals.out')
  res<-read.table(file,comment.char = "#",header=T) %>% filter(data=='catch' & residual != -99.9 ) %>% mutate(Species=paste(Species.n,sp.names[Species.n]))
 
  catch<-ggplot(data = res) + 
    geom_histogram(mapping = aes(x = stand.residual, y=..density..),color="darkblue", fill="lightblue") + 
    facet_wrap(~Species,nrow=4)+
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1))
    #stat_function(fun = dnorm, args = list(mean = mean(res$stand.residual), sd = sd(res$stand.residual)))
 
  # survey
  
  file<-file.path(dir,'catch_survey_residuals.out')
  res<-read.table(file,comment.char = "#",header=T) %>% filter(data=='survey' & residual != -99.9 ) %>% mutate(Species=sp.names[Species.n])
 
  survey<-ggplot(data = res) + 
    geom_histogram(mapping = aes(x = stand.residual, y=..density..),color="darkblue", fill="lightblue") + 
    facet_wrap(~Species,nrow=4)+
    stat_function(fun = dnorm, args = list(mean = 0, sd = 1))
 
  return(list(catchPlot=catch,surveyPlot=survey))
}

# test a<-Standardized.residuals();lapply(a,print)

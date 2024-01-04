source(file.path(rsms.root.prog,"make_rsms_data_function.R"))

dat<-make_rsms_data(dir="S23",annual=F,outDir=rsms.root)

dat$data$spNames
if (FALSE) {
  str(dat)
  dat$data$keyLogFsta
  dat$data$keyVarObsCatch
  
  head(dat$data$keyCatch )
  
  dat$data$nlogFfromTo  
  dat$data$keyVarLogN 
  
  round(ftable(dat$data$seasFprop[[1]][45:47,,]),3)
  dat$data$recAge


  logcatch<-dat$data$logCatch
  catch<-as.data.frame(dat$data$keyCatch)  %>% mutate(catch=exp(logcatch)) %>% as_tibble()
  round(xtabs(catch~ age+year+s,data=catch))
}

#dat$data$seasFprop[[1]][10:20,3,]
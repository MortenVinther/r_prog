source(file.path(rsms.root.prog,"make_rsms_data_function.R"))

dat<-make_rsms_data(dir="S27",annual=T)
dat$data$spNames

str(dat)
dat$data$keyLogFsta
dat$data$keyVarObsCatch

head(dat$data$keyCatch )

dat$data$nlogFfromTo  
dat$data$keyVarLogN 

if (FALSE) {
  logcatch<-dat$data$logCatch
  catch<-as.data.frame(dat$data$keyCatch)  %>% mutate(catch=exp(logcatch)) %>% as_tibble()
  round(xtabs(catch~ age+year+s,data=catch))
}

#dat$data$seasFprop[[1]][10:20,3,]

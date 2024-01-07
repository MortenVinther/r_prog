source(file.path(rsms.root.prog,"make_rsms_data_function.R"))

dat<-make_rsms_data(dir="S21",annual=F,outDir=rsms.root)

#dat<-make_rsms_data(dir="ns_2023_ss_input",annual=F,outDir=rsms.root)
#dat<-make_rsms_data(dir="S16_S17_S21",annual=F,outDir=rsms.root)
dat$data$spNames

data<-dat[['data']]
parameters<-dat['parameters']
save(file=file.path(rsms.root,"rsms_input.Rdata"))
if (FALSE) {
  str(dat)
  dat$data$keyLogFsta
  dat$data$keyVarObsCatch
  dat$data$zeroCatchYear
  head(dat$data$keyCatch )
  
  dat$data$nlogFfromTo  
  dat$data$keyVarLogN 
  
  round(ftable(dat$data$seasFprop[[1]][45:47,,]),3)
  dat$data$recAge


  logcatch<-dat$data$logCatch
  catch<-as.data.frame(dat$data$keyCatch)  %>% mutate(catch=exp(logcatch)) %>% as_tibble()
  round(xtabs(catch~ age+year+s,data=catch))
}


source(file.path(root.prog,"r_prog","rsms","make_rsms_data.R"))

dat<-make_rsms_data(dir="S17",annual=TRUE)
dat$data$spNames

str(dat)

logcatch<-dat$data$logCatch
catch<-as.data.frame(dat$data$keyCatch)  %>% mutate(catch=exp(logcatch)) %>% as_tibble()
catch
summary(catch)

round(xtabs(catch~ age+year+s,data=catch))

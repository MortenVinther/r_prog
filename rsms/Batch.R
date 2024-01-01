# run ini.r in SMS dir as the first

combSp<-c("S16","S17","S18","S19","S20","S21","s22","S23")



#### then
my.comb<-"S16"
library(RTMB)
library(tidyverse)
sam.root<-file.path("~","cod");
rsms.root<-file.path("~","cod","RSMS");

for (my.comb in combSp) {
  source(file.path(root.prog,"r_prog","rsms","make_rsms_data.R"))
  
  make_rsms_data(dir=my.comb,annual=TRUE)
  # makes  save(data,parameters,file=file.path(rsms.root,"rsms_input.Rdata"))
  load(file.path(rsms.root,"rsms_input.Rdata"),verbos=T)
  cat(data$spNames,'\n')
  source(file.path(root.prog,"r_prog","rsms","rsms.R")) 
  save(obj,opt,file=file.path(sam.root,paste0(my.comb,'.Rdata')))
}

for (my.comb in combSp) {
  load(file.path(sam.root,paste0(my.comb,'.Rdata')))
  cat("Comb :",my.comb,'\n')
  cat("objective:",opt$objective,"  convergence:",opt$convergence, "  # 0 indicates successful convergence\n")
  
  rep<-obj$report()
  print(rep$nlls)
  
}

my.comb

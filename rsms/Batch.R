# run ini.r in SMS dir as the first

#### then
my.comb<-"S16"
library(RTMB)
library(tidyverse)
sam.root<-file.path("~","cod");
rsms.root<-file.path("~","cod","RSMS");

source(file.path(root.prog,"r_prog","rsms","make_rsms_data.R"))

make_rsms_data(dir=my.comb,annual=TRUE)
# makes  save(data,parameters,file=file.path(rsms.root,"rsms_input.Rdata"))
load(file.path(rsms.root,"rsms_input.Rdata"),verbos=T)
cat(data$spNames,'\n')


source(file.path(root.prog,"r_prog","rsms","rsms.R"))       

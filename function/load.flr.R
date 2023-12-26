
library(FLCore)
#library(FLAssess)
 #  library(FLEDA)
if (allLibraries) { 
  #ibrary(FLXSA)
  #library(FLBRP)
  #library(FLHCR)
  #library(FLEDA)

 }

do.test<-FALSE


#library(FLSMS)

#overwrite FLSMS library
source(file.path(my.FLR.path,"flsms.control.R"))
source(file.path(my.FLR.path,"flsms.predict.control.R"))

source(file.path(my.FLR.path,"flstockmulti.R"))
source(file.path(my.FLR.path,"flindex.sms.R"))
source(file.path(my.FLR.path,"xsa2sms.R"))
source(file.path(my.FLR.path,"sms2flsmss.R"))

source(file.path(my.FLR.path,"flsms.R"))
source(file.path(my.FLR.path,"flsmss.R"))

source(file.path(my.FLR.path,"flop.control.R"))
source(file.path(my.FLR.path,"flop_msfd.control.R"))




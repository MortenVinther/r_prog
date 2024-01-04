#Remove all objects
#rm(list = ls())

library(RTMB)
library(tidyverse)

# Operating System
OS<- .Platform$OS.type
#  directory for SMS, runs
if (OS=="unix")   home<-"~"  else home<-file.path("C:","_C_drev")

sam.root<-file.path(home,"RSMS");
rsms.root<-file.path(home,"RSMS")
#rsms.root.prog<-file.path(home,'SMS','r_prog','rsms')
rsms.root.prog<-file.path(root.prog,'r_prog','rsms')


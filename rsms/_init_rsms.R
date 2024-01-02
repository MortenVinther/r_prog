#Remove all objects
#rm(list = ls())

library(RTMB)
library(tidyverse)

# Operating System
OS<- .Platform$OS.type
# Harddisk drive for SMS, runs
if (OS=="unix")   home<-"~"  else home<-"C:"

sam.root<-file.path(home,"cod");
rsms.root<-file.path(home,"cod","RSMS")
rsms.root.prog<-file.path(home,'SMS','r_prog','rsms')


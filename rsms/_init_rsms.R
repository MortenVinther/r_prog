#Remove all objects
#rm(list = ls())

library(RTMB)
library(tidyverse)

rsms.root.prog<-file.path(root.prog,'r_prog','rsms')

stom.input<-file.path(root,"rsms_input","stom_input")

source(file.path(rsms.root.prog,"various_functions.R"))
source(file.path(rsms.root.prog,"rsms.control.R"))
source(file.path(rsms.root.prog,"summary_plot.R"))


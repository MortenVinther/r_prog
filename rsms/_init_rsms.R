#Remove all objects
#rm(list = ls())

library(RTMB)
library(tidyverse)
library(ggpubr)
#library(AICcmodavg)

rsms.root.prog<-file.path(root.prog,'r_prog','rsms')

stom.input<-file.path(root,"rsms_input","stom_input")

source(file.path(rsms.root.prog,"various_functions.R"))
source(file.path(rsms.root.prog,"rsms.control.R"))
source(file.path(rsms.root.prog,"parameters.R"))
source(file.path(rsms.root.prog,"plotCompareRuns.R"))
source(file.path(rsms.root.prog,"make_rsms_data_function.R"))
source(file.path(rsms.root.prog,"pick_species.R"))
source(file.path(rsms.root.prog,"quarter2annual.R"))
source(file.path(rsms.root.prog,"rsms_function.R"))
#source(file.path(rsms.root.prog,"summary_plot.R"))
source(file.path(rsms.root.prog,"map_param.R"))
source(file.path(rsms.root.prog,"lowerUpper.R"))
source(file.path(rsms.root.prog,"batch_control.R"))


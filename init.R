#Remove all objects
rm(list = ls())

# Operating System
OS<- .Platform$OS.type
# Harddisk drive for SMS, runs
if (OS=="unix") {
  dosDrive<-"~"
  root<-file.path(dosDrive,"SMS")   # root directory for the SMS package, runs
  root.prog<-root                   # root directory for the SMS package, R-programs (r-prog directory
  sms.exe<-"./sms"
} else {
  dosDrive<-"C:"
  root<-file.path(dosDrive,"_C_drev","SMS-git")
  root.prog<-file.path(dosDrive,"_C_drev","SMS-git")
  sms.exe<-'sms.exe'
}

root.copy<-root
#my.stock.dir<-"NorthSeaKeyRun_2020"

#my.stock.dir<-"NS_2020_APP"
#my.stock.dir<-"NS_2023_01_2020_key_run"


#my.stock.dir<-"NorthSeaKeyRun_2023"
#my.stock.dir<-"NS_2023_key_run_directors_cut"

my.stock.dir<-"rsms_input"

#my.stock.dir<-"rsms_SAN-area-1r"
#my.stock.dir<-"rsms_SAN-area-3r"


my.stock.dir<-"NS_2022_seawise_ver_2024"

# make a backup of the SMS source code
# file.copy(file.path(root.prog,"program","sms.tpl"),file.path(root.prog,"program",paste("sms_",format(Sys.time(), "%Y_%m_%d-%H-%M"),'.tpl',sep='')),overwrite =TRUE)
# file.copy(file.path(root.prog,"program","op.tpl"),file.path(root.prog,"program",paste("op_",format(Sys.time(), "%Y_%m_%d-%H-%M"),'.tpl',sep='')),overwrite =TRUE)

#Installation of FLR, if needed
# install.packages("FLCore", repos="http://R-Forge.R-project.org")

# Sys.unsetenv("GITHUB_PAT")
# remotes::install_github("flr/FLCore")
###################### do not change the code below this line ###############################################################################

makeAllGraphs<-F      # batch job to make "all" graphs and tables after a key-run



# Use all libraries or just a simple configuration (few libraries and limited access to FLR) or full access to
allLibraries<-FALSE

# Path to data directory
data.path<-file.path(root,my.stock.dir)

# Path to R-programme directory, do not change
prog.path<-file.path(root.prog,"r_prog")
my.FLR.path<-file.path(root.prog,"r_prog","flsms")
rsms<-file.path(prog.path,"rsms")

# path to the sms.exe file, used for retrospective runs
#sms.command<-file.path("..","program","sms")

# libraries
library(lattice)
#library(MASS)
#library(gtools)  # for def of macro

if (allLibraries) library(quantreg)
library(tidyverse)
setwd(data.path)

# Path to R-programme directory
prog.path.func<-file.path(root.prog,"r_prog","function")

# for SMS (not RSMS) runs
source(file.path(prog.path.func,"init_r_functions.R"))

cat("active stock directory:",getwd(),"\n");


save(list = ls(all.names = TRUE), file = file.path(data.path,"SMS.Rdata"), envir = .GlobalEnv)



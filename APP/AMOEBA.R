# Script for creating North Sea AMOEBAs
# Modified from a script Written by Jeremy Collie on 19 August 2003


library(tidyverse)
amoeba_dir<-file.path("H:","Amoeba")
data_dir<-file.path(amoeba_dir,'data')


#read as much as possible from Jeremy's S-plus code and data
source(file.path(amoeba_dir,'read_Splus.R')) #returns list, "jeremy" with S-plus code and data 


#Extract (mannually)  the required functions from Jeremy's old S-plus program and put the adapted code in R_functions.R
source(file.path(amoeba_dir,'R_functions.R'))
##########

n.species<-10
n.fleet<-8
eff.change<-1.1 # effort change used to "fit" the model

# fleet names
all_fleet<-read_csv(file=file.path(data_dir,'fleetNames.csv'))
all_fleet
fleet<-filter(all_fleet,order>0)

###########
#
# 1. Calculate the angles of the AMOEBAs.
# The table of yields simulated with status quo effort levels is statusq.xls
# Import this table to Splus and name it "status".  The rows are the fleet abbreviations

jeremy[["STATUS"]]

status<-read_csv(file=file.path(data_dir,'statusq.csv'))  # status quo yield
status<-as.matrix(status)

dimnames(status)[[1]] <- fleet$fleets
round(status)

# Now do the principal components analysis on the correlation matix
status.prc <- princomp(t(status),cor=T)
summary(status.prc)
# Note that the first two principal components explain 55% of the variance
# The fleet angles come from the loadings of the first two principal components
theta <- ang(status.prc$loadings[,1],status.prc$loadings[,2])
# The species angles come from the scores of the first two principal components
omega <- ang(status.prc$scores[,1],status.prc$scores[,2])
# Both theta and omega were jittered slightly so they wouldn't be superimposed.
#


jeremy[["THETA"]]
theta # not the same, but THETHA might be from another task (SSB plot)

jeremy[["OMEGA"]]
omega # not the same

# 2. Calculate the model for spawning stock biomass
# The table of spawning stock biomasses simulated with 10% increase in effort in each
# fleet is ssb110.xls.  The first row of this table is with no effort change (status.quo).
# Import this table to Splus and call it "ssbtab"
ssbtab<-read_csv(file=file.path(data_dir,'ssb110.csv'))
ssbtab<-as.matrix(ssbtab)
rownames(ssbtab)<-all_fleet$fleets
ssbtab
# just checking
jeremy[["SSBTAB"]]

# emat is the (n.fleet+1)*(n.fleet+1) matrix with the corresponding effort combinations
make_effort_matrix<-function(efc) {
  emat<-matrix(nrow=n.fleet+1,ncol=n.fleet+1)
  emat[,]<-1
  diag(emat)<-efc
  emat[1,1]<-1
  return(emat)
}
emat<-make_effort_matrix(efc=eff.change)
# just checking
jeremy[["emat"]]-emat

# einv is its inverse calculated as
 einv <- solve(emat)
 einv
 # just checking
 round(jeremy[["einv"]])
 
 
 ssbparm <- einv%*%ssbtab
 ssbparm
 rownames(ssbparm)<-all_fleet$fleets
 
 # just checking
 round(jeremy[["SSBPARM"]],1)
 round(ssbparm,1)
 round(jeremy[["SSBPARM"]]-ssbparm,3)
 
 
 # just testing
 if (FALSE) {
   eff_chg<-1
   effnew <- c(1, rep(eff_chg,n.fleet)) # effnew is vector of effort multipliers, one per fleet
   ssb <- t(ssbparm) %*% effnew
   ssb- ssbtab[1,] #Fsq should be the same as first row of ssbparam
  
   effnew[n.fleet+1]<-1.1
   ssb <- t(ssbparm) %*% effnew
   ssb- ssbtab[n.fleet+1,] # should be the same as last row of ssbparam
}
# it works !!!
 
 
 
 
### prediction, just using Jeremy's plot 
eff_scale<-0.75
eff <- rep(eff_scale,n.fleet)		# creates a vector of effort levels, one for each fleet
maxlen<-1.5 # MV
ssbplot(eff)            # plots two amoebas, one for effort and one for ssb


#
# 3. Calculate the models for each of the fishing fleets
# The spreadsheet containing the simulations with effort increases by 10% is called yield110.csv

y<-read_csv(file=file.path(data_dir,'yield110.csv'))
filter(y,fleetName=='Fixed gear')
y<-left_join(y,all_fleet, by = c("fleetName"="fleetNames")) %>% mutate(fleetName=paste(formatC(order,width=2),fleetName,sep='_'),fleets=NULL,fletNames=NULL,order=NULL) 
y<-left_join(y,all_fleet, by = c("id"="fleetNames")) %>% mutate(id=paste(formatC(order,width=2),id,sep='_'),fleetNames=NULL,order=NULL,fleets=NULL) 
filter(y,fleetName==" 2_Fixed gear" &species=='COD')

yy<-tapply(y$Yield,list(y$fleetName,y$id,y$species),sum)  # reformat
yy[1,,]

yy[2,,]

#estimates yield parameters for each fleet, 
yield_params<-lapply(1:dim(yy)[[1]],function(fltn){
  cpue<-yy[fltn,,]
  cpue[1 + fltn, ] <- cpue[1 + fltn,]/eff.change
  params <- einv %*% cpue
  rownames(params) <-rownames(cpue)
  return(params)
})


names(yield_params)<-fleet$fleets
yield_params["fix"]
jeremy[["PARAMS"]]$fix

# just an example
eff_chg<-0.78
effnew <- c(1, rep(eff_chg,n.fleet)) # effnew is vector of effort multipliers, one per fleet
names(effnew)<-all_fleet$fleets
effnew

yield <- matrix(nrow = n.species, ncol = n.fleet)
colnames(yield)<-fleet$fleets
rownames(yield)<- colnames(yield_params[[1]])

#catch by one fleet 
my.gear<-'fix'
cpuehat <- t(yield_params[[my.gear]]) %*% effnew
yhat<-cpuehat*effnew[my.gear]
yhat

#catch by all fleets   
for (j in 1:n.fleet) {
  cpuehat <- t(yield_params[[j]]) %*% effnew
  yhat <- cpuehat * effnew[j+1]
  yield[, j ] <- yhat
}
effnew
round(cbind(yield,all=rowSums(yield,na.rm=TRUE)))

##  optimisation
jeremy[["FIND.MSY"]]
find.msy<-function (parm = params) {
  avec <- c(1:n.fleet)
  B <- matrix(rep(avec, n.fleet), nrow = n.fleet)
  for (i in 1:n.fleet) {
    avec[i] <- sum(parm[[i]][1, ],na.rm = T)
  }
  for (i in 1:n.fleet) {
    for (j in 1:n.fleet) {
      B[i, j] <- sum(parm[[i]][j + 1, ],na.rm = T)
      if (i == j) B[i, j] <- B[i, j] * 2
    }
  }
 list(avec=avec, B=B)
}

opt<-find.msy(parm=yield_params)
opt
jeremy[["AVEC"]]
jeremy[["bmat"]]

### MSY effort
Binv <- solve(-opt[["B"]])
Emsy <- Binv %*% opt[["avec"]]

rownames(Emsy)<-fleet$fleets
Emsy
jeremy[["effmsy"]]  # jeremy has constrained fixed [2,] to 0.2, but practically the same values


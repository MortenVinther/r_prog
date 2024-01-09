#
### init: Run: init.R in the SMS directory and _init_rsms.R
source(file.path(rsms.root.prog,"make_rsms_data_function.R"))
source(file.path(rsms.root.prog,"pick_species.R"))
source(file.path(rsms.root.prog,"rsms_function.R"))

### Extract data from SMS

if (FALSE) {  # transform  SMS data into RSMS format 
  inp_all<-make_rsms_data(dir="ns_2023_ss_input",outDir=rsms.root)
  save(inp_all,file=file.path(rsms.root,"rsms_input_all.Rdata"))
  load(file=file.path(rsms.root,"rsms_input_all.Rdata"),verbose=TRUE)
}
load(file=file.path(rsms.root,"rsms_input_all.Rdata"),verbose=TRUE)

# read a data set on SMS format
# inp<-make_rsms_data(dir="S21",outDir=rsms.root)

# select a combination of species from the (full) data set
inp<-pick_species(ps=c(2L), inp=inp_all) 

# inp=inp_all

#  transform quarterly data into to annual data (testing)
if (FALSE) inp<-into_annual(inp)



inp$data$spNames

#### prepare to run

data<-inp[['data']]
parameters<-inp[['parameters']]

data$stockRecruitmentModelCode[]<-1L
data$Debug<-1L


#### Run rsms


#Going back to our objective function f, first step is to check that you can evaluate the function as a normal R function:
# func(parameters)   # KALDET VIL PÅVIRKE KALDET TIL MakeAdFun !!!?
#An error at this point is obviously not due to RTMB.

# adjust if the are species/year combination with zero catches (or assumed F is very low and highly uncertain)
if (data$zeroCatchYearExists==1) {
  UfMap<-matrix(1L:(dim(parameters$Uf)[[1]]*dim(parameters$Uf)[[2]]),nrow=sum(data$nlogF),byrow=TRUE)
  for (s in 1:data$nSpecies) if (length(data$zeroCatchYear[[s]]) >0 ) {
    zy<-data$zeroCatchYear[[s]]
    fromTo<-data$nlogFfromTo[s,]
    UfMap[fromTo[1]:fromTo[2],zy]<-NA
    parameters$Uf[fromTo[1]:fromTo[2],zy]<-log(0.001)
  }
  UfMap<-factor(UfMap)
}

if (data$zeroCatchYearExists==1) my.map<-list(Uf=UfMap) else my.map=list()

  #logSdLogObsSurvey=factor(rep(NA,length(parameters$logSdLogObsSurvey))),
  # logSdLogN =factor(rep(NA,length(parameters$logSdLogN)))
  #rho =factor(rep(NA,length(parameters$rho)))


obj <- MakeADFun(func, parameters, random=c("Un","Uf"),silent=FALSE,map=my.map)

#obj$simulate()
#checkConsistency(obj)
                 
lower <- obj$par*0-Inf
upper <- obj$par*0+Inf
lower["rho"] <- 0.01
upper["rho"] <- 0.99

nl<-names(lower)

lower[nl=="logSdLogObsSurvey"]<-rep(log(0.15),length(parameters$logSdLogObsSurvey))
upper[nl=="logSdLogObsSurvey"]<-rep(log(2.0),length(parameters$logSdLogObsSurvey))

lower[nl=="logSdLogObsCatch"]<-rep(log(0.1),length(parameters$logSdLogObsCatch))
upper[nl=="logSdLogObsCatch"]<-rep(log(2.0),length(parameters$logSdLogObsCatch))


# N.sandeel
#lower[nl=="logSdLogN"]<-log(c(0.1,0.1))
#upper[nl=="logSdLogN"]<-log(c(20,0.3))

#t(rbind(lower,upper))    
opt <- nlminb(obj$par, obj$fn, obj$gr, lower=lower, upper=upper)

cat("\nobjective:",opt$objective,"  convergence:",opt$convergence) # 0 indicates successful convergence.

